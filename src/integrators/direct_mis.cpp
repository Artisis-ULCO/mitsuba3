#include <mitsuba/render/integrator.h>
#include <mitsuba/render/bsdf.h>
#include <mitsuba/render/emitter.h>
#include <mitsuba/core/properties.h>

NAMESPACE_BEGIN(mitsuba)

/**!

.. _integrator-direct:

Direct illumination integrator (:monosp:`direct`) with custom MIS
-----------------------------------------------------------------

.. pluginparameters::

 * - shading_samples
   - |int|
   - This convenience parameter can be used to set both :code:`emitter_samples` and
     :code:`bsdf_samples` at the same time.

 * - emitter_samples
   - |int|
   - Optional more fine-grained parameter: specifies the number of samples that should be generated
     using the direct illumination strategies implemented by the scene's emitters.
     (Default: set to the value of :monosp:`shading_samples`)

 * - bsdf_samples
   - |int|
   - Optional more fine-grained parameter: specifies the number of samples that should be generated
     using the BSDF sampling strategies implemented by the scene's surfaces.
     (Default: set to the value of :monosp:`shading_samples`)

 * - hide_emitters
   - |bool|
   - Hide directly visible emitters.
     (Default: no, i.e. |false|)

.. subfigstart::
.. subfigure:: ../../resources/data/docs/images/render/integrator_direct_bsdf.jpg
   :caption: (**a**) BSDF sampling only
   :label: fig-direct-bsdf
.. subfigure:: ../../resources/data/docs/images/render/integrator_direct_lum.jpg
   :caption: (**b**) Emitter sampling only
   :label: fig-direct-lum
.. subfigure:: ../../resources/data/docs/images/render/integrator_direct_both.jpg
   :caption: (**c**) MIS between both sampling strategies
   :label: fig-direct-both
.. subfigend::
   :width: 0.32
   :label: fig-direct

This integrator implements a direct illumination technique that makes use
of *multiple importance sampling*: for each pixel sample, the
integrator generates a user-specifiable number of BSDF and emitter
samples and combines them using the power heuristic. Usually, the BSDF
sampling technique works very well on glossy objects but does badly
everywhere else (**a**), while the opposite is true for the emitter sampling
technique (**b**). By combining these approaches, one can obtain a rendering
technique that works well in both cases (**c**).

The number of samples spent on either technique is configurable, hence
it is also possible to turn this plugin into an emitter sampling-only
or BSDF sampling-only integrator.

.. note:: This integrator does not handle participating media or indirect illumination.

.. tabs::
    .. code-tab::  xml
        :name: direct-integrator

        <integrator type="directmis"/>

    .. code-tab:: python

        'type': 'direct'

 */

template <typename Float, typename Spectrum>
class DirectIntegratorMIS : public SamplingIntegrator<Float, Spectrum> {
public:
    MI_IMPORT_BASE(SamplingIntegrator, m_hide_emitters)
    // [MIS]: add of MIS Model type from `mitsuba/render/fwd.h`
    MI_IMPORT_TYPES(Scene, Sampler, Medium, Emitter, EmitterPtr, BSDF, BSDFPtr, MISModel)

    DirectIntegratorMIS(const Properties &props) : Base(props) {
        // if (!props.has_property("batch_samples")) {
        //     Throw("Required property 'batch_samples'");
        // }

        // [MIS] current direct integrator parameters has been removed
        // use only batch samples information
        // batch_samples = props.get<size_t>("batch_samples", 10);
    }

    std::pair<Spectrum, Mask> sample(const Scene *scene,
                                     Sampler *sampler,
                                     const RayDifferential3f &ray,
                                     MISModel *mis,
                                     const Medium * /* medium */,
                                     Float * /* aovs */,
                                     Mask active) const override {
        MI_MASKED_FUNCTION(ProfilerPhase::SamplingIntegratorSample, active);

        // [MIS] TODO: check if necessary to sample individually (different ray)
        // Classical balance heuristic gives better results when sample once 
        SurfaceInteraction3f si = scene->ray_intersect(
            ray, +RayFlags::All, /* coherent = */ true, active);
        Mask valid_ray = active && si.is_valid();

        Spectrum result(0.f);

        // ----------------------- Visible emitters -----------------------
        if (!m_hide_emitters) {
            EmitterPtr emitter_vis = si.emitter(scene, active);
            if (dr::any_or<true>(dr::neq(emitter_vis, nullptr)))
                result += emitter_vis->eval(si, active);
        }

        active &= si.is_valid();
        if (dr::none_or<false>(active))
            return { result, valid_ray };

        // [MIS] compute number of samples and weight (only if ray is valid)
        uint32_t batch_samples = mis->n_batch_samples();
        Float alpha_bsdf = mis->get_alpha(0);
        Float alpha_emitter = mis->get_alpha(1);

        uint32_t m_bsdf_samples = (uint32_t)(batch_samples * alpha_bsdf);
        uint32_t m_emitter_samples = batch_samples - m_bsdf_samples;

        ScalarFloat m_weight_bsdf = 1.f / (ScalarFloat) m_bsdf_samples;
        ScalarFloat m_weight_emitter  = 1.f / (ScalarFloat) m_emitter_samples;

        // ----------------------- Emitter sampling -----------------------
        BSDFContext ctx;
        BSDFPtr bsdf = si.bsdf(ray);
        auto flags = bsdf->flags();
        Mask sample_emitter = active && has_flag(flags, BSDFFlags::Smooth);

        if (dr::any_or<true>(sample_emitter)) {
            for (size_t i = 0; i < m_emitter_samples; ++i) {

                mis->update_n_samples();

                // [MIS] add sample for light sampling method
                mis->update_n_samples_method(1);

                Mask active_e = sample_emitter;
                DirectionSample3f ds;
                Spectrum emitter_val;
                std::tie(ds, emitter_val) = scene->sample_emitter_direction(
                    si, sampler->next_2d(active_e), true, active_e);
                active_e &= dr::neq(ds.pdf, 0.f);
                if (dr::none_or<false>(active_e))
                    continue;

                // Query the BSDF for that emitter-sampled direction
                Vector3f wo = si.to_local(ds.d);

                /* Determine BSDF value and probability of having sampled
                   that same direction using BSDF sampling. */
                auto [bsdf_val, bsdf_pdf] = bsdf->eval_pdf(ctx, si, wo, active_e);
                bsdf_val = si.to_world_mueller(bsdf_val, -wo, si.wi);

                // [MIS]: add light sampling luminance and PDFs
                mis->add_sampling_data(1, bsdf_val * emitter_val, {bsdf_pdf, ds.pdf});

                Float mis_w = dr::select(ds.delta, Float(1.f), mis->mis_weight(alpha_emitter,
                    ds.pdf, bsdf_pdf)) * m_weight_emitter;
                result[active_e] += mis_w * bsdf_val * emitter_val;
            }
        }

        // ------------------------ BSDF sampling -------------------------
        for (size_t i = 0; i < m_bsdf_samples; ++i) {

            mis->update_n_samples();

            // [MIS] add sample for bsdf sampling method
            mis->update_n_samples_method(0);

            auto [bs, bsdf_val] = bsdf->sample(ctx, si, sampler->next_1d(active),
                                               sampler->next_2d(active), active);
            bsdf_val = si.to_world_mueller(bsdf_val, -bs.wo, si.wi);

            Mask active_b = active && dr::any(dr::neq(unpolarized_spectrum(bsdf_val), 0.f));

            // Trace the ray in the sampled direction and intersect against the scene
            SurfaceInteraction3f si_bsdf =
                scene->ray_intersect(si.spawn_ray(si.to_world(bs.wo)), active_b);

            // Retain only rays that hit an emitter
            EmitterPtr emitter = si_bsdf.emitter(scene, active_b);
            active_b &= dr::neq(emitter, nullptr);

            if (dr::any_or<true>(active_b)) {
                Spectrum emitter_val = emitter->eval(si_bsdf, active_b);
                Mask delta = has_flag(bs.sampled_type, BSDFFlags::Delta);

                /* Determine probability of having sampled that same
                   direction using Emitter sampling. */
                DirectionSample3f ds(scene, si_bsdf, si);

                Float emitter_pdf =
                    dr::select(delta, 0.f, scene->pdf_emitter_direction(si, ds, active_b));

                // [MIS]: add bsdf sampling luminance and PDFs
                mis->add_sampling_data(0, bsdf_val * emitter_val, {bs.pdf, emitter_pdf});

                result[active_b] +=
                    bsdf_val * emitter_val *
                    mis->mis_weight(alpha_bsdf, bs.pdf, emitter_pdf) * m_weight_bsdf;
            } else {
                // no emitter found, hence no contribution
                mis->add_sampling_data(0, 0.f, {bs.pdf, 0.f});
            }
        }

        // [MIS]: update alphas and number of samples
        mis->update_n_iterations();

        if (mis->number_of_samples() >= mis->start_at_samples()) {
            // std::cout << " -- start update of alphas" << std::endl;
            mis->update_alphas();
        }

        return { result, valid_ray };
    }

    std::string to_string() const override {
        std::ostringstream oss;
        oss << "MISDirectIntegrator[" << std::endl
            // << "  batch_samples = " << batch_samples << "," << std::endl
            << "]";
        return oss.str();
    }

    MI_DECLARE_CLASS()

private:
    // uint32_t batch_samples; // require even number of samples
};

MI_IMPLEMENT_CLASS_VARIANT(DirectIntegratorMIS, SamplingIntegrator)
MI_EXPORT_PLUGIN(DirectIntegratorMIS, "Direct integrator MIS");
NAMESPACE_END(mitsuba)
