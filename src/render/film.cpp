#include <mitsuba/render/film.h>
#include <mitsuba/core/plugin.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/render/mis.h>

NAMESPACE_BEGIN(mitsuba)

MI_VARIANT Film<Float, Spectrum>::Film(const Properties &props) : Object() {
    bool is_m_film = string::to_lower(props.plugin_name()) == "mfilm";

    // Horizontal and vertical film resolution in pixels
    m_size = ScalarVector2u(
        props.get<uint32_t>("width", is_m_film ? 1 : 768),
        props.get<uint32_t>("height", is_m_film ? 1 : 576)
    );

    // Crop window specified in pixels - by default, this matches the full sensor area.
    ScalarVector2u crop_size = ScalarVector2u(
        props.get<uint32_t>("crop_width", m_size.x()),
        props.get<uint32_t>("crop_height", m_size.y())
    );

    ScalarPoint2u crop_offset = ScalarPoint2u(
        props.get<uint32_t>("crop_offset_x", 0),
        props.get<uint32_t>("crop_offset_y", 0)
    );

    set_crop_window(crop_offset, crop_size);

    /* If set to true, regions slightly outside of the film plane will also be
       sampled, which improves the image quality at the edges especially with
       large reconstruction filters. */
    m_sample_border = props.get<bool>("sample_border", false);

    // Use the provided reconstruction filter, if any.
    for (auto &[name, obj] : props.objects(false)) {
        auto *rfilter = dynamic_cast<ReconstructionFilter *>(obj.get());
        if (rfilter) {
            if (m_filter)
                Throw("A film can only have one reconstruction filter.");

            m_filter = rfilter;
            props.mark_queried(name);
        }
    }

    // Use a Gaussian reconstruction filter if none has been specified
    if (!m_filter)
        m_filter =
            PluginManager::instance()->create_object<ReconstructionFilter>(
                Properties("gaussian"));

    // [MIS] get expected type
    std::string mis_model_type = props.get<std::string>("mis", "power");

    // [MIS] initialize the number of sampling data
    for(uint32_t i = 0; i < crop_size.x() * crop_size.y(); i++) {

        // [MIS] TODO specify the sampling method (factory: use of film property?)
        // Fixed n_methods currently
        std::unique_ptr<MISModel> mis_div;
        
        if (mis_model_type == "balance")
            mis_div = std::make_unique<MISBalance<Float, Spectrum>>(2);
        else if (mis_model_type == "power")
            mis_div = std::make_unique<MISPower<Float, Spectrum>>(2);
        else if (mis_model_type == "linear1")
            mis_div = std::make_unique<MISLinear1<Float, Spectrum>>(2);
        else if (mis_model_type == "linear2")
            mis_div = std::make_unique<MISLinear2<Float, Spectrum>>(2);
        else if (mis_model_type == "linear3")
            mis_div = std::make_unique<MISLinear3<Float, Spectrum>>(2);
        else if (mis_model_type == "tsallis") {
            Float gamma = props.get<Float>("gamma", 1.f);
            uint32_t batch_samples = props.get<uint32_t>("batch", 10);
            mis_div = std::make_unique<MISTsallis<Float, Spectrum>>(2, gamma, batch_samples);
        }
        else
            mis_div = std::make_unique<MISPower<Float, Spectrum>>(2);

        mis_models.push_back(std::move(mis_div));
    }
}

MI_VARIANT Film<Float, Spectrum>::~Film() { }

MI_VARIANT void Film<Float, Spectrum>::traverse(TraversalCallback *callback) {
    callback->put_parameter("size", m_size, +ParamFlags::NonDifferentiable);
    callback->put_parameter("crop_size", m_crop_size, +ParamFlags::NonDifferentiable);
    callback->put_parameter("crop_offset", m_crop_offset, +ParamFlags::NonDifferentiable);
}

MI_VARIANT void Film<Float, Spectrum>::parameters_changed(const std::vector<std::string> &keys) {
    ScalarVector2u crop_size = m_crop_size;
    ScalarPoint2u crop_offset = m_crop_offset;

    // Reset the crop window to match the full sensor area if necessary
    if (string::contains(keys, "size")) {
        if (!string::contains(keys, "crop_size"))
            crop_size = ScalarPoint2u(m_size);
        if (!string::contains(keys, "crop_offset"))
            crop_offset = 0;
    }

    set_crop_window(crop_offset, crop_size);
}

MI_VARIANT void
Film<Float, Spectrum>::prepare_sample(const UnpolarizedSpectrum & /* spec */,
                                      const Wavelength & /* wavelengths */,
                                      Float * /* aovs */,
                                      Float /* weight */,
                                      Float /* alpha */,
                                      Mask /* active */) const {
    NotImplementedError("prepare_sample");
}

MI_VARIANT const typename Film<Float, Spectrum>::Texture *
Film<Float, Spectrum>::sensor_response_function() {
    return m_srf.get();
}

MI_VARIANT void Film<Float, Spectrum>::set_crop_window(const ScalarPoint2u &crop_offset,
                                                        const ScalarVector2u &crop_size) {
    if (dr::any(crop_offset + crop_size > m_size))
        Throw("Invalid crop window specification: crop_offset(%u, %u) + "
              "crop_size(%u, %u) > size(%u, %u)", crop_offset.x(), crop_offset.y(),
              crop_size.x(), crop_size.y(), m_size.x(), m_size.y());

    m_crop_size   = crop_size;
    m_crop_offset = crop_offset;
}

MI_VARIANT void Film<Float, Spectrum>::set_size(const ScalarPoint2u &size) {
    m_size = size;
    // Reset the crop window to match the full sensor area
    set_crop_window(ScalarVector2u(0, 0), m_size);
}

MI_VARIANT std::string Film<Float, Spectrum>::to_string() const {
    std::ostringstream oss;
    oss << "Film[" << std::endl
        << "  size = "        << m_size        << "," << std::endl
        << "  crop_size = "   << m_crop_size   << "," << std::endl
        << "  crop_offset = " << m_crop_offset << "," << std::endl
        << "  sample_border = " << m_sample_border << "," << std::endl
        << "  m_filter = " << m_filter << std::endl
        << "]";
    return oss.str();
}


MI_IMPLEMENT_CLASS_VARIANT(Film, Object, "film")
MI_INSTANTIATE_CLASS(Film)
NAMESPACE_END(mitsuba)
