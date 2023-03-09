#include <mitsuba/render/mis.h>

#include <limits>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name MISModel implementations
// =======================================================================

MI_VARIANT MISModel<Float, Spectrum>::MISModel(uint32_t n_methods) 
    : Object(), n_methods(n_methods), n_samples(0) {

    for (uint32_t i = 0; i < n_methods; i++) {

        // init MIS data
        n_samples_methods.insert({i, 0.f});
        luminance_sum.insert({i, 0.f});
        squared_sum.insert({i, 0.f});
    
        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.f);
        pdf_sum.insert({i, methods_pdf_sum});
    }

    // init alphas (depending of MIS method)
    this->init_alphas();
};

MI_VARIANT MISModel<Float, Spectrum>::~MISModel() { }

MI_VARIANT void MISModel<Float, Spectrum>::init_alphas() {
    // [MIS]
    for (uint32_t i = 0; i < this->n_methods; i++) {
        this->alphas.insert({i, 1. / (double)this->n_methods});
    }
}

MI_VARIANT void MISModel<Float, Spectrum>::update_n_samples() {
    n_samples++;
}


MI_VARIANT Float MISModel<Float, Spectrum>::get_alpha(uint32_t sampling_method_id) const {
    return alphas.at(sampling_method_id);
}

MI_VARIANT void MISModel<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {

    // [MIS] TODO: check if required use of scalar_RGB mode only
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();
    luminance_sum[sampling_method_id] += lum;

    // cumulate PDF sum
    std::vector<Float> pdf_method = pdf_sum[sampling_method_id];

    for (uint32_t i = 0; i < n_methods; i++)
        pdf_method[i] += pdfs[i];

    // increase number of sample for this method
    n_samples_methods[sampling_method_id] += 1;

    Float mean_luminance = (luminance_sum[sampling_method_id] / n_samples_methods[sampling_method_id]);
    squared_sum[sampling_method_id] += (mean_luminance - lum) * (mean_luminance - lum);

    // [MIS]: Debug
    // std::cout << "-----------------------------" << std::endl;
    // std::cout << sampling_method_id << " has now " << n_samples_methods[sampling_method_id] << " samples" << std::endl;
    // std::cout << sampling_method_id << " has now " << luminance_sum[sampling_method_id] << " luminance (sum)" << std::endl;
    // std::cout << sampling_method_id << " has now pdfs (sum) : [ ";

    // for (uint32_t i = 0; i < n_methods; i++)
    //     std::cout << pdf_method[i] << " ";
    // std::cout << "]" << std::endl;
}   

//! @}
// =======================================================================

// =======================================================================
//! @{ \name MISBalance implementations
// =======================================================================

MI_VARIANT MISBalance<Float, Spectrum>::MISBalance(uint32_t n_methods) : Base(n_methods) {
}

MI_VARIANT void MISBalance<Float, Spectrum>::update_alphas() {
    // [MIS]: balance do nothing
}

MI_VARIANT Float MISBalance<Float, Spectrum>::mis_weight(Float /* alpha */, Float pdf_a, Float pdf_b) const {
    Float w = pdf_a / (pdf_a + pdf_b);
    return dr::select(dr::isfinite(w), w, 0.f);
}


//! @}
// =======================================================================


// =======================================================================
//! @{ \name MISPower implementations
// =======================================================================

MI_VARIANT MISPower<Float, Spectrum>::MISPower(uint32_t n_methods) : Base(n_methods) {
}

MI_VARIANT void MISPower<Float, Spectrum>::update_alphas() {
    // [MIS]: power do nothing
}

MI_VARIANT Float MISPower<Float, Spectrum>::mis_weight(Float /* alpha */, Float pdf_a, Float pdf_b) const {
    pdf_a *= pdf_a;
    pdf_b *= pdf_b;
    Float w = pdf_a / (pdf_a + pdf_b);
    return dr::select(dr::isfinite(w), w, 0.f);
}

//! @}
// =======================================================================

// =======================================================================
//! @{ \name MISDivergence implementations
// =======================================================================

MI_VARIANT MISDivergence<Float, Spectrum>::MISDivergence(uint32_t n_methods) : Base(n_methods) {
}

MI_VARIANT Float MISDivergence<Float, Spectrum>::mis_weight(Float alpha, Float pdf_a, Float pdf_b) const {
    // [MIS]: Fixed to two sampling methods only
    Float w = (pdf_a * alpha) / (pdf_a * alpha + pdf_b * (1.f - alpha));
    return dr::select(dr::isfinite(w), w, 0.f);
}

//! @}
// =======================================================================

// =======================================================================
//! @{ \name MISLinear1 implementations
// =======================================================================

MI_VARIANT MISLinear1<Float, Spectrum>::MISLinear1(uint32_t n_methods) : Base(n_methods) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISLinear1<Float, Spectrum>::update_alphas() {
    // [MIS]: update alphas using first linear heuristic
    Float eps = std::numeric_limits<Float>::epsilon();
    Float f1 = (luminance_sum[0] / (n_samples_methods[0] + eps));
    Float f2 = (luminance_sum[1] / (n_samples_methods[1] + eps));
            
    Float p11 = (pdf_sum[0][0] / (n_samples_methods[0] + eps));
    Float p21 = (pdf_sum[0][1] / (n_samples_methods[0] + eps));
    Float p12 = (pdf_sum[1][0] / (n_samples_methods[1] + eps));
    Float p22 = (pdf_sum[1][1] / (n_samples_methods[1] + eps));

    // Float f1 = luminance_sum[0];
    // Float f2 = luminance_sum[1];
            
    // Float p11 = pdf_sum[0][0];
    // Float p21 = pdf_sum[0][1];
    // Float p12 = pdf_sum[1][0];
    // Float p22 = pdf_sum[1][1];

    Float nominator = p22 * f1 - p21 * f2;
    Float denominator = p11 * f2 - p21 * f2 - p12 * f1 + p22 * f1;

    // alpha[0] is for Emitter sampling 
    // alpha[1] is for BSDF sampling
    alphas[0] = nominator / (denominator + eps);

    // compute variances
    std::vector<Float> variances;

    for (uint32_t i = 0; i < this->n_methods; i++) {
        variances.push_back(squared_sum[i] / (n_samples_methods[i] + eps));
        // std::cout << "Method n°" << i << " has luminance variance of " << variances[i] << std::endl;
    }

    // if out of expected bounds then check contribution variance
    if (alphas[0] <= 0 || alphas[0] >= 1)
        alphas[0] = variances[0] >= variances[1] ? 0.99f : 0.01f;
        
    alphas[1] = 1.f - alphas[0];

    // [MIS] Debug
    // std::cout << "[Update] MIS has now alphas : [ ";
    // for (uint32_t i = 0; i < this->n_methods; i++)
    //     std::cout << alphas[i] << " ";
    // std::cout << "]" << std::endl;
}

//! @}
// =======================================================================


MI_IMPLEMENT_CLASS_VARIANT(MISModel, Object, "MIS")
MI_IMPLEMENT_CLASS_VARIANT(MISBalance, MISModel, "MIS balance")
MI_IMPLEMENT_CLASS_VARIANT(MISPower, MISModel, "MIS power")
MI_IMPLEMENT_CLASS_VARIANT(MISDivergence, MISModel, "MIS Divergence")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear1, MISDivergence, "MIS Linear 1")

MI_INSTANTIATE_CLASS(MISModel)
MI_INSTANTIATE_CLASS(MISBalance)
MI_INSTANTIATE_CLASS(MISPower)
MI_INSTANTIATE_CLASS(MISDivergence)
MI_INSTANTIATE_CLASS(MISLinear1)

NAMESPACE_END(mitsuba)