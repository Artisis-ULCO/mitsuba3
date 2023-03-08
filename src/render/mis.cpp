#include <mitsuba/render/mis.h>

NAMESPACE_BEGIN(mitsuba)

// -----------------------------------------------------------------------------

MI_VARIANT MISModel<Float, Spectrum>::MISModel(uint32_t n_methods) 
    : n_methods(n_methods), n_samples(0) {

    for (uint32_t i = 0; i < n_methods; i++) {

        // init MIS data
        n_samples_methods.insert({i, 0});
        luminance_sum.insert({i, 0.});

        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.);
        pdf_sum.insert({i, methods_pdf_sum});
    }

    // init alphas (depending of MIS method)
    this->init_alphas();
};

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

    // TODO: check if required use of scalar_RGB mode only
    Float lum = 0.2126 * luminance[0] + 0.7152 * luminance[1] + 0.0722 * luminance[2];
    luminance_sum[sampling_method_id] += lum;

    // cumulate PDF sum
    std::vector<Float> pdf_method = pdf_sum[sampling_method_id];
    // std::transform(pdf_method.begin(), pdf_method.end(), 
    //     pdfs.cbegin(), pdfs.cend(), [](Float a, const Float b) { return a + b;});

    // increase number of sample for this method
    n_samples_methods[sampling_method_id] += 1;
}

// -----------------------------------------------------------------------------

// Specific Tsallis type
MI_VARIANT MISTsallis<Float, Spectrum>::MISTsallis(uint32_t n_methods) : Base(n_methods) {
}

MI_VARIANT void MISTsallis<Float, Spectrum>::update_alphas() {
    // TODO: [MIS]
}

MI_VARIANT MISTsallis<Float, Spectrum>::~MISTsallis() { }

// -----------------------------------------------------------------------------

MI_IMPLEMENT_CLASS_VARIANT(MISModel, Object, "mis")
MI_IMPLEMENT_CLASS_VARIANT(MISTsallis, MISModel)

MI_INSTANTIATE_CLASS(MISModel)
MI_INSTANTIATE_CLASS(MISTsallis)

NAMESPACE_END(mitsuba)