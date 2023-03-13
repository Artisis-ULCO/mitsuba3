#include <mitsuba/render/mis.h>

#include <limits>
#include <cmath>

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

MI_VARIANT uint32_t MISModel<Float, Spectrum>::number_of_methods() const {
    return n_methods;
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

// =======================================================================
//! @{ \name MISLinear2 implementations
// =======================================================================

MI_VARIANT MISLinear2<Float, Spectrum>::MISLinear2(uint32_t n_methods) : Base(n_methods) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISLinear2<Float, Spectrum>::update_alphas() {
    // [MIS]: update alphas using first linear heuristic
    Float eps = std::numeric_limits<Float>::epsilon();

    Float n1 = n_samples_methods[0];
    Float n2 = n_samples_methods[1];
    Float f1 = luminance_sum[0];
    Float f2 = luminance_sum[1];
            
    Float p11 = pdf_sum[0][0];
    Float p21 = pdf_sum[0][1];
    Float p12 = pdf_sum[1][0];
    Float p22 = pdf_sum[1][1];

    Float nominator = (n1 * f1 * p22) - (n2 * f2 * p21);
    Float denominator = (n2 * f2 * p11) - (n1 * f1 * p12) - (n2 * f2 * p21) + (n1 * f1 * p22);

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

// =======================================================================
//! @{ \name MISLinear3 implementations
// =======================================================================

MI_VARIANT MISLinear3<Float, Spectrum>::MISLinear3(uint32_t n_methods) : Base(n_methods) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 

    for (uint32_t i = 0; i < n_methods; i++) {
        sum_mu[i] = 0.f;
    }
}

MI_VARIANT void MISLinear3<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {
    
    // pdfs[0] is BSDF PDF
    // pdfs[1] is Light PDF
    Base::add_sampling_data(sampling_method_id, luminance, pdfs);

    // TODO: improve (recomputed lum: from Base)
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();

    Float pdf_a, pdf_b = 0.f;

    // BSDF sampling
    if (sampling_method_id == 0) {
        pdf_a = pdfs[0];
        pdf_b = pdfs[1];
    } else { // Light sampling
        pdf_a = pdfs[1];
        pdf_b = pdfs[0];
    }
    Float w = pdf_a / (pdf_a + pdf_b);
    w = dr::select(dr::isfinite(w), w, 0.f);
    sum_mu[sampling_method_id] += lum * w;
}

MI_VARIANT void MISLinear3<Float, Spectrum>::update_alphas() {
    // [MIS]: update alphas using first linear heuristic
    Float eps = std::numeric_limits<Float>::epsilon();

    Float n1 = n_samples_methods[0];
    Float n2 = n_samples_methods[1];
    Float f1 = luminance_sum[0];
    Float f2 = luminance_sum[1];
            
    Float p11 = pdf_sum[0][0];
    Float p21 = pdf_sum[0][1];
    Float p12 = pdf_sum[1][0];
    Float p22 = pdf_sum[1][1];

    // std::cout << " @sample " << n_samples << std::endl;
    // std::cout << " -- N1 (BSDF): " << n1 << std::endl;
    // std::cout << " -- N2 (light): " << n2 << std::endl;
    // std::cout << " -- \\mu (BSDF): " << sum_mu[0] << std::endl;
    // std::cout << " -- \\mu (light): " << sum_mu[1] << std::endl;

    Float commonDen = n2 * p11 - n2 * p21 - n1 * p12 + n1 * p22;

    // TODO: check sum of number of samples from each method?
    Float mu = (sum_mu[0] + sum_mu[1]) / (Float)(n_samples + eps);

    // std::cout << " -- \\mu: " << mu << std::endl;
    // std::cout << " -- \\den: " << commonDen << std::endl;

    Float firstTerm = (n2 * f1 - n1 * f2) / (mu * commonDen + eps);
    Float secondTerm = (n2 * p21 + n1 * p22) / (commonDen + eps);

    // alpha[0] is for Emitter sampling 
    // alpha[1] is for BSDF sampling
    alphas[0] = firstTerm - secondTerm;

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



// =======================================================================
//! @{ \name MISTsallis implementations
// =======================================================================

MI_VARIANT MISTsallis<Float, Spectrum>::MISTsallis(uint32_t n_methods, Float gamma, uint32_t batch_samples) 
    : Base(n_methods), gamma(gamma), batch_samples(batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT Float MISTsallis<Float, Spectrum>::get_alpha(uint32_t sampling_method_id) const {
    
    // [MIS] 5 samples for first method / 5 samples for second method
    if (sampling_method_id == 0) {
        if (n_samples < (uint32_t)(batch_samples / 2)) {
            return 0.99;
        } else if (n_samples < batch_samples) {
            return 0.01;
        }
    } else {
        // reverse if sampling method 1
        if (n_samples < (uint32_t)(batch_samples / 2)) {
            return 0.01;
        } else if (n_samples < batch_samples) {
            return 0.99;
        }
    }

    // after we always get the computed MIS alpha
    return alphas.at(sampling_method_id);
}

MI_VARIANT void MISTsallis<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {

    // Classical track
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();
    luminance_sum[sampling_method_id] += lum;

    // increase number of sample for this method
    n_samples_methods[sampling_method_id] += 1;

    Float mean_luminance = (luminance_sum[sampling_method_id] / n_samples_methods[sampling_method_id]);
    squared_sum[sampling_method_id] += (mean_luminance - lum) * (mean_luminance - lum);

    // Tsallis track
    // pdfs[0] is BSDF PDF
    // pdfs[1] is Light PDF

    // We expect:
    // p_1 is BSDF PDF
    // p_2 is Light PDF
    Float f_pdf = pdfs.at(0);
    Float g_pdf = pdfs.at(1);

    // Compute xi_sum and xi_prime_sum
    Float lum_tsallis = std::pow(lum, gamma);

    // [MIS]: update xi and xi' with respect to equations 30-31
    Float alpha_probs = alphas.at(0) * f_pdf + alphas.at(1) * g_pdf;

    // check not null value
    if (alpha_probs <= 0)
        alpha_probs += std::numeric_limits<Float>::epsilon();

    xi_sum += (lum_tsallis / std::pow(alpha_probs, gamma + 1)) 
                    * (f_pdf - g_pdf);

    xi_prime_sum += (lum_tsallis / std::pow(alpha_probs, gamma + 2)) 
                    * (std::pow(f_pdf - g_pdf, 2));
}   

MI_VARIANT void MISTsallis<Float, Spectrum>::update_alphas() {
    
    // TODO: compute new alphas (depending of n_samples and batch)
    if ((n_samples % batch_samples) != 0)
        return;

    Float eps = std::numeric_limits<Float>::epsilon();

    // Need to update UpdateProbsMIS, in order to take into account 
    // luminance without balance heuristic and then update internal data 
    // in order to compute xi_{\alpha}
    Float xi_alpha = xi_sum / batch_samples;

    Float xi_prime_alpha = xi_prime_sum * (-gamma / batch_samples);

    if (xi_prime_alpha <= 0)
        xi_prime_alpha += eps;

    // std::cout << " -- Current alpha: " << alphas[0] << std::endl;
    // std::cout << " -- xi_alpha: " << xi_alpha << std::endl;
    // std::cout << " -- xi_prime_alpha: " << xi_prime_alpha << std::endl;
    // std::cout << " -- Gradient step: " << (xi_alpha / xi_prime_alpha) << std::endl;
    // std::cout << " -- new Alpha MIS: " << alphas[0] - (xi_alpha / xi_prime_alpha) << std::endl;

    // [MIS]: update xi and xi' with respect to equation 32
    alphas[0] = alphas[0] - (xi_alpha / xi_prime_alpha);
    // std::cout << "[spp : " << n_samples << " ]:: computed alpha: " << alphas[0] << std::endl;

    // reset for next batch (number of generated paths)
    reset();

    // compute variances
    std::vector<Float> variances;

    for (uint32_t i = 0; i < this->n_methods; i++)
        variances.push_back(squared_sum[i] / (n_samples_methods[i] + eps));

    // if out of expected bounds then check contribution variance
    if (alphas[0] <= 0 || alphas[0] >= 1)
        alphas[0] = variances[0] >= variances[1] ? 0.99f : 0.01f;
        
    // std::cout << " -- Real new Alpha MIS: " << alphas[0] << std::endl;
    alphas[1] = 1.f - alphas[0];
}

MI_VARIANT void MISTsallis<Float, Spectrum>::reset() {

    // reset tracked data for new samples batch
    xi_sum = 0.f;
    xi_prime_sum = 0.f;
}

//! @}
// =======================================================================

MI_IMPLEMENT_CLASS_VARIANT(MISModel, Object, "MIS")
MI_IMPLEMENT_CLASS_VARIANT(MISBalance, MISModel, "MIS balance")
MI_IMPLEMENT_CLASS_VARIANT(MISPower, MISModel, "MIS power")
MI_IMPLEMENT_CLASS_VARIANT(MISDivergence, MISModel, "MIS Divergence")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear1, MISDivergence, "MIS Linear 1")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear2, MISDivergence, "MIS Linear 2")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear3, MISDivergence, "MIS Linear 3")
MI_IMPLEMENT_CLASS_VARIANT(MISTsallis, MISDivergence, "MIS Tsallis")

MI_INSTANTIATE_CLASS(MISModel)
MI_INSTANTIATE_CLASS(MISBalance)
MI_INSTANTIATE_CLASS(MISPower)
MI_INSTANTIATE_CLASS(MISDivergence)
MI_INSTANTIATE_CLASS(MISLinear1)
MI_INSTANTIATE_CLASS(MISLinear2)
MI_INSTANTIATE_CLASS(MISLinear3)
MI_INSTANTIATE_CLASS(MISTsallis)

NAMESPACE_END(mitsuba)