#include <mutex>

#include <drjit/morton.h>
#include <mitsuba/render/mis.h>

#include <limits>
#include <cmath>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name MISModel implementations
// =======================================================================

MI_VARIANT MISModel<Float, Spectrum>::MISModel(uint32_t n_methods, uint32_t batch_samples) 
    : Object(), n_methods(n_methods), batch_samples(batch_samples), n_samples(0) {

    Float eps = std::numeric_limits<Float>::epsilon();

    for (uint32_t i = 0; i < n_methods; i++) {

        // init MIS data
        n_samples_methods.insert({i, 0.f});
        luminance_sum.insert({i, 0.f});
        squared_sum.insert({i, 0.f});
        max_pdf.insert({i, eps});
    
        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.f);
        pdf_sum.insert({i, methods_pdf_sum});
    }

    // init alphas (depending of MIS method)
    this->init_alphas();
};

MI_VARIANT MISModel<Float, Spectrum>::~MISModel() { }

MI_VARIANT void MISModel<Float, Spectrum>::init_alphas() {
    // [MIS] by default half
    for (uint32_t i = 0; i < this->n_methods; i++) {
        this->alphas.insert({i, 1. / (double)this->n_methods});
    }
}

MI_VARIANT void MISModel<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {

    // [MIS] TODO: check if required use of scalar_RGB mode only
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();
    luminance_sum[sampling_method_id] += lum;

    // cumulate PDF sum
    for (uint32_t i = 0; i < n_methods; i++) {
        if (pdfs[i] > max_pdf[i])
            max_pdf[i] = pdfs[i];

        pdf_sum[sampling_method_id][i] +=  pdfs[i];
    }

    // increase number of sample for this method
    n_samples_methods[sampling_method_id] += 1;

    Float mean_luminance = (luminance_sum[sampling_method_id] / n_samples_methods[sampling_method_id]);
    squared_sum[sampling_method_id] += (mean_luminance - lum) * (mean_luminance - lum);

    // [MIS]: Debug
    // std::cout << "-----------------------------" << std::endl;
    // std::cout << " -- " << sampling_method_id << " has now " << n_samples_methods[sampling_method_id] << " samples" << std::endl;
    // std::cout << " -- Added PDFs are: [";
    // for (uint32_t i = 0; i < n_methods; i++)
    //     std::cout << pdfs[i] << " ";
    // std::cout << "]" << std::endl;

    // std::cout << " -- " << sampling_method_id << " has now " << luminance_sum[sampling_method_id] << " luminance (sum)" << std::endl;
    // std::cout << " -- " << sampling_method_id << " has now pdfs (sum) : [ ";

    // for (uint32_t i = 0; i < n_methods; i++)
    //     std::cout << pdf_sum[i];
    // std::cout << "]" << std::endl;
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

MI_VARIANT uint32_t MISModel<Float, Spectrum>::number_of_samples() const { 
    return n_samples; 
}

MI_VARIANT uint32_t MISModel<Float, Spectrum>::n_batch_samples() const { 
    return batch_samples; 
}

MI_VARIANT uint32_t MISModel<Float, Spectrum>::number_of_samples_method(uint32_t sampling_method_id) const { 
    return n_samples_methods.at(sampling_method_id); 
}

MI_VARIANT std::vector<Float> MISModel<Float, Spectrum>::get_pdfs(uint32_t sampling_method_id) const {
    return pdf_sum.at(sampling_method_id);
}

//! @}
// =======================================================================

// =======================================================================
//! @{ \name MISBalance implementations
// =======================================================================

MI_VARIANT MISBalance<Float, Spectrum>::MISBalance(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
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

MI_VARIANT MISPower<Float, Spectrum>::MISPower(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
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

MI_VARIANT MISDivergence<Float, Spectrum>::MISDivergence(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
}

MI_VARIANT Float MISDivergence<Float, Spectrum>::mis_weight(Float alpha, Float pdf_a, Float pdf_b) const {
    // [MIS]: Fixed to two sampling methods only
    Float w = (pdf_a * alpha) / (pdf_a * alpha + pdf_b * (1.f - alpha));
    return dr::select(dr::isfinite(w), w, 0.f);
}

//! @}
// =======================================================================

// =======================================================================
//! @{ \name MISLight implementations
// =======================================================================

MI_VARIANT MISLight<Float, Spectrum>::MISLight(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISLight<Float, Spectrum>::update_alphas() {
    
    // Importance sampling using Light
    alphas[0] = 0.f;
    alphas[1] = 1.f;
}

MI_VARIANT Float MISLight<Float, Spectrum>::mis_weight(Float alpha, Float pdf_a, Float pdf_b) const {
    // [MIS]: Fixed to two sampling methods only
    Float w = (pdf_a * alpha) / (pdf_a * alpha + pdf_b * (1.f - alpha));
    return dr::select(dr::isfinite(w), w, 0.f);
}

// =======================================================================
//! @{ \name MISBSDF implementations
// =======================================================================

MI_VARIANT MISBSDF<Float, Spectrum>::MISBSDF(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISBSDF<Float, Spectrum>::update_alphas() {
    
    // Importance sampling using BSDF
    alphas[0] = 1.f;
    alphas[1] = 0.f;
}

MI_VARIANT Float MISBSDF<Float, Spectrum>::mis_weight(Float alpha, Float pdf_a, Float pdf_b) const {
    // [MIS]: Fixed to two sampling methods only
    Float w = (pdf_a * alpha) / (pdf_a * alpha + pdf_b * (1.f - alpha));
    return dr::select(dr::isfinite(w), w, 0.f);
}

// =======================================================================
//! @{ \name MISLinear1 implementations
// =======================================================================

MI_VARIANT MISLinear1<Float, Spectrum>::MISLinear1(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISLinear1<Float, Spectrum>::reset() {

    // reset tracked data for new samples batch
    for (uint32_t i = 0; i < n_methods; i++) {

        n_samples_methods[i] = 0.f;
        luminance_sum[i] = 0.f;
        squared_sum[i] = 0.f;
    
        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.f);
        pdf_sum[i] = methods_pdf_sum;
    }
}

MI_VARIANT void MISLinear1<Float, Spectrum>::update_alphas() {

    // TODO: compute new alphas (depending of n_samples and batch)
    // if ((n_samples % batch_samples) != 0)
    //     return;

    // [MIS]: update alphas using first linear heuristic
    Float eps = std::numeric_limits<Float>::epsilon();
    // Float f1 = luminance_sum[0] / (n_samples_methods[0] + eps);
    // Float f2 = luminance_sum[1] / (n_samples_methods[1] + eps);
            
    // Float p11 = pdf_sum[0][0] / (n_samples_methods[0] + eps);
    // Float p21 = pdf_sum[0][1] / (n_samples_methods[0] + eps);
    // Float p12 = pdf_sum[1][0] / (n_samples_methods[1] + eps);
    // Float p22 = pdf_sum[1][1] / (n_samples_methods[1] + eps);

    Float f1 = luminance_sum[0];
    Float f2 = luminance_sum[1];
            
    Float p11 = pdf_sum[0][0];
    Float p21 = pdf_sum[0][1];
    Float p12 = pdf_sum[1][0];
    Float p22 = pdf_sum[1][1];

    Float nominator = p22 * f1 - p21 * f2;
    Float denominator = p11 * f2 - p21 * f2 - p12 * f1 + p22 * f1;

    // std::cout << " @sample " << n_samples << std::endl;
    // std::cout << " -- f1: " << f1 << std::endl;
    // std::cout << " -- f2: " << f2 << std::endl;

    // std::cout << " -- p11: " << p11 << std::endl;
    // std::cout << " -- p21: " << p21 << std::endl;
    // std::cout << " -- p12: " << p12 << std::endl;
    // std::cout << " -- p22: " << p22 << std::endl;

    // std::cout << " -- nominator: " << nominator << std::endl;
    // std::cout << " -- denominator: " << denominator << std::endl;

    // std::cout << " -- previous alpha: " << alphas[0] << std::endl;
    // std::cout << " -- alpha: " << nominator / (denominator + eps) << std::endl;

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
    auto check_alpha = (alphas.at(0) <= 0) or (alphas.at(0) >= 1);
    if (dr::any_or<true>(check_alpha))
        alphas[0] = dr::any_or<true>(variances.at(0) >= variances.at(1)) ? 0.9f : 0.1f;
        
    alphas[1] = 1.f - alphas[0];

    // reset();
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

MI_VARIANT MISLinear2<Float, Spectrum>::MISLinear2(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISLinear2<Float, Spectrum>::reset() {

    // reset tracked data for new samples batch
    for (uint32_t i = 0; i < n_methods; i++) {

        n_samples_methods[i] = 0.f;
        luminance_sum[i] = 0.f;
        squared_sum[i] = 0.f;
    
        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.f);
        pdf_sum[i] = methods_pdf_sum;
    }
}

MI_VARIANT void MISLinear2<Float, Spectrum>::update_alphas() {

    // TODO: compute new alphas (depending of n_samples and batch)
    // if ((n_samples % batch_samples) != 0)
    //     return;

    // [MIS]: update alphas using second linear heuristic
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

    // MIS [Debug]
    // std::cout << " @sample " << n_samples << std::endl;
    // std::cout << " -- N1 (BSDF): " << n1 << std::endl;
    // std::cout << " -- N2 (light): " << n2 << std::endl;
    // std::cout << " -- f1: " << f1 << std::endl;
    // std::cout << " -- f2: " << f2 << std::endl;

    // std::cout << " -- p11: " << p11 << std::endl;
    // std::cout << " -- p21: " << p21 << std::endl;
    // std::cout << " -- p12: " << p12 << std::endl;
    // std::cout << " -- p22: " << p22 << std::endl;

    // std::cout << " -- nominator: " << nominator << std::endl;
    // std::cout << " -- denominator: " << denominator << std::endl;

    // std::cout << " -- previous alpha: " << alphas[0] << std::endl;
    // std::cout << " -- alpha: " << nominator / (denominator + eps) << std::endl;

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
    auto check_alpha = (alphas.at(0) <= 0) or (alphas.at(0) >= 1);
    if (dr::any_or<true>(check_alpha))
        alphas[0] = dr::any_or<true>(variances.at(0) >= variances.at(1)) ? 0.9f : 0.1f;

    alphas[1] = 1.f - alphas[0];

    reset();
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

MI_VARIANT MISLinear3<Float, Spectrum>::MISLinear3(uint32_t n_methods, uint32_t batch_samples) : Base(n_methods, batch_samples) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 

    for (uint32_t i = 0; i < n_methods; i++) {
        sum_mu[i] = 0.f;
    }
}

MI_VARIANT void MISLinear3<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {
    
    // pdfs[0] is Light PDF
    // pdfs[1] is BSDF PDF
    Base::add_sampling_data(sampling_method_id, luminance, pdfs);

    // TODO: improve (recomputed lum: from Base)
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();

    Float pdf_a, pdf_b = 0.f;

    // use of number of samples per method in order to balance
    // and compute \mu
    Float alpha_bsdf = n_samples_methods[0] / ((Float)n_samples_methods[0] + n_samples_methods[1]);
    Float alpha_lum = n_samples_methods[1] / ((Float)n_samples_methods[0] + n_samples_methods[1]);

    // Light sampling
    if (sampling_method_id == 0) {
        pdf_a = pdfs[0] * alpha_bsdf;
        pdf_b = pdfs[1] * alpha_lum;
    } else { // BSDF sampling
        pdf_a = pdfs[1] * alpha_lum;
        pdf_b = pdfs[0] * alpha_bsdf;
    }
    Float w = pdf_a / (pdf_a + pdf_b);
    w = dr::select(dr::isfinite(w), w, 0.f);
    sum_mu[sampling_method_id] += lum * w;
}

MI_VARIANT void MISLinear3<Float, Spectrum>::reset() {

    // reset tracked data for new samples batch
    for (uint32_t i = 0; i < n_methods; i++) {

        n_samples_methods[i] = 0.f;
        luminance_sum[i] = 0.f;
        squared_sum[i] = 0.f;
    
        auto methods_pdf_sum = std::vector<Float>(n_methods, 0.f);
        pdf_sum[i] = methods_pdf_sum;
        sum_mu[i] = 0;
    }
}

MI_VARIANT void MISLinear3<Float, Spectrum>::update_alphas() {

    // TODO: compute new alphas (depending of n_samples and batch)
    // if ((n_samples % batch_samples) != 0)
    //     return;

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
    Float common_den = n2 * p11 - n2 * p21 - n1 * p12 + n1 * p22;

    uint32_t total_samples = 0;
    for (uint32_t i = 0; i < this->n_methods; i++)
        total_samples += n_samples_methods[i];

    // TODO: check sum of number of samples from each method? (or half?)
    Float mu = (sum_mu[0] + sum_mu[1]) / (Float)(total_samples + eps);

    Float first_term = (n2 * f1 - n1 * f2) / (mu * common_den + eps);
    Float second_term = (n2 * p21 + n1 * p22) / (common_den + eps);

    // MIS [Debug]
    // std::cout << " @sample " << n_samples << std::endl;
    // std::cout << " -- N1 (BSDF): " << n1 << std::endl;
    // std::cout << " -- N2 (light): " << n2 << std::endl;
    // std::cout << " -- \\mu (BSDF): " << sum_mu[0] << std::endl;
    // std::cout << " -- \\mu (light): " << sum_mu[1] << std::endl;

    // std::cout << " -- \\mu: " << mu << std::endl;
    // std::cout << " -- \\den: " << common_den << std::endl;

    // std::cout << " -- First term: " << first_term << std::endl;
    // std::cout << " -- Second term: " << second_term << std::endl;
    // std::cout << " -- Previous alpha: " << alphas[0] << std::endl;
    // std::cout << " -- alpha: " << first_term - second_term << std::endl;

    // alpha[0] is for Emitter sampling 
    // alpha[1] is for BSDF sampling
    alphas[0] = first_term - second_term;

    // compute variances
    std::vector<Float> variances;

    for (uint32_t i = 0; i < this->n_methods; i++) {
        variances.push_back(squared_sum[i] / (n_samples_methods[i] + eps));
        // std::cout << "Method n°" << i << " has luminance variance of " << variances[i] << std::endl;
    }

    // if out of expected bounds then check contribution variance
    auto check_alpha = (alphas.at(0) <= 0) or (alphas.at(0) >= 1);
    if (dr::any_or<true>(check_alpha))
        alphas[0] = dr::any_or<true>(variances.at(0) >= variances.at(1)) ? 0.9f : 0.1f;

    alphas[1] = 1.f - alphas[0];

    // reset();

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

MI_VARIANT MISTsallis<Float, Spectrum>::MISTsallis(uint32_t n_methods, uint32_t batch_samples, Float gamma) 
    : Base(n_methods, batch_samples), gamma(gamma) {
    // expected only 2 sampling methods
    Assert(n_methods == 2); 
}

MI_VARIANT void MISTsallis<Float, Spectrum>::add_sampling_data(uint32_t sampling_method_id, 
    const Spectrum &luminance, 
    const std::vector<Float> &pdfs) {

    Float eps = std::numeric_limits<Float>::epsilon();

    for (uint32_t i = 0; i < n_methods; i++)
        if (pdfs[i] > max_pdf[i])
            max_pdf[i] = pdfs[i];

    // Classical track
    Float lum = 0.2126f * luminance.x() + 0.7152f * luminance.y() + 0.0722f * luminance.z();
    luminance_sum[sampling_method_id] += lum;

    // increase number of sample for this method
    n_samples_methods[sampling_method_id] += 1;

    Float mean_luminance = (luminance_sum[sampling_method_id] / n_samples_methods[sampling_method_id]);
    squared_sum[sampling_method_id] += (mean_luminance - lum) * (mean_luminance - lum);

    // Tsallis track
    // pdfs[0] is Light PDF
    // pdfs[1] is BSDF PDF

    // We expect:
    // p_1 is Light PDF
    // p_2 is BSDF PDF
    Float f_pdf = pdfs[0];
    Float g_pdf = pdfs[1];

    // Compute xi_sum and xi_prime_sum
    Float lum_tsallis = dr::pow(lum, gamma);

    // [MIS]: update xi and xi' with respect to equations 30-31
    Float alpha_probs = alphas[0] * f_pdf + alphas[1] * g_pdf;

    // check not null value
    if (dr::any_or<true>(alpha_probs <= 0.f))
        alpha_probs += eps;

    xi_sum += (lum_tsallis / dr::pow(alpha_probs, gamma + 1.f)) 
                    * (f_pdf - g_pdf);

    xi_prime_sum += (lum_tsallis / dr::pow(alpha_probs, gamma + 2.f)) 
                    * ((f_pdf - g_pdf) * (f_pdf - g_pdf));

    // std::cout << "---------------------" << std::endl;
    // std::cout << "@sample: " << n_samples << std::endl;
    // std::cout << " -- XiSum: " << xi_sum << std::endl;
    // std::cout << " -- XiPrimeSum: " << xi_prime_sum << std::endl;
}   

MI_VARIANT void MISTsallis<Float, Spectrum>::update_alphas() {
    
    // TODO: compute new alphas (depending of n_samples and batch)
    // if ((n_samples % batch_samples) != 0)
    //     return;

    Float eps = std::numeric_limits<Float>::epsilon();

    // Need to update UpdateProbsMIS, in order to take into account 
    // luminance without balance heuristic and then update internal data 
    // in order to compute xi_{\alpha}
    // uint32_t total_samples = 0;
    // for (uint32_t i = 0; i < this->n_methods; i++)
    //     total_samples += n_samples_methods[i];

    Float xi_alpha = xi_sum / batch_samples;

    Float xi_prime_alpha = xi_prime_sum * (-gamma / batch_samples);

    if (dr::any_or<true>(xi_prime_alpha <= 0))
        xi_prime_alpha += eps;

    // std::cout << "Alpha update @" << n_samples << std::endl;
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
    auto check_alpha = (alphas.at(0) <= 0) or (alphas.at(0) >= 1);
    if (dr::any_or<true>(check_alpha))
        alphas[0] = dr::any_or<true>(alphas.at(0) >= alphas.at(1)) ? 0.9f : 0.1f;
        
    // std::cout << " -- Real new Alpha MIS: " << alphas[0] << std::endl;
    alphas[1] = 1.f - alphas[0];

    // std::cout << " -- Adapted Alpha MIS: " << alphas[0] << std::endl;
}

MI_VARIANT void MISTsallis<Float, Spectrum>::reset() {

    // reset tracked data for new samples batch
    xi_sum = 0.f;
    xi_prime_sum = 0.f;

    for (uint32_t i = 0; i < this->n_methods; i++)
        n_samples_methods[i] = 0;
}

//! @}
// =======================================================================

MI_IMPLEMENT_CLASS_VARIANT(MISModel, Object, "MIS")
MI_IMPLEMENT_CLASS_VARIANT(MISBalance, MISModel, "MIS balance")
MI_IMPLEMENT_CLASS_VARIANT(MISPower, MISModel, "MIS power")
MI_IMPLEMENT_CLASS_VARIANT(MISDivergence, MISModel, "MIS Divergence")
MI_IMPLEMENT_CLASS_VARIANT(MISLight, MISModel, "MIS Light")
MI_IMPLEMENT_CLASS_VARIANT(MISBSDF, MISModel, "MIS BSDF")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear1, MISDivergence, "MIS Linear 1")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear2, MISDivergence, "MIS Linear 2")
MI_IMPLEMENT_CLASS_VARIANT(MISLinear3, MISDivergence, "MIS Linear 3")
MI_IMPLEMENT_CLASS_VARIANT(MISTsallis, MISDivergence, "MIS Tsallis")

MI_INSTANTIATE_CLASS(MISModel)
MI_INSTANTIATE_CLASS(MISBalance)
MI_INSTANTIATE_CLASS(MISPower)
MI_INSTANTIATE_CLASS(MISDivergence)
MI_INSTANTIATE_CLASS(MISLight)
MI_INSTANTIATE_CLASS(MISBSDF)
MI_INSTANTIATE_CLASS(MISLinear1)
MI_INSTANTIATE_CLASS(MISLinear2)
MI_INSTANTIATE_CLASS(MISLinear3)
MI_INSTANTIATE_CLASS(MISTsallis)

NAMESPACE_END(mitsuba)