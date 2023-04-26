#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/spectrum.h>
#include <map>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISModel : public Object {
public:

    MI_IMPORT_TYPES()

    void update_n_samples();

    uint32_t number_of_samples() const;

    uint32_t number_of_methods() const;

    virtual Float get_alpha(uint32_t sampling_method_id) const;

    virtual void add_sampling_data(uint32_t sampling_method_id, 
        const Spectrum &luminance, 
        const std::vector<Float> &pdfs); 

    virtual void update_alphas() = 0;

    virtual Float mis_weight(Float alpha, Float pdf_a, Float pdf_b) const = 0; 

    // TODO: Check (expected protected) 
    // Virtual destructor
    virtual ~MISModel();

    std::vector<Float> get_pdfs(uint32_t sampling_method_id) const;

    uint32_t number_of_samples_method(uint32_t sampling_method_id) const;

    MI_DECLARE_CLASS()

protected:
    MISModel(uint32_t n_methods);

    virtual void init_alphas();

protected:
    uint32_t n_methods;
    uint32_t n_samples; // total n_samples
    std::map<uint32_t, Float> alphas;
    std::map<uint32_t, uint32_t> n_samples_methods;
    std::map<uint32_t, Float> luminance_sum;
    std::map<uint32_t, Float> squared_sum;
    std::map<uint32_t, std::vector<Float>> pdf_sum;
};

// Specific MIS types: Balance
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISBalance : public MISModel<Float, Spectrum> {

    MI_IMPORT_BASE(MISModel)
    MI_IMPORT_TYPES()

public:
    MISBalance(uint32_t n_methods);

    virtual void update_alphas() override;
    Float mis_weight(Float alpha, Float pdf_a, Float pdf_b) const override;

    MI_DECLARE_CLASS()
};

// Specific MIS types: Power
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISPower : public MISModel<Float, Spectrum> {

    MI_IMPORT_BASE(MISModel)
    MI_IMPORT_TYPES()

public:
    MISPower(uint32_t n_methods);

    void update_alphas() override;
    Float mis_weight(Float alpha, Float pdf_a, Float pdf_b) const override;

    MI_DECLARE_CLASS()
};

// Specific MIS types: Light only
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISLight : public MISBalance<Float, Spectrum> {

    MI_IMPORT_BASE(MISBalance, alphas)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISLight(uint32_t n_methods);

    void update_alphas() override;
};

// Specific MIS types: BSDF only
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISBSDF : public MISBalance<Float, Spectrum> {

    MI_IMPORT_BASE(MISBalance, alphas)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISBSDF(uint32_t n_methods);

    void update_alphas() override;
};

// Specific MIS types: abstract divergence
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISDivergence : public MISModel<Float, Spectrum> {

    MI_IMPORT_BASE(MISModel)
    MI_IMPORT_TYPES()

public:
    virtual Float mis_weight(Float alpha, Float pdf_a, Float pdf_b) const override;

    MI_DECLARE_CLASS()
    
protected:
    MISDivergence(uint32_t n_methods, uint32_t batch_samples);

    virtual void update_alphas() = 0;
    virtual void reset() = 0;

protected:
    uint32_t batch_samples;
};

// Specific MIS types: Linear1 divergence
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISLinear1 : public MISDivergence<Float, Spectrum> {

    MI_IMPORT_BASE(MISDivergence, luminance_sum, pdf_sum, n_samples_methods, 
                squared_sum, alphas, n_samples, batch_samples)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISLinear1(uint32_t n_methods, uint32_t batch_samples);

    void update_alphas() override;

protected:
    void reset() override {};
};

// Specific MIS types: Linear2 divergence
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISLinear2 : public MISDivergence<Float, Spectrum> {

    MI_IMPORT_BASE(MISDivergence, luminance_sum, pdf_sum, n_samples_methods, 
                squared_sum, alphas, n_samples, batch_samples)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISLinear2(uint32_t n_methods, uint32_t batch_samples);

    void update_alphas() override;

protected:
    void reset() override {};
};

// Specific MIS types: Linear3 divergence
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISLinear3 : public MISDivergence<Float, Spectrum> {

    MI_IMPORT_BASE(MISDivergence, luminance_sum, pdf_sum, n_samples_methods, 
                squared_sum, alphas, n_samples, batch_samples)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISLinear3(uint32_t n_methods, uint32_t batch_samples);

    void add_sampling_data(uint32_t sampling_method_id, 
        const Spectrum &luminance, 
        const std::vector<Float> &pdfs);

    void update_alphas() override;

protected:
    void reset() override {};

private:
    std::map<uint32_t, Float> sum_mu;
};

// Specific MIS types: Tsallis divergence
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISTsallis : public MISDivergence<Float, Spectrum> {

    MI_IMPORT_BASE(MISDivergence, luminance_sum, pdf_sum, n_samples_methods, 
                squared_sum, alphas, n_samples, batch_samples)
    MI_IMPORT_TYPES()

public:
    MI_DECLARE_CLASS()
    
    MISTsallis(uint32_t n_methods, uint32_t batch_samples, Float gamma);

    Float get_alpha(uint32_t sampling_method_id) const;

    void add_sampling_data(uint32_t sampling_method_id, 
        const Spectrum &luminance, 
        const std::vector<Float> &pdfs);

    void update_alphas() override;

protected:
    void reset() override;

private:
    Float gamma;
    Float xi_sum;
    Float xi_prime_sum;
};

MI_EXTERN_CLASS(MISModel)
MI_EXTERN_CLASS(MISBalance)
MI_EXTERN_CLASS(MISPower)
MI_EXTERN_CLASS(MISDivergence)
MI_EXTERN_CLASS(MISLight)
MI_EXTERN_CLASS(MISBSDF)
MI_EXTERN_CLASS(MISLinear1)
MI_EXTERN_CLASS(MISLinear2)
MI_EXTERN_CLASS(MISLinear3)
MI_EXTERN_CLASS(MISTsallis)
NAMESPACE_END(mitsuba)