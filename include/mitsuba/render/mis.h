#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/object.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/spectrum.h>
#include <mitsuba/render/fwd.h>
#include <map>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISModel : public Object {
public:

    void update_n_samples();

    Float get_alpha(uint32_t sampling_method_id) const;

    void add_sampling_data(uint32_t sampling_method_id, 
        const Spectrum &luminance, 
        const std::vector<Float> &pdfs); 

    virtual void update_alphas() = 0;

    // Virtual destructor
    virtual ~MISModel() { }

    MI_DECLARE_CLASS()

protected:

    MISModel(uint32_t n_methods);

    virtual void init_alphas();

    uint32_t n_methods;
    uint32_t n_samples; // total n_samples
    std::map<uint32_t, Float> alphas;
    std::map<uint32_t, uint32_t> n_samples_methods;
    std::map<uint32_t, Float> luminance_sum;
    std::map<uint32_t, std::vector<Float>> pdf_sum;
};

// Specific MIS types
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB MISTsallis : public MISModel<Float, Spectrum> {

    MI_IMPORT_BASE(MISModel)

    public:
        MISTsallis(uint32_t n_methods);

        void update_alphas() override;

        virtual ~MISTsallis();

        MI_DECLARE_CLASS()
};

MI_EXTERN_CLASS(MISModel)
MI_EXTERN_CLASS(MISTsallis)
NAMESPACE_END(mitsuba)