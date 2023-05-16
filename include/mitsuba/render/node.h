#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/spectrum.h>

#include <mitsuba/core/object.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GNNNode : public Object {
public:

    MI_IMPORT_CORE_TYPES()
    
    Point3f get_position() const;
    Vector3f get_normal() const;
    Spectrum get_radiance() const;
    bool is_primary() const;

    virtual std::vector<Float> get_properties() const;

    // Virtual destructor
    virtual ~GNNNode();

    bool operator==(const GNNNode &other) const{
        return position == other.position and normal == other.normal;
    }

    GNNNode(Point3f position, Vector3f normal, Spectrum radiance);
    GNNNode(Point3f position, Vector3f normal, Spectrum radiance, bool primary);

    MI_DECLARE_CLASS()

protected:
    Point3f position;
    Vector3f normal;
    std::vector<Spectrum> radiances; // accumulate radiance
    bool primary;
};

MI_EXTERN_CLASS(GNNNode)
NAMESPACE_END(mitsuba)