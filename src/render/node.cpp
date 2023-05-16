#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/node.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GNNNode implementations
// =======================================================================

MI_VARIANT GNNNode<Float, Spectrum>::GNNNode(Point3f position, Vector3f normal, Spectrum radiance) 
    : Object(), position(position), normal(normal), primary(false) {

    radiances.push_back(radiance);
}

MI_VARIANT GNNNode<Float, Spectrum>::GNNNode(Point3f position, Vector3f normal, Spectrum radiance, bool primary) 
    : Object(), position(position), normal(normal), primary(primary) {

    radiances.push_back(radiance);
}

MI_VARIANT typename GNNNode<Float, Spectrum>::Point3f GNNNode<Float, Spectrum>::get_position() const {
    return position;
};

MI_VARIANT typename GNNNode<Float, Spectrum>::Vector3f GNNNode<Float, Spectrum>::get_normal() const {
    return normal;
};

MI_VARIANT Spectrum GNNNode<Float, Spectrum>::get_radiance() const {
    
    Spectrum radiance = 0.f;

    for (auto r : radiances)
        radiance += r;

    if (radiances.size() > 0)
        return radiance / radiances.size();
    else 
        return radiance;
        
};

MI_VARIANT std::vector<Float> GNNNode<Float, Spectrum>::get_properties() const {
    std::vector<Float> values;
    return values;
}

MI_VARIANT GNNNode<Float, Spectrum>::~GNNNode() {
}
//! @}
// =======================================================================


MI_IMPLEMENT_CLASS_VARIANT(GNNNode, Object, "GNNNode")

MI_INSTANTIATE_CLASS(GNNNode)

NAMESPACE_END(mitsuba)
