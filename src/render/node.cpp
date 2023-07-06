#include <mitsuba/mitsuba.h>
#include <mitsuba/json.hpp>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/node.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GNNNode implementations
// =======================================================================

MI_VARIANT GNNNode<Float, Spectrum>::GNNNode(Point3f position, Vector3f normal, Spectrum radiance, uint32_t depth) 
    : Object(), position(position), normal(normal), depth(depth), primary(false) {

    radiances.push_back(radiance);
}

MI_VARIANT GNNNode<Float, Spectrum>::GNNNode(Point3f position, Vector3f normal, Spectrum radiance, uint32_t depth, bool primary) 
    : Object(), position(position), normal(normal), depth(depth), primary(primary) {

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

MI_VARIANT bool GNNNode<Float, Spectrum>::is_primary() const {
    return primary;
};

MI_VARIANT std::vector<Float> GNNNode<Float, Spectrum>::get_properties() const {
    std::vector<Float> values;

    // add position
    values.push_back(position.x());
    values.push_back(position.y());
    values.push_back(position.z());

    // add normal
    values.push_back(normal.x());
    values.push_back(normal.y());
    values.push_back(normal.z());

    // get radiance and push back channels
    Spectrum radiance = get_radiance();

    values.push_back(radiance.x());
    values.push_back(radiance.y());
    values.push_back(radiance.z());

    // also add depth
    values.push_back((Float)depth);

    return values;
}

MI_VARIANT nlohmann::json GNNNode<Float, Spectrum>::to_json() const {

    nlohmann::json json;

    auto properties = get_properties();

    // properties extraction
    nlohmann::json v_properties = nlohmann::json::array();
    for (auto prop : properties)
        v_properties.push_back(prop);

    json["attr"] = v_properties;

    // position extraction
    nlohmann::json v_position = nlohmann::json::array();
    v_position.push_back(position.x());
    v_position.push_back(position.y());
    v_position.push_back(position.z());

    json["pos"] = v_position;
    json["primary"] = primary;

    return json;
}

MI_VARIANT GNNNode<Float, Spectrum>::~GNNNode() {
}
//! @}
// =======================================================================


MI_IMPLEMENT_CLASS_VARIANT(GNNNode, Object, "GNNNode")

MI_INSTANTIATE_CLASS(GNNNode)

NAMESPACE_END(mitsuba)
