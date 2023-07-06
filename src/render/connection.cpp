#include <mitsuba/json.hpp>
#include <mitsuba/render/fwd.h>
// #include <mitsuba/render/node.h>
#include <mitsuba/render/connection.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GNNConnection implementations
// =======================================================================

MI_VARIANT GNNConnection<Float, Spectrum>::GNNConnection(GNNNode* from_node, GNNNode* to_node, std::vector<Float> data) 
    : Object(), from_node(from_node), to_node(to_node), data(data), built(false) {

};

MI_VARIANT GNNConnection<Float, Spectrum>::GNNConnection(GNNNode* from_node, GNNNode* to_node, std::vector<Float> data, bool built) 
    : Object(), from_node(from_node), to_node(to_node), data(data), built(built) {

};

MI_VARIANT ref<typename GNNConnection<Float, Spectrum>::GNNNode> GNNConnection<Float, Spectrum>::from() const {
    return from_node;
};

MI_VARIANT ref<typename GNNConnection<Float, Spectrum>::GNNNode> GNNConnection<Float, Spectrum>::to() const {
    return to_node;
};

MI_VARIANT std::vector<Float> GNNConnection<Float, Spectrum>::get_properties() const {
    return data;
};

MI_VARIANT nlohmann::json GNNConnection<Float, Spectrum>::to_json() const {

    nlohmann::json json;

    auto properties = get_properties();

    // properties extraction
    nlohmann::json v_properties = nlohmann::json::array();
    
    for (auto prop : properties)
        v_properties.push_back(prop);
    json["attr"] = v_properties;
    json["built"] = built;

    return json;
}

MI_VARIANT GNNConnection<Float, Spectrum>::~GNNConnection() {
}

//! @}
// =======================================================================

MI_IMPLEMENT_CLASS_VARIANT(GNNConnection, Object, "GNNConnection")

MI_INSTANTIATE_CLASS(GNNConnection)


NAMESPACE_END(mitsuba)
