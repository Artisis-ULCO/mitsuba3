#include <mitsuba/render/fwd.h>
// #include <mitsuba/render/node.h>
#include <mitsuba/render/graph.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GNNConnection implementations
// =======================================================================

MI_VARIANT GNNConnection<Float, Spectrum>::GNNConnection(GNNNode* from_node, GNNNode* to_node, std::vector<Float> data) 
    : Object(), from_node(from_node), to_node(to_node), data(data) {

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

MI_VARIANT GNNConnection<Float, Spectrum>::~GNNConnection() {
}

//! @}
// =======================================================================

MI_IMPLEMENT_CLASS_VARIANT(GNNConnection, Object, "GNNConnection")

MI_INSTANTIATE_CLASS(GNNConnection)


NAMESPACE_END(mitsuba)
