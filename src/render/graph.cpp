#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/render/graph.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GNNGraph implementations
// =======================================================================

MI_VARIANT GNNGraph<Float, Spectrum>::GNNGraph() : Object() {
};

MI_VARIANT GNNGraph<Float, Spectrum>::GNNGraph(Point3f origin, std::vector<Float> targets) : Object(), origin(origin), targets(targets) {
};

MI_VARIANT std::vector<typename GNNGraph<Float, Spectrum>::GNNNode> GNNGraph<Float, Spectrum>::get_nodes() const {
    return nodes;
};

MI_VARIANT std::vector<typename GNNGraph<Float, Spectrum>::GNNConnection> GNNGraph<Float, Spectrum>::get_connections() const {
    return connections;
};

MI_VARIANT typename GNNGraph<Float, Spectrum>::Point3f GNNGraph<Float, Spectrum>::get_origin() const {
    return origin;
};

MI_VARIANT std::vector<Float> GNNGraph<Float, Spectrum>::get_targets() const {
    return targets;
};

MI_VARIANT typename GNNGraph<Float, Spectrum>::GNNNode GNNGraph<Float, Spectrum>::get_node_by_index(uint32_t index) const {
    
    // TODO: check out of bounds
    return nodes.at(index);
};

MI_VARIANT int GNNGraph<Float, Spectrum>::get_node_index(const GNNNode &node) const {
    
    auto it = std::find(nodes.cbegin(), nodes.cend(), node);
  
    // If element was found
    if (it != nodes.end()) 
        return it - nodes.begin();
    else
        return -1;
};

MI_VARIANT void GNNGraph<Float, Spectrum>::set_origin(Point3f origin) {
    
    this->origin = origin;
};

MI_VARIANT void GNNGraph<Float, Spectrum>::set_targets(std::vector<Float> data) {
    
    this->targets = data;
};

MI_VARIANT bool GNNGraph<Float, Spectrum>::add_node(GNNNode node) {
    
    if (get_node_index(node) != -1) {
        nodes.push_back(node);
        return true;
    }
    else
        return false;
};

MI_VARIANT bool GNNGraph<Float, Spectrum>::add_connection(GNNConnection connection) {
    
    auto it = std::find(connections.cbegin(), connections.cend(), connection);
        
    if (it != connections.end()) {
        connections.push_back(connection);  
        return true;
    }
    
    return false;
};


MI_VARIANT std::vector<Float> GNNGraph<Float, Spectrum>::get_properties() const {
    return targets;
};

MI_VARIANT GNNGraph<Float, Spectrum>::~GNNGraph() {
};

//! @}
// =======================================================================

MI_IMPLEMENT_CLASS_VARIANT(GNNGraph, Object, "GNNGraph")

MI_INSTANTIATE_CLASS(GNNGraph)

NAMESPACE_END(mitsuba)
