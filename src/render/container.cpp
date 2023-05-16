#include <mitsuba/render/graph.h>
#include <mitsuba/render/container.h>

#include <random>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GraphContainer implementations
// =======================================================================

MI_VARIANT GraphContainer<Float, Spectrum>::GraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors) 
    : Object(), build_at(build_at), n_nodes(n_nodes), n_neighbors(n_neighbors), n_samples(0) {

};

MI_VARIANT uint32_t GraphContainer<Float, Spectrum>::number_of_samples() const {
    return n_samples;
};

MI_VARIANT void GraphContainer<Float, Spectrum>::update_n_samples() {
    n_samples++;
};

MI_VARIANT bool GraphContainer<Float, Spectrum>::can_track() const {
    return n_samples < build_at;
};

MI_VARIANT bool GraphContainer<Float, Spectrum>::can_build() const {
    return n_samples == build_at;
};

MI_VARIANT uint32_t GraphContainer<Float, Spectrum>::number_of_connections() const {
    return connections.size();
};

MI_VARIANT void GraphContainer<Float, Spectrum>::add_graph(GNNGraph* graph) {

    for (auto &node : graph->get_nodes()) 
        add_node(node);

    for (auto &connection : graph->get_connections()) 
        add_connection(connection);
};

MI_VARIANT bool GraphContainer<Float, Spectrum>::add_node(GNNNode* node) {
    
    auto it = std::find(nodes.cbegin(), nodes.cend(), node);

    if (it == nodes.cend()) {
        nodes.push_back(node);
        return true;
    }

    return false;
};

MI_VARIANT bool GraphContainer<Float, Spectrum>::add_connection(GNNConnection* connection) {
    
    auto it = std::find(connections.cbegin(), connections.cend(), connection);
        
    if (it == connections.cend()) {
        connections.push_back(connection);  
        return true;
    }
    
    return false;
};

MI_VARIANT GraphContainer<Float, Spectrum>::~GraphContainer() {
}

//! @}
// =======================================================================


// =======================================================================
//! @{ \name SimpleGraphContainer implementations
// =======================================================================

MI_VARIANT SimpleGraphContainer<Float, Spectrum>::SimpleGraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors) 
    : Base(build_at, n_nodes, n_neighbors) {

};

MI_VARIANT void SimpleGraphContainer<Float, Spectrum>::build_connections(const Scene *scene) {
    
    if (this->nodes.size() < 2)
        return;

    // select a number of nodes for this graph
    std::vector<ref<GNNNode>> selected_nodes;

    for (uint32_t i = 0; i < this->n_nodes; i++) {
        auto node = nodes[rand() % nodes.size()];
        selected_nodes.push_back(node);
    }

    if (selected_nodes.size() > 0) {

        // try to create connection
        for (auto node : selected_nodes) {

            // select a number of neighbors nodes
            std::vector<ref<GNNNode>> neighbor_nodes;

            for (uint32_t i = 0; i < this->n_neighbors; i++) {
                auto neighbor_node = nodes[rand() % nodes.size()];

                if (neighbor_node != node)
                    neighbor_nodes.push_back(neighbor_node);
            }

            for (auto c_neighbor_node : neighbor_nodes) {
                // do necessary to create connection from near origin point
                if (node->is_primary() and c_neighbor_node->is_primary())
                    continue;

                Point3f origin = node->get_position();
                Point3f point = c_neighbor_node->get_position();

                Vector3f direction = point - origin;
                Vector3f normalized_d = direction / dr::sqrt(dr::sum(direction * direction));
                auto ray = Ray3f(origin, normalized_d);

                SurfaceInteraction3f si = scene->ray_intersect(ray, +RayFlags::All, true, /* active */ true);
        
                Float distance = dr::sqrt(dr::sum((point - origin) * (point - origin)));

                if (si.is_valid() and si.t >= distance) {
                    // add connection into current graph (from -> to)
                    auto from_to_connection = new GNNConnection(node, c_neighbor_node, {si.t});
                    
                    // add connection into neighbor graph (to -> from)
                    auto to_from_connection = new GNNConnection(c_neighbor_node, node, {si.t});

                    // also inside current container
                    connections.push_back(from_to_connection);
                    connections.push_back(to_from_connection);
                }
            }
        }
    }
};


//! @}
// =======================================================================



MI_IMPLEMENT_CLASS_VARIANT(GraphContainer, Object, "GraphContainer")
MI_IMPLEMENT_CLASS_VARIANT(SimpleGraphContainer, GraphContainer, "SimpleGraphContainer")

MI_INSTANTIATE_CLASS(GraphContainer)
MI_INSTANTIATE_CLASS(SimpleGraphContainer)

NAMESPACE_END(mitsuba)
