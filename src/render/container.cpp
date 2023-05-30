#include <mitsuba/json.hpp>
#include <mitsuba/render/container.h>

#include <random>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GraphContainer implementations
// =======================================================================

MI_VARIANT GraphContainer<Float, Spectrum>::GraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors) 
    : Object(), build_at(build_at), n_nodes(n_nodes), n_neighbors(n_neighbors), n_samples(0), export_done(false), target(0.f) {

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

MI_VARIANT nlohmann::json GraphContainer<Float, Spectrum>::json() const {
    return  nlohmann::json::parse(json_data_str);
}

MI_VARIANT uint32_t GraphContainer<Float, Spectrum>::number_of_connections() const {
    return connections.size();
};

MI_VARIANT void GraphContainer<Float, Spectrum>::add_radiance(Spectrum radiance) {

    // keep track of obtained radiances
    Float w = 1.f / (Float)(n_samples + 1);
    c_radiance = (1.f - w) * c_radiance + w * radiance;
};

MI_VARIANT void GraphContainer<Float, Spectrum>::accum_target(Spectrum radiance) {

    // keep track of obtained radiance
    Float w = 1.f / (Float)(n_samples + 1);
    target = (1.f - w) * target + w * radiance;
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

MI_VARIANT int GraphContainer<Float, Spectrum>::get_node_index(const GNNNode* node) const {
    
    auto it = std::find(nodes.cbegin(), nodes.cend(), node);
  
    // If element was found
    if (it != nodes.end()) 
        return it - nodes.begin();
    else
        return -1;
};

MI_VARIANT Spectrum GraphContainer<Float, Spectrum>::get_target() const {
    
    return target;
};

MI_VARIANT Spectrum GraphContainer<Float, Spectrum>::get_radiance() const {
    
    return c_radiance;    
};

MI_VARIANT void GraphContainer<Float, Spectrum>::clear() {

    // clear data
    connections.clear();
    nodes.clear();

    // reinit string
    this->json_data_str = "";

    // TODO: ensure remove of container once its done
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
    
    // if export already done, not necessary to continue
    if (export_done or nodes.size() < 2)
        return;

    // // filter primary node (only one use)
    // std::vector<ref<GNNNode>> primary_nodes;

    // // expected at least `n_samples` primary nodes
    // for (auto c_node : nodes)
    //     if (c_node->is_primary())
    //         primary_nodes.push_back(c_node);

    // // select primary node
    // auto c_primary_node = primary_nodes[rand() % nodes.size()];

    // // TODO: switch connection using this primary nodes
    // // TODO: remove all others primary nodes

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

                // only use the current kept primary node

                // not necessary to create connection from near origin point
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
                    // set to built
                    auto from_to_connection = new GNNConnection(node, c_neighbor_node, {si.t}, true);
                    
                    // add connection into neighbor graph (to -> from)
                    // set to built
                    auto to_from_connection = new GNNConnection(c_neighbor_node, node, {si.t}, true);

                    // also inside current container
                    connections.push_back(from_to_connection);
                    connections.push_back(to_from_connection);
                }
            }
        }
    }
};

MI_VARIANT void SimpleGraphContainer<Float, Spectrum>::prepare_export() {

    // create export data to json
    nlohmann::json data;

    // export nodes data
    nlohmann::json nodes_x;
    nlohmann::json nodes_pos;
    nlohmann::json nodes_primary;
    for (auto node : nodes) {
        auto node_data = node->to_json();
        nodes_x.push_back(node_data["attr"]);
        nodes_pos.push_back(node_data["pos"]);
        nodes_primary.push_back(node_data["primary"]);
    }

    data["x"] = nodes_x;
    data["x_primary"] = nodes_primary;
    data["pos"] = nodes_pos;

    // export edges data
    nlohmann::json edges_indices;
    nlohmann::json edges_attr;
    nlohmann::json edges_built;
    for (auto connection : connections) {

        // get nodes indices
        int node_index_from = this->get_node_index(connection->from());
        int node_index_to = this->get_node_index(connection->to());

        edges_indices.push_back(nlohmann::json::array({node_index_from, node_index_to}));

        // retrieve connection data 
        auto connection_data = connection->to_json();
        edges_attr.push_back(connection_data["attr"]);
        edges_built.push_back(connection_data["built"]);
    }

    data["edge_index"] = edges_indices;
    data["edge_attr"] = edges_attr;
    data["edge_built"] = edges_built;

    // store current json representation before clearing
    this->json_data_str = data.dump();

    // mark export as done
    export_done = true;
};


//! @}
// =======================================================================



MI_IMPLEMENT_CLASS_VARIANT(GraphContainer, Object, "GraphContainer")
MI_IMPLEMENT_CLASS_VARIANT(SimpleGraphContainer, GraphContainer, "SimpleGraphContainer")

MI_INSTANTIATE_CLASS(GraphContainer)
MI_INSTANTIATE_CLASS(SimpleGraphContainer)

NAMESPACE_END(mitsuba)
