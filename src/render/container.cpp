#include <mitsuba/render/graph.h>
#include <mitsuba/render/container.h>

#include <random>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GraphContainer implementations
// =======================================================================

MI_VARIANT GraphContainer<Float, Spectrum>::GraphContainer(uint32_t n_graphs, uint32_t n_nodes_per_graphs, uint32_t n_neighbors) 
    : Object(), n_graphs(n_graphs), n_nodes_per_graphs(n_nodes_per_graphs), n_neighbors(n_neighbors) {

};

MI_VARIANT uint32_t GraphContainer<Float, Spectrum>::number_of_samples() const {
    return n_samples;
};

MI_VARIANT void GraphContainer<Float, Spectrum>::update_n_samples() {
    n_samples++;
};

MI_VARIANT void GraphContainer<Float, Spectrum>::add_graph(GNNGraph graph) {

    this->graphs.push_back(graph);
};

MI_VARIANT void GraphContainer<Float, Spectrum>::add_graphs(std::vector<GNNGraph> graphs) {

    for (auto graph : graphs)
        this->graphs.push_back(graph);
};

MI_VARIANT GraphContainer<Float, Spectrum>::~GraphContainer() {
}

//! @}
// =======================================================================


// =======================================================================
//! @{ \name SimpleGraphContainer implementations
// =======================================================================

MI_VARIANT SimpleGraphContainer<Float, Spectrum>::SimpleGraphContainer(uint32_t n_graphs, uint32_t n_nodes_per_graphs, uint32_t n_neighbors) 
    : Base(n_graphs, n_nodes_per_graphs, n_neighbors) {

};

MI_VARIANT void SimpleGraphContainer<Float, Spectrum>::build_connections(const Scene *scene) {
    
    // randomly select graphs
    std::vector<typename SimpleGraphContainer<Float, Spectrum>::GNNGraph> selected_graphs;

    for (uint32_t i = 0; i < this->n_graphs; i++) {
        auto selected = this->graphs[rand() % this->graphs.size()];
        selected_graphs.push_back(selected);
    }


    // for graph in random.choices(pos_graphs, k=n_graphs):
                
    //     selected_nodes = random.choices(graph.nodes, k=n_nodes_per_graphs)
    //     potential_neighbors = [ g for g in pos_graphs if g is not graph ]
        
    //     # check if there is at least 1 potential neighbor
    //     if len(potential_neighbors) > 0:
    //         neighbors_graphs = random.choices([ g for g in pos_graphs if g is not graph], \
    //                                 k=n_neighbors)

    //         # try now to create connection
    //         for node in selected_nodes:
                
    //             # select randomly one neighbor graph
    //             selected_graph = random.choice(neighbors_graphs)

    //             # randomly select current neighbor graph node for the connection
    //             neighbor_selected_node = random.choice(selected_graph.nodes)

    //             # do not continue the selected nodes are primary rays or origin
    //             if node.primary and neighbor_selected_node.primary:
    //                 continue
                
    //             # create Ray from current node
    //             origin, point = mi.Vector3f(node.position), mi.Vector3f(neighbor_selected_node.position)

    //             # get direction and create new ray
    //             direction = point - origin
    //             normalized_d = direction / np.sqrt(np.sum(direction ** 2))
    //             ray = mi.Ray3f(origin, normalized_d)

    //             # try intersect using this ray
    //             si = scene.ray_intersect(ray)
                
    //             # if connections exists, then the node is also attached to the graph
    //             # new bi-directionnal connections are created between `node` and 
    //             #  `neighbor_selected_node` with distance data
    //             if si.is_valid() and si.t >= math.dist(point, origin):

    //                 # add connection into current graph (from -> to)
    //                 connection = RayConnection(node, neighbor_selected_node, \
    //                     {'distance': si.t}, ConnectionTag.BUILT)
    //                 self._n_built_nodes += graph.add_node(neighbor_selected_node)
    //                 self._n_built_connections += graph.add_connection(connection)

    //                 # add connection into neighbor graph (to -> from)
    //                 connection = RayConnection(neighbor_selected_node, node, \
    //                     {'distance': si.t}, ConnectionTag.BUILT)
    //                 self._n_built_nodes += selected_graph.add_node(node)
    //                 self._n_built_connections += selected_graph.add_connection(connection)
                    
    // return True
};


//! @}
// =======================================================================



MI_IMPLEMENT_CLASS_VARIANT(GraphContainer, Object, "GraphContainer")
MI_IMPLEMENT_CLASS_VARIANT(SimpleGraphContainer, GraphContainer, "SimpleGraphContainer")

MI_INSTANTIATE_CLASS(GraphContainer)
MI_INSTANTIATE_CLASS(SimpleGraphContainer)

NAMESPACE_END(mitsuba)
