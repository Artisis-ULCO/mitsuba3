#include <mitsuba/render/graph.h>
#include <mitsuba/render/container.h>

NAMESPACE_BEGIN(mitsuba)

// =======================================================================
//! @{ \name GraphContainer implementations
// =======================================================================

MI_VARIANT GraphContainer<Float, Spectrum>::GraphContainer(uint32_t n_graphs, uint32_t n_nodes_per_graphs, uint32_t n_neighbors) 
    : Object(), n_graphs(n_graphs), n_nodes_per_graphs(n_nodes_per_graphs), n_neighbors(n_neighbors) {

};

MI_VARIANT void GraphContainer<Float, Spectrum>::add_graphs(std::vector<GNNGraph> graphs) {

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

MI_VARIANT void SimpleGraphContainer<Float, Spectrum>::build_connections() {
    
};


//! @}
// =======================================================================



MI_IMPLEMENT_CLASS_VARIANT(GraphContainer, Object, "GraphContainer")
MI_IMPLEMENT_CLASS_VARIANT(SimpleGraphContainer, GraphContainer, "SimpleGraphContainer")

MI_INSTANTIATE_CLASS(GraphContainer)
MI_INSTANTIATE_CLASS(SimpleGraphContainer)

NAMESPACE_END(mitsuba)
