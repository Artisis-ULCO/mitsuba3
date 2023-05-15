#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/object.h>

#include <mitsuba/render/connection.h>
#include <mitsuba/render/graph.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GraphContainer : public Object {
public:

    MI_IMPORT_TYPES(GNNGraph, GNNConnection)

    void add_graphs(std::vector<GNNGraph> graphs); 

    virtual void build_connections() = 0;

    // Virtual destructor
    virtual ~GraphContainer();

    MI_DECLARE_CLASS()

protected:
    GraphContainer(uint32_t n_graphs, uint32_t n_nodes_per_graphs, uint32_t n_neighbors);

protected:
    uint32_t n_graphs;
    uint32_t n_nodes_per_graphs;
    uint32_t n_neighbors;
    std::string reference;

    // Store also graphs and connections
    std::vector<GNNGraph> graphs;
    std::vector<GNNConnection> connections;
};

// Specific Container types: SimpleGraphContainer
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB SimpleGraphContainer : public GraphContainer<Float, Spectrum> {

    MI_IMPORT_BASE(GraphContainer)
    MI_IMPORT_TYPES()

public:
    SimpleGraphContainer(uint32_t n_graphs, uint32_t n_nodes_per_graphs, uint32_t n_neighbors);

    virtual void build_connections();

    MI_DECLARE_CLASS()
};


MI_EXTERN_CLASS(GraphContainer)
MI_EXTERN_CLASS(SimpleGraphContainer)
NAMESPACE_END(mitsuba)