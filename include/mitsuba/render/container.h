#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/object.h>

#include <mitsuba/render/connection.h>
#include <mitsuba/render/graph.h>
#include <mitsuba/render/scene.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GraphContainer : public Object {
public:

    MI_IMPORT_TYPES(GNNNode, GNNGraph, GNNConnection, Scene)

    uint32_t number_of_samples() const;
    void update_n_samples();
    bool can_track() const;
    bool can_build() const;

    uint32_t number_of_connections() const;

    void add_graph(GNNGraph* graph);

    virtual void build_connections(const Scene *scene) = 0;

    // Virtual destructor
    virtual ~GraphContainer();

    MI_DECLARE_CLASS()

private:
    bool add_node(GNNNode* node);
    bool add_connection(GNNConnection* connection);

protected:
    GraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors);

protected:
    uint32_t build_at;
    uint32_t n_nodes;
    uint32_t n_neighbors;
    uint32_t n_samples;
    std::string reference;

    // Store also graphs and connections
    std::vector<ref<GNNNode>> nodes;
    std::vector<ref<GNNConnection>> connections;
};

// Specific Container types: SimpleGraphContainer
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB SimpleGraphContainer : public GraphContainer<Float, Spectrum> {

    MI_IMPORT_BASE(GraphContainer, connections, nodes)
    MI_IMPORT_TYPES(Scene, GNNNode, GNNGraph, GNNConnection)

public:
    SimpleGraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors);

    virtual void build_connections(const Scene *scene);

    MI_DECLARE_CLASS()
};


MI_EXTERN_CLASS(GraphContainer)
MI_EXTERN_CLASS(SimpleGraphContainer)
NAMESPACE_END(mitsuba)