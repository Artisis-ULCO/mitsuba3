#pragma once

#include <mitsuba/mitsuba.h>
// #include <mitsuba/render/fwd.h>
#include <mitsuba/core/object.h>
// #include <mitsuba/core/spectrum.h>

#include <mitsuba/render/node.h>
#include <mitsuba/render/connection.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GNNGraph : public Object {
public:

    MI_IMPORT_TYPES(GNNNode, GNNConnection)

    // some accessors
    std::vector<GNNNode> get_nodes() const;
    std::vector<GNNConnection> get_connections() const;

    std::vector<Float> get_targets() const; 
    GNNNode get_node_by_index(uint32_t index) const;
    uint32_t get_node_index(const GNNNode &node) const;

    std::vector<GNNConnection> get_connections_from(const GNNNode &node) const {

        std::vector<GNNConnection> found;

        // TODO

        return found;
    }
    std::vector<GNNConnection> get_connections_to(const GNNNode &node) const {
        std::vector<GNNConnection> found;

        // TODO

        return found;
    }
    
    // set or add informations into graph
    void set_targets(std::vector<Float> data); 
    bool add_node(GNNNode node); 
    bool add_connection(GNNConnection connection); 

    std::vector<Float> get_properties() const; 

    // Virtual destructor
    virtual ~GNNGraph();

    MI_DECLARE_CLASS()

protected:
    GNNGraph();
    GNNGraph(std::vector<Float> targets);

protected:
    std::vector<GNNNode> nodes;
    std::vector<GNNConnection> connections;
    std::vector<Float> targets;
};

MI_EXTERN_CLASS(GNNGraph)
NAMESPACE_END(mitsuba)