#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/fwd.h>
#include <mitsuba/core/vector.h>
#include <mitsuba/core/spectrum.h>

#include <mitsuba/core/object.h>

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

    Point3f get_origin() const;
    std::vector<Float> get_targets() const; 
    GNNNode get_node_by_index(uint32_t index) const;
    int get_node_index(const GNNNode &node) const;

    std::vector<GNNConnection> get_connections(const GNNNode &node) const {

        std::vector<GNNConnection> found;

        if (get_node_index(node) != -1) {
            
            std::copy_if(connections.cbegin(), connections.cend(), std::back_inserter(found), [&node](const GNNConnection &c){
                return c.from() == node or c.to() == node;
            });
        }

        return found;

    }

    std::vector<GNNConnection> get_connections_from(const GNNNode &node) const {

        std::vector<GNNConnection> found;

        if (get_node_index(node) != -1) {
            
            std::copy_if(connections.cbegin(), connections.cend(), std::back_inserter(found), [&node](const GNNConnection &c){
                return c.from() == node;
            });
        }

        return found;
    }
    
    std::vector<GNNConnection> get_connections_to(const GNNNode &node) const {
        std::vector<GNNConnection> found;

        if (get_node_index(node) != -1) {
            
            std::copy_if(connections.cbegin(), connections.cend(), std::back_inserter(found), [&node](const GNNConnection &c){
                return c.to() == node;
            });
        }

        return found;
    }
    
    // set or add informations into graph
    void set_origin(Point3f origin); 
    void set_targets(std::vector<Float> data); 
    bool add_node(GNNNode node); 
    bool add_connection(GNNConnection connection); 

    std::vector<Float> get_properties() const; 

    // Virtual destructor
    virtual ~GNNGraph();

    GNNGraph();
    GNNGraph(Point3f origin, std::vector<Float> targets);

    MI_DECLARE_CLASS()

protected:
    Point3f origin;
    std::vector<GNNNode> nodes;
    std::vector<GNNConnection> connections;
    std::vector<Float> targets;
};

MI_EXTERN_CLASS(GNNGraph)
NAMESPACE_END(mitsuba)