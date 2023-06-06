#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/core/object.h>

#include <mitsuba/render/node.h>
#include <mitsuba/render/connection.h>
#include <mitsuba/render/scene.h>

#include <mitsuba/json.hpp>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GraphContainer : public Object {
public:

    MI_IMPORT_TYPES(GNNNode, GNNConnection, Scene)

    uint32_t number_of_samples() const;
    Spectrum get_direct_radiance() const;
    Spectrum get_indirect_radiance() const;
    Spectrum get_direct_target() const;
    Spectrum get_indirect_target() const;

    void update_n_samples();
    bool can_track() const;
    bool can_build() const;
    nlohmann::json json() const;

    uint32_t number_of_connections() const;
    int get_node_index(const GNNNode* node) const;

    bool add_node(GNNNode* node);
    bool add_connection(GNNConnection* connection);
    void add_direct_radiance(Spectrum radiance);
    void add_indirect_radiance(Spectrum radiance);
    void accum_direct_target(Spectrum radiance);
    void accum_indirect_target(Spectrum radiance);

    void clear();

    virtual void build_connections(const Scene *scene) = 0;
    virtual void prepare_export() = 0;

    // Virtual destructor
    virtual ~GraphContainer();

    MI_DECLARE_CLASS()

protected:
    GraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors);

protected:
    uint32_t build_at;
    uint32_t n_nodes;
    uint32_t n_neighbors;
    uint32_t n_samples;
    bool export_done;
    std::string json_data_str;

    // accumulated radiances at this point
    Spectrum c_direct_radiance; // GNN radiance
    Spectrum c_indirect_radiance; // GNN indirect radiance
    Spectrum direct_target; // total accumulated radiance
    Spectrum indirect_target; // total accumulated radiance

    // Store also graphs, nodes and connections
    std::vector<ref<GNNNode>> nodes;
    std::vector<ref<GNNConnection>> connections;
};

// Specific Container types: SimpleGraphContainer
template <typename Float, typename Spectrum>
class MI_EXPORT_LIB SimpleGraphContainer : public GraphContainer<Float, Spectrum> {

    MI_IMPORT_BASE(GraphContainer, connections, nodes, export_done)
    MI_IMPORT_TYPES(Scene, GNNNode, GNNConnection)

public:
    SimpleGraphContainer(uint32_t build_at, uint32_t n_nodes, uint32_t n_neighbors);

    virtual void build_connections(const Scene *scene) override;
    virtual void prepare_export() override;

    MI_DECLARE_CLASS()
};


MI_EXTERN_CLASS(GraphContainer)
MI_EXTERN_CLASS(SimpleGraphContainer)
NAMESPACE_END(mitsuba)