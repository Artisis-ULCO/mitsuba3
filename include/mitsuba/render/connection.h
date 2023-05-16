#pragma once

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/fwd.h>
#include <mitsuba/core/object.h>
// #include <mitsuba/core/spectrum.h>
#include <mitsuba/render/node.h>

NAMESPACE_BEGIN(mitsuba)

template <typename Float, typename Spectrum>
class MI_EXPORT_LIB GNNConnection : public Object {
public:

    // MI_IMPORT_RENDER_BASIC_TYPES()
    MI_IMPORT_TYPES(GNNNode)

    GNNNode from() const; 
    GNNNode to() const; 
    virtual std::vector<Float> get_properties() const; 

    // Virtual destructor
    virtual ~GNNConnection();

    bool operator==(const GNNConnection &other) const{
        return from_node == other.from_node and to_node == other.to_node;
    }
   
    GNNConnection(GNNNode from_node, GNNNode to_node, std::vector<Float> data);

    MI_DECLARE_CLASS()

protected:
    GNNNode from_node;
    GNNNode to_node;
    std::vector<Float> data;
};


MI_EXTERN_CLASS(GNNConnection)
NAMESPACE_END(mitsuba)