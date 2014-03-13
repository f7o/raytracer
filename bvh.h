#ifndef BVH_H_
#define BVH_H_

#include <vector>

#include "mesh.h"
#include "triangle.h"
#include "ray.h"
#include "aabb.h"
#include "bvh_node.h"
#include "cl_bvh.h"
// Bounding Volume Hierarchy (BVH) Interface.
class BVH
{
public:
    BVH(void);
    ~BVH(void);

    // Constructs the BVH from the given mesh.
    void buildBVH(Mesh const& mesh);

    // Computes an intersection between a ray and the BVH.
    // In case of a successful intersection, the method returns true and
    // the result is stored in the argument 'intersection'.
    bool intersect(Ray const& ray, Intersection* intersection, float const maxDistance) const;
    
    // Traverses the BVH, collects and prints statistics to the given stream.
    int printStatistics(std::ostream& out);

    // Samples the edges of AABBs for visualization and adds it to the mesh.
    // Warning: Can produce huge meshes.
    void addAABBsToMesh(Mesh* mesh, int nrOfDots) const;

    int parseBVH(BVHNode const* node, struct cl_bvh_node* scene_array, int &counter) const;

    void getBVH(cl_bvh_node* scene_array) const;



private:
    BVHNode* root;
};

/* ------------------------ Implementation ------------------------ */

inline BVH::BVH()
{
    this->root = NULL;
}

inline BVH::~BVH()
{
    delete this->root;
}
#endif // BVH_H_
