#include <limits>
#include <stdexcept>
#include <iostream>

#include "visual.h"
#include "bvh_node.h"

namespace
{
    Vec3f getSplitPoint(AABB const& aabb)
    {
        switch (aabb.getLongestAxis())
        {
        case 0:
            return Vec3f((aabb.getAABBMax()[0] + aabb.getAABBMin()[0]) / 2.0f,
                 aabb.getAABBMax()[1],
                 aabb.getAABBMax()[2]);

        case 1:
            return Vec3f(aabb.getAABBMax()[0],
                (aabb.getAABBMax()[1] + aabb.getAABBMin()[1]) / 2.0f,
                aabb.getAABBMax()[2]);

        case 2:
            return Vec3f(aabb.getAABBMax()[0],
                aabb.getAABBMax()[1],
                (aabb.getAABBMax()[2] + aabb.getAABBMin()[2]) / 2.0f);

        default:
            break;
        }

        throw std::runtime_error("Invalid longest axis");
        return Vec3f();
    }
} // namespace

BVHNode::BVHNode()
{
    this->right = NULL;
    this->left = NULL;
}

BVHNode::~BVHNode()
{
    delete this->left;
    delete this->right;
}

void BVHNode::insert(Mesh const& mesh, std::vector<unsigned int>* faceIDs)
{
    // Compute AABB for both triangles and triangle centroids.
    AABB centroidAABB;
    for (std::size_t i = 0; i < faceIDs->size(); ++i)
    {
        unsigned int v1 = mesh.faces[(*faceIDs)[i] * 3 + 0];
        unsigned int v2 = mesh.faces[(*faceIDs)[i] * 3 + 1];
        unsigned int v3 = mesh.faces[(*faceIDs)[i] * 3 + 2];

        centroidAABB.merge((mesh.vertices[v1] + mesh.vertices[v2]
            + mesh.vertices[v3]) / 3.0f);
        this->aabb.merge(mesh.vertices[v1]);
        this->aabb.merge(mesh.vertices[v2]);
        this->aabb.merge(mesh.vertices[v3]);
    }

    // If there are only N triangles left, make this a leaf node.
    if (faceIDs->size() <= MAX_LEAF_TRIANGLES)
    {
        for (std::size_t i = 0; i < faceIDs->size(); ++i)
        {
            Triangle tri(&mesh, (*faceIDs)[i]);
            this->triangles.push_back(tri);
        }
        return;
    }

    Vec3f const sp = getSplitPoint(centroidAABB);
    AABB const aabbLeft(centroidAABB.getAABBMin(), sp);

    // Split index list into two index lists.
    std::vector<unsigned int> faceIDsLeft;
    std::vector<unsigned int> faceIDsRight;
    for (unsigned int i = 0; i < faceIDs->size(); ++i)
    {
        Triangle tri(&mesh, (*faceIDs)[i]);
        // Check if triangle belongs to left or right bounding box.
        if (aabbLeft.inside(tri.getCentroid()))
            faceIDsLeft.push_back((*faceIDs)[i]);
        else
            faceIDsRight.push_back((*faceIDs)[i]);
    }

    // The old index list is no longer needed. Clear it.
    faceIDs->clear();

    // Recurse into left and right child node.
    this->left = new BVHNode();
    this->left->insert(mesh, &faceIDsLeft);
    this->right = new BVHNode();
    this->right->insert(mesh, &faceIDsRight);
}

//test if a ray intersects this node and find the triangle with shortest distance
bool BVHNode::intersect(Ray const& ray, Intersection* intersection, float const maxDistance) const {
	if(!this->aabb.inside(ray.position)) {
		float distance= 0.0f;
		if(!this->aabb.intersect(ray, &distance))
			return false;
		if(distance > maxDistance)
			return false;
		if(distance > intersection->distance) {
			return false;
		}
	}
	if (!this->triangles.empty())
	{
		bool success = false;
		for (std::size_t i = 0; i < this->triangles.size(); ++i)
			success |= this->triangles[i].intersect(ray, intersection);
		return success;
	}

	float distanceLeft= 0.0f;
	if(!this->left->aabb.intersect(ray, &distanceLeft))
		return this->right->intersect(ray, intersection, maxDistance);
	float distanceRight= 0.0f;
	if(!this->right->aabb.intersect(ray, &distanceRight))
		return this->left->intersect(ray, intersection, maxDistance);
	if(distanceLeft < distanceRight) {
		bool const hitL = this->left->intersect(ray, intersection, maxDistance);
		bool const hitR = this->right->intersect(ray, intersection, maxDistance);
		return hitL || hitR;
	} else {
		bool const hitR = this->right->intersect(ray, intersection, maxDistance);
		bool const hitL = this->left->intersect(ray, intersection, maxDistance);
		return hitL || hitR;
	}
}
