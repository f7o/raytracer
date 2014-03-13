#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "vec3.h"
#include "aabb.h"
#include "ray.h"
#include "mesh.h"

// Representation of a mesh/ray intersection.
//
// - mesh    : mesh that was intersected
// - faceID  : ID of the triangle that was intersected
// - bary    : barycentric coordinates of the intersection
// - position: The position of the intersection point
// - distance: The distance from the ray origin to the intersection point.
//     This field is also used to test whether an intersection happened.
//     In case of failure, it is set to 0.0f.
struct Intersection
{
    Mesh const* mesh;
    unsigned int faceID;
    Vec3f bary;
    Vec3f position;
    float distance;
};

// Representation of a triangle which contains an intersection test
// and a few convenience functions.
class Triangle
{
public:
    Triangle(void);
    Triangle(Mesh const* mesh, unsigned int faceID);

    Vec3f getCentroid(void) const;
    Vec3f getNormalVector(void) const;
    Vec3f getAABBMin(void) const;
    Vec3f getAABBMax(void) const;
    AABB getAABB(void) const;

    // Computes an intersection between a ray and the triangle.
    // Returns true iff the ray hits the triangle. If this case, additional
    // information is returned in the intersection field.
    bool intersect(Ray const& ray, Intersection* intersection) const;

    Vec3f const& operator[] (int index) const;
    Vec3f& operator[] (int index);
private:
    Mesh const* mesh;
    unsigned int faceID;
    // Three points of this triangle in the mesh
    Vec3f p1;
    Vec3f p2;
    Vec3f p3;
    // Precomputations
    Vec3f u, v, n;
    float uu, uv, vv, D;
};
#endif /* TRIANGLE_H_ */
