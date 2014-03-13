#ifndef RAY_TRACER_H
#define RAY_TRACER_H

#include <vector>

#include "mesh.h"
#include "ray.h"
#include "triangle.h"
#include "bvh.h"
#include "options_ray_tracer.h"


class RayTracer
{
public:




public:
    // Constructor with options argument
    RayTracer(RayTracerOptions const& opts);
    
    // Traces the scene contained in a BVH sturcture and
    // writes the result into the 'image' vector
    void trace(BVH const& scene, std::vector<unsigned char>* image);


private:
    // compute shading as angle between ray direction and normal
    float shading(Ray const& ray, Vec3f const& normal);
    
    // compute actual ambient occlusion
    float ambientOcclusion(BVH const& scene,
        Vec3f const& point, Vec3f const& normal, int const pID);
    
    // compute smooth normal at a specific intersection point.
    // This is interpolated between the surrounding vertex normals
    // and depends on the actual intersection point inside the triangle.
    Vec3f getSmoothNormal(Intersection const& intersection);
    
    // compute flat normal at a specific intersection point.
    // This is the same normal for a whole triangle, no matter where
    // the ray actually hits.
    Vec3f getFlatNormal(Intersection const& intersection);

private:
    RayTracerOptions opts;
};

#endif // RAY_TRACER_H

