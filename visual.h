#ifndef VISUAL_H_
#define VISUAL_H_

#include "vec3.h"
#include "mesh.h"
#include "ray.h"
#include "aabb.h"

// Debugging facility to sample a line from start to end.
void sampleLineToMesh(Mesh* mesh, int nrOfDots, Vec3f const& start, Vec3f const& end);

// Debugging facility to sample a AABB.
void sampleAABBToMesh(Mesh* mesh, int nrOfDots, AABB const& aabb);

// Debugging facility to sample a ray.
void sampleRayToMesh(Mesh* mesh, Ray const& ray, int nrOfDots, float length);

#endif /* VISUAL_H_ */
