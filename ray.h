#ifndef RAY_H_
#define RAY_H_

#include <xmmintrin.h>
#include "vec3.h"

// Simple representation of a ray.
struct Ray
{
    Vec3f position;
    Vec3f direction;
    Vec3f div;
};

#endif /* RAY_H_ */
