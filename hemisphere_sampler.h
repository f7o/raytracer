#ifndef HEMISPHERE_SAMPLER_H
#define HEMISPHERE_SAMPLER_H

#include "vec3.h"
#include "random.h"

#define PI 3.1415926535897932384626433832795


// Hemisphere Sampler
// Samples random cosine-weighted directions on a hemisphere
// for a given normal direction
class HemisphereSampler
{
public:
    // Constructor with normal direction.
    // Sets up a basis for the given intersection normal
    // that is used for consecutive calls of 'Vec3f sample(void);'
    HemisphereSampler(Vec3f const& normal, const int pID);
    
    // for each call returns one random direction on the hemisphere
    // based on the normal direction given in the constructor
    Vec3f sample(void);
    Vec3f randomFloat(void);

private:
    Vec3f basisX;
    Vec3f basisY;
    Vec3f basisZ;
    
    static unsigned int seed;
    RandomGenerator randGenrator;
};

#endif // HEMISPHERE_SAMPLER_H

