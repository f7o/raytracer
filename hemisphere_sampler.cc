#include "hemisphere_sampler.h"
#include <iostream>

HemisphereSampler::HemisphereSampler(Vec3f const& normal, const int pID)
{
    this->basisY = normal.normalized();
    Vec3f h = this->basisY;
    if (fabs(h[0])<=fabs(h[1]) && fabs(h[0])<=fabs(h[2]))
        h[0]= 1.0;
    else if (fabs(h[1])<=fabs(h[0]) && fabs(h[1])<=fabs(h[2]))
        h[1]= 1.0;
    else
        h[2]= 1.0;
    
    this->basisX = h.cross(this->basisY);
    this->basisX = this->basisX.normalized();
    this->basisZ = this->basisX.cross(this->basisY);
    this->basisZ = this->basisZ.normalized();
    
    this->randGenrator.initialize(536870923u * pID);
}

Vec3f
HemisphereSampler::sample()
{
    /* use better random generator */
    float xi1 = this->randGenrator.randFloat();
    float xi2 = this->randGenrator.randFloat();
    
    float theta = acos(sqrt(1.0-xi1));
    float phi = 2.0 * PI * xi2;
    
    float xs = sinf(theta) * cosf(phi);
    float ys = cosf(theta);
    float zs = sinf(theta) * sinf(phi);
    
    Vec3f direction = this->basisX * xs + this->basisY * ys + this->basisZ * zs;
    return direction.normalized();
}
Vec3f
HemisphereSampler::randomFloat() {
	/* use better random generator */
	    float xi1 = this->randGenrator.randFloat();
	    float xi2 = this->randGenrator.randFloat();

	    float theta = acos(sqrt(1.0-xi1));
	    float phi = 2.0 * PI * xi2;

	    float xs = sinf(theta) * cosf(phi);
	    float ys = cosf(theta);
	    float zs = sinf(theta) * sinf(phi);
	    return Vec3f(xs,ys,zs);
}

