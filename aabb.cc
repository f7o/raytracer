#include "aabb.h"
#include <iostream>

void AABB::merge(AABB const& bb)
{
    for (int i = 0; i < 3; ++i)
    {
        this->min[i] = std::min(this->min[i], bb.min[i]);
        this->max[i] = std::max(this->max[i], bb.max[i]);
    }
}

void AABB::merge(Vec3f const& vec)
{
    for (int i = 0; i < 3; ++i)
    {
        this->min[i] = std::min(this->min[i], vec[i]);
        this->max[i] = std::max(this->max[i], vec[i]);
    }
}

int AABB::getLongestAxis() const
{
    Vec3f diff = this->max - this->min;
    if (diff[0] >= diff[1] && diff[0] >= diff[2])
        return 0;
    if (diff[1] >= diff[0] && diff[1] >= diff[2])
        return 1;
    return 2;
}

bool AABB::inside(Vec3f const& point) const
{
    for (int i = 0; i < 3; ++i)
        if (point[i] > this->max[i] || point[i] < this->min[i])
            return false;
    return true;
}

bool AABB::intersect(Ray const& ray, float* intersectDistance) const
{
    // Smits ray-box intersection test using slabs
    // http://www.cs.utah.edu/~awilliam/box/box.pdf
	// *** SSE ***
	__m128 min_ps= _mm_set_ps(0.0f, min[2], min[1], min[0]);
	__m128 max_ps= _mm_set_ps(0.0f, max[2], max[1], max[0]);
	__m128 position_ps= _mm_set_ps(0.0f, ray.position[2], ray.position[1], ray.position[0]);
	__m128 div_ps= _mm_set_ps(0.0f, ray.div[2], ray.div[1], ray.div[0]);

	// Compute intersections
    __m128 t1_ps= _mm_mul_ps(_mm_sub_ps(min_ps, position_ps), div_ps);
    __m128 t2_ps= _mm_mul_ps(_mm_sub_ps(max_ps, position_ps), div_ps);

    __m128 tmin_ps= _mm_min_ps(t1_ps, t2_ps);
    __m128 tmax_ps= _mm_max_ps(t1_ps, t2_ps);
	float *tmin_a= (float *) &tmin_ps;
	float *tmax_a= (float *) &tmax_ps;

	if (tmin_a[0] > tmax_a[1] || tmin_a[1] > tmax_a[0])
		return false;

	float tmin = std::max(tmin_a[0], tmin_a[1]);
	float tmax = std::min(tmax_a[0], tmax_a[1]);

	if (tmin > tmax_a[2] || tmin_a[2] > tmax)
		return false;

	tmin = std::max(tmin, tmin_a[2]);
	tmax = std::min(tmax, tmax_a[2]);

    *intersectDistance= fabsf(tmin);
    return tmin < 100000.0f && tmax > 0.0f;
}
