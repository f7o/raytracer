#include "triangle.h"
#include <iostream>

Triangle::Triangle(void) {
	this->mesh= NULL;
	this->faceID= -1;
}

Triangle::Triangle(Mesh const* mesh, unsigned int faceID) {
	this->mesh= mesh;
	this->faceID= faceID;
	this->p1= this->mesh->vertices[mesh->faces[faceID * 3 + 0]];
	this->p2= this->mesh->vertices[mesh->faces[faceID * 3 + 1]];
	this->p3= this->mesh->vertices[mesh->faces[faceID * 3 + 2]];

    // Precomputations used in intersection test
    u = p2 - p1;
    v = p3 - p1;
    n = u.cross(v);
    uu = u.dot(u);
    uv = u.dot(v);
    vv = v.dot(v);
    D = uv * uv - uu * vv;
}

Vec3f const& Triangle::operator[] (int index) const {
	if(index == 0)
		return p1;
	if(index == 1)
		return p2;
	if(index == 2)
		return p3;
	return this->mesh->vertices[this->mesh->faces[this->faceID * 3 + index]];
}

Vec3f Triangle::getCentroid() const
{
    return (p1 + p2 + p3) / 3.0f;
}

Vec3f Triangle::getAABBMin() const
{
    Vec3f min;
    min[0] = std::min(p1[0], std::min(p2[0], p3[0]));
    min[1] = std::min(p1[1], std::min(p2[1], p3[1]));
    min[2] = std::min(p1[2], std::min(p2[2], p3[2]));
    return min;
}

Vec3f Triangle::getAABBMax() const
{
    Vec3f max;
    max[0] = std::max(p1[0], std::max(p2[0], p3[0]));
    max[1] = std::max(p1[1], std::max(p2[1], p3[1]));
    max[2] = std::max(p1[2], std::max(p2[2], p3[2]));
    return max;
}

AABB Triangle::getAABB() const
{
    return AABB(this->getAABBMin(), this->getAABBMax());
}

Vec3f Triangle::getNormalVector() const
{
    Vec3f n = (p2 - p1).cross(p3 - p1);
    return n / n.length();
}

/*
 * Src: http://geomalgorithms.com/a06-_intersect-2.html
 */
bool Triangle::intersect(Ray const& ray, Intersection* intersection) const
{

	float const EPSILON2 = 0.000001;

    Vec3f w0 = ray.position - p1;
    float a = -n.dot(w0);
    float b = n.dot(ray.direction);

    // Check if ray is parallel to triangle plane.
    if (fabs(b) < EPSILON2)
        return false;

    // get intersect point of ray with triangle plane
    float r = a / b;
    if (r < 0.0) // ray goes away from triangle
        return false;

    // intersect point of ray and plane
    Vec3f intersectionPoint = ray.position + (ray.direction * r);

    // is I inside T?
    Vec3f w = intersectionPoint - p1;
    float wu = w.dot(u);
    float wv = w.dot(v);

    // get and test parametric coords
    float s = (uv * wv - vv * wu) / D;
    if (s < -0.00001f || s > 1.00001)         // I is outside T
        return false;
    float t = (uv * wu - uu * wv) / D;
    if (t < -0.00001f || (s + t) > 1.00001)  // I is outside T
        return false;

    // Intersection looks good. Fill result.
    float const distance = (intersectionPoint - ray.position).length();
    if (intersection->distance > distance)
    {
        intersection->mesh = this->mesh;
        intersection->faceID = this->faceID;
        intersection->bary = Vec3f(1.0f - s - t, s, t);
        intersection->position = intersectionPoint;
        intersection->distance = distance;
    }

    return true;
}
