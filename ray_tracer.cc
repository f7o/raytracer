#include <iostream>
#include <cmath>
#include <omp.h>
#include <sstream>
#include <fstream>
#include <cerrno>
#include <stdexcept>


#include "ray_tracer.h"
#include "hemisphere_sampler.h"
#include "timer.h"
#include "bvh_node.h"

RayTracer::RayTracer(RayTracerOptions const& opts) {
	this->opts = opts;
}

void RayTracer::trace(BVH const& scene, std::vector<unsigned char>* image) {
	//this->opencl_create(scene, image);

	//return;
	// Setup viewing parameters.
	int const height = this->opts.height;
	int const width = this->opts.width;
	float const focalLength = this->opts.focalLength;
	Vec3f cameraPosition(0.0f, 0.0f, 2.0f);
	float const ax = focalLength * (width > height ? width : height);
	float const ay = focalLength * (width > height ? width : height);

	// Generate image.
	image->resize(width * height, 0);

	std::cout << std::endl;
	// Split pixels into numThreads horizontal tiles
	if(this->opts.nThreads != 0)
		omp_set_num_threads(this->opts.nThreads);
	int threads= omp_get_max_threads();
	int tileRangeY = height / threads;
#pragma omp parallel for
	for(int tile=0; tile < threads; ++tile) {
		int tileStartY= tile*tileRangeY;
		std::stringstream ss;
		ss << ("Thread started, id: ") << tile <<
				", Associated pixel rows: " << tileStartY << " - " << tileStartY+tileRangeY-1
				<< std::endl;
		std::cout << ss.str();
		// Loop over every pixel of this tile
		for (int y = tileStartY; y < tileStartY+tileRangeY; ++y) {
			//std::cout << "\rProcessing " << (100 * y / height) << "%..." << std::flush;
			for (int x = 0; x < width; ++x) {
				int const index = x + y * width;
				// *** Main ray trace from camera position ***
				Ray ray;
				ray.position = cameraPosition;
				ray.direction = Vec3f((x + 0.5f) / ax - width / (2.0f * ax),
						(y + 0.5f) / ay - height / (2.0f * ay), 1.0f);
				ray.direction[1] *= -1.0f;  // Invert y-axis for image.
				ray.direction[2] *= -1.0f;  // Look along negative z-axis.
				ray.direction = ray.direction.normalized();
				ray.div= Vec3f(1.0f/ray.direction[0], 1.0f/ray.direction[1], 1.0f/ray.direction[2]);
				// Compute main intersection through the center of the pixel
				Intersection mainIntersection;
				mainIntersection.distance= 1e6;
				scene.intersect(ray, &mainIntersection, 1e6);
				// Check if intersection occurred
				if (mainIntersection.distance == 0.0f)
					continue;

				// *** Compute shading value for this pixel ***
				float value = 0.0f;
				if(this->opts.shading) {
					// Supersampling (1 is equal to no super sampling, hence, just one ray with its shading)
					float n = static_cast<float>(this->opts.nSuperSamples);
					for(int s=0; s < n; s++) {
						// Compute direction
						float offset= (s + 0.5f)/n;
						Ray ray;
						ray.position = cameraPosition;
						ray.direction= Vec3f(
								(x + offset)/ax - width/(2.0f * ax),
								(y + offset)/ay - height/(2.0f * ay),
								1.0f);
						ray.direction[1] *= -1.0f;  // Invert y-axis for image.
						ray.direction[2] *= -1.0f;  // Look along negative z-axis.
						ray.direction = ray.direction.normalized();
						ray.div= Vec3f(1.0f/ray.direction[0], 1.0f/ray.direction[1], 1.0f/ray.direction[2]);
						// Compute intersection
						Intersection intersection;
						intersection.distance= 1e6;
						//ray.direction = ray.direction.normalized();
						scene.intersect(ray, &intersection, 1e6);

						// Compute and add shading value of this spuer sample ray
						if (intersection.distance == 0.0f)
							continue;
						Vec3f normal = this->getSmoothNormal(intersection);

						float shading =	this->shading(ray, normal);
						value += shading;
					}
					value/= n; // Average shading contributions from super sampling

				} else {
					value= 1.0f; // No shading -> all white
				}
				// Add ambient occlusion factor to the shading value
				if (this->opts.ambientOcclusion) {
					Vec3f normal = this->getSmoothNormal(mainIntersection);
					value *= this->ambientOcclusion(scene, mainIntersection.position, normal, index);
				}
				// Set final shading for this pixel
				image->at(index) = static_cast<unsigned char>(value * 255.0f);
			}
		}
	}
}

float RayTracer::shading(Ray const& ray, Vec3f const& normal) {
	return std::min(1.0f, std::abs(normal.dot(ray.direction)));
}

float RayTracer::ambientOcclusion(BVH const& scene, Vec3f const& point,
		Vec3f const& normal, int const pID) {
	HemisphereSampler hSampler = HemisphereSampler(normal, pID);
	unsigned int aoCounter = 0;
	Vec3f interPoint = point + normal*1e-6;
	// Cast ao rays from intersection point into different directions
	for (int i = 0; i < this->opts.aoNumSamples; ++i) {
		Intersection intersection;
		intersection.distance= 1e6;
		Ray ray;
		ray.position = interPoint;
		ray.direction = hSampler.sample();
		ray.div= Vec3f(1.0f/ray.direction[0], 1.0f/ray.direction[1], 1.0f/ray.direction[2]);
		scene.intersect(ray, &intersection, this->opts.aoMaxDistance);
		if (intersection.distance == 0.0f)
			continue;

		if (intersection.distance <= this->opts.aoMaxDistance)
			aoCounter++;

	}
	return 1.0f
			- (static_cast<float>(aoCounter)
					/ static_cast<float>(this->opts.aoNumSamples));
}

Vec3f RayTracer::getSmoothNormal(Intersection const& intersection) {
	unsigned int v0 = intersection.mesh->faces[intersection.faceID * 3 + 0];
	unsigned int v1 = intersection.mesh->faces[intersection.faceID * 3 + 1];
	unsigned int v2 = intersection.mesh->faces[intersection.faceID * 3 + 2];
	return (intersection.mesh->vnormals[v0] * intersection.bary[0]
			+ intersection.mesh->vnormals[v1] * intersection.bary[1]
			+ intersection.mesh->vnormals[v2] * intersection.bary[2]).normalized();
}

Vec3f RayTracer::getFlatNormal(Intersection const& intersection) {
	unsigned int v0id = intersection.mesh->faces[intersection.faceID * 3 + 0];
	unsigned int v1id = intersection.mesh->faces[intersection.faceID * 3 + 1];
	unsigned int v2id = intersection.mesh->faces[intersection.faceID * 3 + 2];
	Vec3f const& v0 = intersection.mesh->vertices[v0id];
	Vec3f const& v1 = intersection.mesh->vertices[v1id];
	Vec3f const& v2 = intersection.mesh->vertices[v2id];
	return (v1 - v0).cross(v2 - v0).normalized();
}

