#ifndef CL_RAY_TRACER_H
#define CL_RAY_TRACER_H

#include <vector>

#include "mesh.h"
#include "ray.h"
#include "triangle.h"
#include "bvh.h"
//#include "CL/cl.h"
#include "cl_1.1.hpp"
#include "cl_bvh.h"
#include "options_ray_tracer.h"

class CLRayTracer
{
protected:
	// GPU boilerplate stuff
			int _device;
			bool _noDevice;
			cl::Platform _useClPlatform;
			cl::Device _useClDevice;
			cl::Context _useClContext;
			cl::Program _useClProgram;
			cl::CommandQueue _useClQueue;
			cl::Kernel _useClKernelRays;
			cl::Kernel _useClKernelIntersect;
			cl::Kernel _useClKernelShade;
			cl::Kernel _useClKernelAoRays;
			cl::Kernel _useClKernelAoPointRays;
			cl::Kernel _useClKernelAoFaktor;
			cl::Buffer _imageBuffer;
			cl::Buffer _imageAoBuffer;
			cl::Buffer _rayBuffer;
			cl::Buffer _intersectionBuffer;
			cl::Buffer _aoRayBuffer;
			cl::Buffer _aoRandomBuffer;
			cl::Buffer _nodeBuffer;


public:





public:
    // Constructor with options argument
    CLRayTracer(RayTracerOptions const& opts);
    
    bool noDevice();
    int opencl_create(BVH const& scene, unsigned char* image, Mesh const& mesh, int size_scene_nodes);



private:
    int
    opencl_init();
    int compileKernel(std::string const& filename);
    void read_file (std::string const& filename, std::string* data);
    int compute_rays();
    int compute_ao_rays();
    int compute_intersections(cl_bvh_node* cl_bvh);
    int compute_ao_intersections(cl_bvh_node* cl_bvhs);
    int compute_ao_point_intersections(cl_bvh_node* cl_bvhs);
    void compute_ao(cl_bvh_node* cl_bvhs);
    void compute_shading();
    void compute_ao_factor(unsigned char* image);


private:
    RayTracerOptions opts;
};

inline bool CLRayTracer::noDevice() {return this->_noDevice;}

#endif // CL_RAY_TRACER_H

