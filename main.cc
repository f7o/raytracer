#include <fstream>
#include <vector>
#include <limits>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "timer.h"
#include "mesh.h"
#include "vec3.h"
#include "bvh.h"
#include "triangle.h"
#include "visual.h"
#include "ray_tracer.h"
#include "cl_ray_tracer.h"
#include "options_ray_tracer.h"

// Main routine.
int main (int argc, char** argv)
{

	// Check arguments.
    if (argc != 3)
    {
        std::cerr << "Syntax: " << argv[0] << " <IN_MESH.OFF> <OUT_IMAGE.PGM>" << std::endl;
        return 1;
    }
    std::string inMesh(argv[1]);
    std::string outImage(argv[2]);

    // Read input mesh.
    std::cout << "Reading input mesh..." << std::endl;
    Mesh mesh;
    load_off_mesh(inMesh, &mesh);
    compute_vertex_normals(&mesh);
    std::cout << "  Vertices: " << mesh.vertices.size() << std::endl;
    std::cout << "  Triangles: " << (mesh.faces.size() / 3) << std::endl;

    // Build BVH.
    std::cout << "Building BVH..." << std::flush;
    BVH bvh;
    {
        ClockTimer timer;
        bvh.buildBVH(mesh);
        std::cout << " took " << timer.get_elapsed() << "ms." << std::endl;
    }
    int node_size = bvh.printStatistics(std::cout);

    // Render options for output image.
    RayTracerOptions rtOptions;
    rtOptions.width = 600;
    rtOptions.height = 600;
    rtOptions.focalLength = 1.0f;
    rtOptions.nSuperSamples = 4;
    rtOptions.shading = true;
    rtOptions.ambientOcclusion = true;
    rtOptions.aoMaxDistance = 0.2f;
    rtOptions.aoNumSamples = 1000;
    rtOptions.nThreads = 16; // (0 is auto selection on runtime)
    //*************switch between openCL and CPU calculating ******************//
    rtOptions.openCL = true;
    rtOptions.nodeCount = node_size;

    std::cout << "Rendering image..." << std::flush;
    std::vector<unsigned char> image;

    bool no_gpu = false;
    unsigned char* image2;
    if(rtOptions.openCL) {
    	image2 = new unsigned char[rtOptions.width*rtOptions.height];
    	CLRayTracer clTracer(rtOptions);
    	if(clTracer.noDevice()) {
    		no_gpu = true;
    	} else {
    		clTracer.opencl_create(bvh, image2, mesh, node_size);
    	}
    } else {
    	ClockTimer timer;
    	RayTracer tracer(rtOptions);
    	tracer.trace(bvh, &image);
    	std::cout << " took " << timer.get_elapsed() << "ms." << std::endl;
    }

    // Write output image.
    std::ofstream out(outImage.c_str());
    if (!out.good())
    {
        std::cerr << "Error opening output file!" << std::endl;
        return 1;
    }
    out << "P5 " << rtOptions.width << " " << rtOptions.height << " 255\n";
    if(no_gpu || !rtOptions.openCL) {
    	out.write(reinterpret_cast<char const*>(&image[0]), rtOptions.width * rtOptions.height);
    } else {
    	std::cout << "Writing OpenCl output..." << std::endl;
    	out.write((char*)(image2), rtOptions.width * rtOptions.height);
    }
    out.close();

    return 0;
}
