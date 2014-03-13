#include <iostream>
#include <cmath>
//#include <omp.h>
#include <sstream>
#include <fstream>
#include <cerrno>
#include <stdexcept>

//#include "cl_1.2.hpp"
#include "cl_ray_tracer.h"
#include "hemisphere_sampler.h"
#include "timer.h"
#include "bvh_node.h"

#define CHECK(v) if(v != CL_SUCCESS) {std::cout << "CL error " << get_error_string(v) << " at " << __FILE__ << ":" << __LINE__ << std::endl; std::exit(-1);}

// utility function to put out string errors
const char * get_error_string(cl_int err) {
	switch (err) {
	case 0:
		return "CL_SUCCESS";
	case -1:
		return "CL_DEVICE_NOT_FOUND";
	case -2:
		return "CL_DEVICE_NOT_AVAILABLE";
	case -3:
		return "CL_COMPILER_NOT_AVAILABLE";
	case -4:
		return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
	case -5:
		return "CL_OUT_OF_RESOURCES";
	case -6:
		return "CL_OUT_OF_HOST_MEMORY";
	case -7:
		return "CL_PROFILING_INFO_NOT_AVAILABLE";
	case -8:
		return "CL_MEM_COPY_OVERLAP";
	case -9:
		return "CL_IMAGE_FORMAT_MISMATCH";
	case -10:
		return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
	case -11:
		return "CL_BUILD_PROGRAM_FAILURE";
	case -12:
		return "CL_MAP_FAILURE";

	case -30:
		return "CL_INVALID_VALUE";
	case -31:
		return "CL_INVALID_DEVICE_TYPE";
	case -32:
		return "CL_INVALID_PLATFORM";
	case -33:
		return "CL_INVALID_DEVICE";
	case -34:
		return "CL_INVALID_CONTEXT";
	case -35:
		return "CL_INVALID_QUEUE_PROPERTIES";
	case -36:
		return "CL_INVALID_COMMAND_QUEUE";
	case -37:
		return "CL_INVALID_HOST_PTR";
	case -38:
		return "CL_INVALID_MEM_OBJECT";
	case -39:
		return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
	case -40:
		return "CL_INVALID_IMAGE_SIZE";
	case -41:
		return "CL_INVALID_SAMPLER";
	case -42:
		return "CL_INVALID_BINARY";
	case -43:
		return "CL_INVALID_BUILD_OPTIONS";
	case -44:
		return "CL_INVALID_PROGRAM";
	case -45:
		return "CL_INVALID_PROGRAM_EXECUTABLE";
	case -46:
		return "CL_INVALID_KERNEL_NAME";
	case -47:
		return "CL_INVALID_KERNEL_DEFINITION";
	case -48:
		return "CL_INVALID_KERNEL";
	case -49:
		return "CL_INVALID_ARG_INDEX";
	case -50:
		return "CL_INVALID_ARG_VALUE";
	case -51:
		return "CL_INVALID_ARG_SIZE";
	case -52:
		return "CL_INVALID_KERNEL_ARGS";
	case -53:
		return "CL_INVALID_WORK_DIMENSION";
	case -54:
		return "CL_INVALID_WORK_GROUP_SIZE";
	case -55:
		return "CL_INVALID_WORK_ITEM_SIZE";
	case -56:
		return "CL_INVALID_GLOBAL_OFFSET";
	case -57:
		return "CL_INVALID_EVENT_WAIT_LIST";
	case -58:
		return "CL_INVALID_EVENT";
	case -59:
		return "CL_INVALID_OPERATION";
	case -60:
		return "CL_INVALID_GL_OBJECT";
	case -61:
		return "CL_INVALID_BUFFER_SIZE";
	case -62:
		return "CL_INVALID_MIP_LEVEL";
	case -63:
		return "CL_INVALID_GLOBAL_WORK_SIZE";
	default:
		return "Unknown OpenCL error";
	}
}

CLRayTracer::CLRayTracer(RayTracerOptions const& opts) {
	this->_noDevice = false;
	this->_device = 0;
	this->opts = opts;
	this->opencl_init();
}

void CLRayTracer::read_file(std::string const& filename, std::string* data) {
	std::ifstream in(filename.c_str(), std::ios::binary);
	if (!in)
		throw std::runtime_error(::strerror(errno));

	in.seekg(0, std::ios::end);
	data->resize(in.tellg());
	in.seekg(0, std::ios::beg);
	in.read(&data->at(0), data->size());
	in.close();

}
int CLRayTracer::opencl_init() {
	std::cout << "Beginning OpenCL initialization..." << std::endl;
	// for error checking
	cl_int err;

	// query platforms
	std::vector<cl::Platform> platforms;
	CHECK(cl::Platform::get(&platforms));


	// get info
	std::cout << "OpenCL platforms found: " << platforms.size() << std::endl;
	int id = 0;
	for (std::vector<cl::Platform>::iterator it = platforms.begin();
			it != platforms.end(); ++it) {
		std::string vendor;
		std::string name;
		std::string extensions;

		it->getInfo(CL_PLATFORM_VENDOR, &vendor);
		it->getInfo(CL_PLATFORM_NAME, &name);
		it->getInfo(CL_PLATFORM_EXTENSIONS, &extensions);

		std::cout << "Platform " << id << ":" << std::endl;
		std::cout << "Vendor: " << vendor << std::endl;
		std::cout << "Name: " << name << std::endl;
		std::cout << "Extensions: " << extensions << std::endl;

		++id;
	}
	_useClPlatform = platforms[0];

	// select device
	std::vector<cl::Device> devices;
	CHECK(_useClPlatform.getDevices(CL_DEVICE_TYPE_ALL, &devices));
	if (devices.size() == 0) {
		this->_noDevice = true;
		std::cout << "No OpenCL devices found! " << std::endl;
		return 0;
	}
	this->_noDevice = false;
	std::cout << "OpenCL devices found: " << devices.size() << std::endl;
	id = 0;
	for (std::vector<cl::Device>::iterator it = devices.begin();
			it != devices.end(); ++it) {

		std::string name;

		it->getInfo(CL_DEVICE_NAME, &name);

		std::cout << "Device " << id << ": " << name << std::endl;

		++id;
	}
	_useClDevice = devices[_device];

	// only work with selected device
	devices.clear();
	devices.push_back(_useClDevice);

	// create context
	_useClContext = cl::Context(devices, NULL, NULL, NULL, &err);
	CHECK(err);

	// create single command queue for used device
	// remark: as of today, CUDA only supports one
	// CL command queue simultaneously on NVIDIA
	// devices
	_useClQueue = cl::CommandQueue(_useClContext, _useClDevice, 0, &err);
	CHECK(err);
	return 1;
}

int CLRayTracer::opencl_create(BVH const& scene, unsigned char* image,
		Mesh const& mesh, int size_scene_nodes) {

	/* compile kernels */

	try {
		compileKernel("ray_tracer.cl");
	} catch (std::exception & e) {
		std::cout << e.what() << std::endl;
		return 0;
	}
	ClockTimer timer;
	std::cout << "RayTracing with OpenCL ... " << std::endl;
	this->compute_rays();
	/* define matrix dimensions */

	struct cl_bvh_node* noteScene = new cl_bvh_node[size_scene_nodes];

	scene.getBVH(noteScene);

	this->compute_intersections(noteScene);

	if(this->opts.shading){

		this->compute_shading();

	}
	this->_useClQueue.enqueueReadBuffer(this->_imageBuffer, true, 0, this->opts.height*this->opts.width * sizeof(unsigned char), image);
	if(this->opts.ambientOcclusion) {
		this->compute_ao_rays();
		this->compute_ao_intersections(noteScene);
		this->compute_ao(noteScene);
		this->compute_ao_point_intersections(noteScene);
		this->compute_ao_factor(image);
		this->_useClQueue.enqueueReadBuffer(this->_imageAoBuffer, true, 0, this->opts.height*this->opts.width * sizeof(unsigned char), image);
	}

	//this->_useClQueue.enqueueReadBuffer(this->_imageBuffer, true, 0, this->opts.height*this->opts.width * sizeof(unsigned char), image);

	/* finish queue */
	this->_useClQueue.finish();

	std::cout << " took " << timer.get_elapsed() << "ms." << std::endl;

	return 1;
}

int CLRayTracer::compileKernel(std::string const& filename) {
	std::cout << "Compile OpenCL Kernel ..."
			<< std::endl;

	std::ifstream kernelFile(filename.c_str());
	if (!kernelFile.is_open())
		throw("Cannot open source file for compiling!");

	// create option string for compiling program
		std::string options;
		std::stringstream s_options;

		s_options << "-D WIDTH=" << this->opts.width << " "
				  << "-D HEIGHT=" << this->opts.height << " "
				  << "-D SAMPLES=" << this->opts.nSuperSamples << " "
				  << "-D NODE_COUNT=" << this->opts.nodeCount << " "
				  << "-D AO_MAX_DISTANCE=" << this->opts.aoMaxDistance << " "
				  << "-D AO_SAMPLES=" << this->opts.aoNumSamples << " "
				  << "-D MAX_LEAF_TRIANGLES=" << 10;
		options = s_options.str();
		std::cout << "Compilation options: " << options << std::endl;


	std::string source;
	read_file(filename, &source);
	cl::Program::Sources sources(1,
			std::make_pair(source.c_str(), source.length()));
	this->_useClProgram = cl::Program(this->_useClContext, sources);

	std::vector<cl::Device> devices;
	devices.push_back(this->_useClDevice);
	if (this->_useClProgram.build(devices, options.c_str()) != 0)
		std::cout << "Error: "
				<< this->_useClProgram.getBuildInfo<CL_PROGRAM_BUILD_LOG>(
						this->_useClDevice) << std::endl;

	this->_useClKernelRays = cl::Kernel(this->_useClProgram, "compute_rays");
	this->_useClKernelIntersect = cl::Kernel(this->_useClProgram, "compute_intersections");
	this->_useClKernelShade = cl::Kernel(this->_useClProgram, "compute_shading");
	this->_useClKernelAoRays = cl::Kernel(this->_useClProgram, "compute_ao_rays");
	this->_useClKernelAoFaktor = cl::Kernel(this->_useClProgram, "compute_ao_factor");
	this->_useClKernelAoPointRays = cl::Kernel(this->_useClProgram, "compute_ao_point_rays");
	return 1;
}



int CLRayTracer::compute_rays() {
	int const height = this->opts.height;
	int const width = this->opts.width;
	int const nSamples = this->opts.nSuperSamples;
	float const focalLength = this->opts.focalLength;
	float const cameraPosition[3] = { 0.0f, 0.0f, 2.0f };

	float const ax = focalLength * (width > height ? width : height);
	float const ay = focalLength * (width > height ? width : height);

	this->_rayBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
			(height*width*nSamples)*48, NULL);

	this->_useClKernelRays.setArg(0, this->_rayBuffer);
	this->_useClKernelRays.setArg(1, cameraPosition[0]);
	this->_useClKernelRays.setArg(2, cameraPosition[1]);
	this->_useClKernelRays.setArg(3, cameraPosition[2]);
	this->_useClKernelRays.setArg(4, ax);
	this->_useClKernelRays.setArg(5, ay);


	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width, height);

	/* launch kernel */
	cl::Event event;
	std::cout << "Rays Computing start ..." << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelRays,
			cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "Ray Computing ende ... " << std::endl;
	/* finish queue */
	//this->_useClQueue.finish();
	return 1;
}
int CLRayTracer::compute_ao_rays() {
	int const height = this->opts.height;
	int const width = this->opts.width;

	float const focalLength = this->opts.focalLength;
	float const cameraPosition[3] = { 0.0f, 0.0f, 2.0f };

	float const ax = focalLength * (width > height ? width : height);
	float const ay = focalLength * (width > height ? width : height);

	this->_rayBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
			(height*width)*48, NULL);

	this->_useClKernelAoRays.setArg(0, this->_rayBuffer);
	this->_useClKernelAoRays.setArg(1, cameraPosition[0]);
	this->_useClKernelAoRays.setArg(2, cameraPosition[1]);
	this->_useClKernelAoRays.setArg(3, cameraPosition[2]);
	this->_useClKernelAoRays.setArg(4, ax);
	this->_useClKernelAoRays.setArg(5, ay);


	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width, height);


	cl::Event event;
	std::cout << "Rays AO Computing start ..." << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelAoRays,
			cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "Ray AO Computing ende ... " << std::endl;

	//this->_useClQueue.finish();
	return 1;
}

int CLRayTracer::compute_intersections(cl_bvh_node* cl_bvh) {
	int const height = this->opts.height;
		int const width = this->opts.width;
		int const nSamples = this->opts.nSuperSamples;
		int const rayCount = height * width * nSamples;

		this->_intersectionBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
					(height*width*nSamples)*80, NULL);
		cl::Buffer nodes = cl::Buffer(this->_useClContext,
					CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
					this->opts.nodeCount * sizeof(struct cl_bvh_node), cl_bvh);

		this->_useClKernelIntersect.setArg(0, this->_intersectionBuffer);
		this->_useClKernelIntersect.setArg(1, rayCount);
		this->_useClKernelIntersect.setArg(2, this->_rayBuffer);
		this->_useClKernelIntersect.setArg(3, nodes);

		cl::NDRange localWorkSize(8, 8);
		cl::NDRange globalWorkSize(rayCount);

		/* launch kernel */
		cl::Event event;
		std::cout << "inter start" << std::endl;
		this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelIntersect,
					cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
		event.wait();
		std::cout << "inter ende" << std::endl;
		/* finish queue */
		//this->_useClQueue.finish();

		return 1;
}

int CLRayTracer::compute_ao_intersections(cl_bvh_node* cl_bvh) {
	int const height = this->opts.height;
	int const width = this->opts.width;
	int const nSamples = this->opts.nSuperSamples;

	this->_intersectionBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
				(height*width)*80, NULL);
	cl::Buffer nodes = cl::Buffer(this->_useClContext,
						CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
						this->opts.nodeCount * sizeof(struct cl_bvh_node), cl_bvh);
	this->_useClKernelIntersect.setArg(0, this->_intersectionBuffer);
	this->_useClKernelIntersect.setArg(2, this->_rayBuffer);
	this->_useClKernelIntersect.setArg(3, nodes);

	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width*height);


	cl::Event event;
	std::cout << "AO Intersection start ..." << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelIntersect,
				cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "AO Intersection ende ..." << std::endl;

	//this->_useClQueue.finish();

	return 1;
}
int CLRayTracer::compute_ao_point_intersections(cl_bvh_node* cl_bvh) {
	int const height = this->opts.height;
	int const width = this->opts.width;


	this->_intersectionBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
				(height*width*this->opts.aoNumSamples)*80, NULL);
	cl::Buffer nodes = cl::Buffer(this->_useClContext,
							CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
							this->opts.nodeCount * sizeof(struct cl_bvh_node), cl_bvh);

	this->_useClKernelIntersect.setArg(0, this->_intersectionBuffer);
	this->_useClKernelIntersect.setArg(2, this->_aoRayBuffer);
	this->_useClKernelIntersect.setArg(3, nodes);

	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width*height*this->opts.aoNumSamples);


	cl::Event event;
	std::cout << "AO Intersection start ..." << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelIntersect,
				cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "AO Intersection ende ..." << std::endl;



	return 1;
}
void CLRayTracer::compute_shading() {
	int const height = this->opts.height;
	int const width = this->opts.width;

	this->_imageBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
					height*width * sizeof(unsigned char), NULL);

	this->_useClKernelShade.setArg(0, this->_intersectionBuffer);
	this->_useClKernelShade.setArg(1, this->_rayBuffer);
	this->_useClKernelShade.setArg(2, this->_imageBuffer);

	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width, height);

	/* launch kernel */
	cl::Event event;
	std::cout << "Shading start ..." << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelShade,
				cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "Shading end ..." << std::endl;




}

void CLRayTracer::compute_ao(cl_bvh_node* cl_bvh) {
	float* aoRandoms = new float[this->opts.aoNumSamples*3];
	int const height = this->opts.height;
	int const width = this->opts.width;
		  	HemisphereSampler hSampler = HemisphereSampler(Vec3f(0.0f,0.0f,0.0f), 234245);

		   	for(int i =0; i< this->opts.aoNumSamples;i++) {

		   		Vec3f tmp = hSampler.randomFloat();

		   		aoRandoms[i*3+0] = tmp[0];

		   		aoRandoms[i*3+1] = tmp[1];

		   		aoRandoms[i*3+2] = tmp[2];
		}
		   	this->_aoRandomBuffer = cl::Buffer(this->_useClContext,
		   					CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		   					this->opts.aoNumSamples * sizeof(float), aoRandoms);
		   	this->_aoRayBuffer = cl::Buffer(this->_useClContext,
		   					CL_MEM_WRITE_ONLY,
		   			   		width*height*this->opts.aoNumSamples * 48, NULL);
		   	cl::Buffer nodes = cl::Buffer(this->_useClContext,
		   								CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		   								this->opts.nodeCount * sizeof(struct cl_bvh_node), cl_bvh);

		   	this->_useClKernelAoPointRays.setArg(0, this->_aoRayBuffer);
		   	this->_useClKernelAoPointRays.setArg(1, nodes);
		   	this->_useClKernelAoPointRays.setArg(2, this->_rayBuffer);
		   	this->_useClKernelAoPointRays.setArg(3, this->_intersectionBuffer);
		   	this->_useClKernelAoPointRays.setArg(4, this->_aoRandomBuffer);


		   	cl::NDRange localWorkSize(8, 8);
		   	cl::NDRange globalWorkSize(width, height);


		   	cl::Event event;
		   	std::cout << "ao Rays start" << std::endl;
		   	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelAoPointRays,
		   					cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);

		   	event.wait();
		   	std::cout << "ao Rays ende" << std::endl;
}

void CLRayTracer::compute_ao_factor(unsigned char* image) {
	int const height = this->opts.height;
	int const width = this->opts.width;

	this->_imageAoBuffer = cl::Buffer(this->_useClContext, CL_MEM_WRITE_ONLY,
						height*width * sizeof(unsigned char), NULL);
	this->_imageAoBuffer = cl::Buffer(this->_useClContext,
					CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
					width*height * sizeof(unsigned char), image);
	this->_useClKernelAoFaktor.setArg(0, this->_intersectionBuffer);
	this->_useClKernelAoFaktor.setArg(1, this->_imageBuffer);
	this->_useClKernelAoFaktor.setArg(2, this->_imageAoBuffer);

	cl::NDRange localWorkSize(8, 8);
	cl::NDRange globalWorkSize(width, height);

	cl::Event event;
	std::cout << "ao faktor compute start" << std::endl;
	this->_useClQueue.enqueueNDRangeKernel(this->_useClKernelAoFaktor,
	cl::NDRange(), globalWorkSize, localWorkSize, NULL, &event);
	event.wait();
	std::cout << "ao faktor compute ende" << std::endl;
}
