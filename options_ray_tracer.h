/*
 * options_ray_tracer.h
 *
 *  Created on: Jan 22, 2014
 *      Author: flo
 */

#ifndef OPTIONS_RAY_TRACER_H_
#define OPTIONS_RAY_TRACER_H_

// Structure with basic options for raytracer
    // - width : image witdh
    // - focalLength      : focal length of virtual camera
    // - nSuperSamples  : Number of Samples !per pixel!
    //                      Runtime basically explodes if you turn this too high
    // - smoothShading    : switch to turn on/off smooth shading
    // - ambientOcclusion : switch to turn on/off ambient occlusion
    // - aoMaxDistance    : Ambient Occlusion distance
    //                      Should be about 10% of the max scene dimension
    // - aoNumSamples     : Number of samples for each ambient occlusion
    //                      evaluation
    struct RayTracerOptions
    {
        int width;
        int height;
        float focalLength;
        int nSuperSamples; // // a number of n^2
        bool shading;
        bool ambientOcclusion;
        float aoMaxDistance;
        int aoNumSamples;
        int nThreads; // a number of n^2
        bool openCL;
        int nodeCount;
    };


#endif /* OPTIONS_RAY_TRACER_H_ */
