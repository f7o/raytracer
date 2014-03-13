/*
 * cl_bvh.h
 *
 *  Created on: Jan 21, 2014
 *      Author: flo
 */

#ifndef CL_BVH_H_
#define CL_BVH_H_
#define MAX_LEAF_TRIANGLES 10

typedef struct cl_triangle {
    	int id;
    	float a[3], b[3], c[3];
    }cl_triangle;
typedef struct cl_aabb {
    	float min[3];
    	float max[3];
    }cl_aabb;
typedef struct cl_bvh_node {
    	int leftIndex;
    	int rightIndex;
    	int tri_counter;
    	float min[3];
    	float max[3];
    	float objects[MAX_LEAF_TRIANGLES*9];
    	//struct cl_triangle objects[MAX_LEAF_TRIANGLES];

    }cl_bvh_node;



#endif /* CL_BVH_H_ */
