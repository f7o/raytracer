//#include "ray_tracer.h"
//#include "hemisphere_sampler.h"
// 

//float3 getNormal(struct Intersection *inter, __global float* scene) ;
//float shading(struct Ray *ray, float3 *normal);
//bool intersectTriangle(__global float* scene, struct Ray* ray, struct Intersection* intersection);
#define MAX_LEAF_TRIANGLES 10
struct Ray {
	float3 position;
	float3 direction;
};

struct Triangle {
	int id;
	float3 a, b, c;
};
struct AABB {
	float3 min;
	float3 max;
};
struct cl_triangle {
    	int id;
    	float a[3], b[3], c[3];
    };
struct cl_aabb {
    	float min[3];
    	float max[3];
    };
struct cl_bvh_node {
    	int leftIndex;
    	int rightIndex;
    	int tri_counter;
    	float min[3];
    	float max[3];
    	float objects[MAX_LEAF_TRIANGLES*9];
    	//struct cl_triangle objects[10];
    };
struct Intersection {
	float3 bary;
	float3 position;
	float distance;
	struct Triangle *obj;
};

float3 my_cross(float3 a, float3 b) {
	return (float3)(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}
float3 getNormal(struct Intersection *inter) {
	
	float3 v0 = inter->obj->a;
	float3 v1 = inter->obj->b;
	float3 v2 = inter->obj->c;

return normalize(cross(v1-v0,v2-v0));

}
float3 getSmoothNormal(struct Intersection *inter, __global float* scene) {
	float3 v0 = inter->obj->a;
	float3 v1 = inter->obj->b;
	float3 v2 = inter->obj->c;
	/*unsigned int v0 = intersection.mesh->faces[intersection.faceID * 3 + 0];
	unsigned int v1 = intersection.mesh->faces[intersection.faceID * 3 + 1];
	unsigned int v2 = intersection.mesh->faces[intersection.faceID * 3 + 2];*/
	return (float3)(0.0f,0.0f,0.0f);
	//return (intersection.mesh->vnormals[v0] * intersection.bary[0]
	//		+ intersection.mesh->vnormals[v1] * intersection.bary[1]
	//		+ intersection.mesh->vnormals[v2] * intersection.bary[2]).normalized();
}

float shading(struct Ray *ray, float3 *normal) {
	float doted = 0.5f;//
	doted = dot(*normal, ray->direction);
	//doted = normal.x*(ray->direction).x+
	//		normal.y*ray->direction.y+
	//		normal.z*ray->direction.z;
	return fmin(1.0f, fabs(doted));
	
}



bool intersectAABB(struct Ray* ray, struct Intersection* intersection, struct AABB* aabb, float* dis)
{
	// Smits ray-box intersection test using slabs
	    // http://www.cs.utah.edu/~awilliam/box/box.pdf
	    float tmin, tmax, tymin, tymax, tzmin, tzmax;
	    float div = 1.0f / ray->direction.x;
	    if (div >= 0)
	    {
	        tmin = (aabb->min.x - ray->position.x) * div;
	        tmax = (aabb->max.x - ray->position.x) * div;
	    }
	    else
	    {
	        tmin = (aabb->max.x - ray->position.x) * div;
	        tmax = (aabb->min.x - ray->position.x) * div;
	    }

	    div = 1 / ray->direction.y;
	    if (div >= 0)
	    {
	        tymin = (aabb->min.y - ray->position.y) * div;
	        tymax = (aabb->max.y - ray->position.y) * div;
	    }
	    else
	    {
	        tymin = (aabb->max.y  - ray->position.y) * div;
	        tymax = (aabb->min.y - ray->position.y) * div;
	    }

	    if (tmin > tymax || tymin > tmax)
	        return false;

	    tmin = max(tmin, tymin);
	    tmax = min(tmax, tymax);

	    div = 1 / ray->direction.z;
	    if (div >= 0)
	    {
	        tzmin = (aabb->min.z - ray->position.z) * div;
	        tzmax = (aabb->max.z - ray->position.z) * div;
	    }
	    else
	    {
	        tzmin = (aabb->max.z - ray->position.z) * div;
	        tzmax = (aabb->min.z - ray->position.z) * div;
	    }

	    if (tmin > tzmax || tzmin > tmax)
	        return false;

	    tmin = max(tmin, tzmin);
	    tmax = min(tmax, tzmax);
	    *dis = fabs(tmin);
	    return tmin < 100000.0f && tmax > 0.0f;

}
bool intersectTriangle(struct Ray* ray, struct Intersection* intersection, struct Triangle* tri)
{
	
	
		 float const EPSILON2 = 0.000001;
		 //printf("%d+%i\n",j,size_scene);
		    // get triangle edge vectors and plane normal
		 	float3 a = tri->a;
		 	float3 b = tri->b;
		 	float3 c = tri->c;
		 	
		    float3 u = b - a;
		    float3 v = c - a;
		    float3 n = cross(u,v);              // cross product

		    float3 w0 = ray->position - a;
		    float a0 = (-1.0f)*dot(n,w0);
		    float b0 = dot(n,ray->direction);
		    
		    // Check if ray is parallel to triangle plane.
		    if (fabs(b0) < EPSILON2)
		    	//continue;
		        return false;
		    
		    // get intersect pointmy of ray with triangle plane
		    float r = a0 / b0;
		    //printf("%f:rrrr\n",r);
		    if (r < 0.0) // ray goes away from triangle
		    	//continue;
		        return false;

		    // intersect point of ray and plane
		    float3 intersectionPoint = ray->position + (ray->direction * r);

		    // is I inside T?
		    float uu = dot(u,u);
		    float uv = dot(u,v);
		    float vv = dot(v,v);
		    float3 w = intersectionPoint - a;
		    float wu = dot(w,u);
		    float wv = dot(w,v);
		    float D = uv * uv - uu * vv;

		    // get and test parametric coords
		    float s = (uv * wv - vv * wu) / D;
		    if (s < -0.00001f || s > 1.00001)         // I is outside T
		        return false;
		    	//continue;
		    float t = (uv * wu - uu * wv) / D;
		    if (t < -0.00001f || (s + t) > 1.00001)  // I is outside T
		        return false;
		    	//continue;
	
		    // Intersection looks good. Fill result.
		    float inter_div = distance(intersectionPoint,ray->position);
			
		    if (intersection->distance > inter_div)
		    {
		    	//printf("success");
		        //intersection->mesh = this->mesh;
		        //intersection->faceID = this->faceID;
		        intersection->bary = (float3)(1.0f - s - t, s, t);
		        intersection->position = intersectionPoint;
		        intersection->distance = inter_div;
		        struct Triangle tri_inter;
		        tri_inter.a = a;
		        tri_inter.b = b;
		        tri_inter.c = c;
		        tri_inter.id = tri->id;
		        intersection->obj = &tri_inter;
		        //if(get_global_id(0)==2030) {
		        		//printf("dist: %f\n",intersection->distance);
		        //		//float3 a  = intersection->obj->a;
		          								 	
		    }	    
	
	
	return true;

}
bool pointInside(struct AABB aabb, float3 point) {
	
		if (point.x > aabb.max.x || point.x < aabb.min.x)
			return false;
		if (point.y > aabb.max.y || point.y < aabb.min.y)
					return false;
		if (point.z > aabb.max.z || point.z < aabb.min.z)
					return false;
	return true;
}

bool intersect(struct Ray* ray, struct Intersection* intersection, __global struct cl_bvh_node* nodes, int size_scene_nodes, float maxDistance)
{
	printf("%d\n",size_scene_nodes);
	const int p = size_scene_nodes;
	int fringe[515151];
	int f = 1;
	int currentIndex = 0;
	fringe[0] = currentIndex;
	while(currentIndex < f) {
		struct cl_bvh_node act_node = nodes[fringe[currentIndex]];
		currentIndex++;
		
		struct AABB aabb;
		aabb.min = (float3)(act_node.min[0],act_node.min[1],act_node.min[2]);
		aabb.max = (float3)(act_node.max[0],act_node.max[1],act_node.max[2]);
		float distance= 0.0f;
		//if(!intersectAABB(ray, intersection, &aabb, &distance))
		//	continue;
		
			if(!pointInside(aabb, ray->position)) {
				float distance= 0.0f;
				if(!intersectAABB(ray, intersection, &aabb, &distance))
					continue;
				if(distance > maxDistance)
					continue;
				if(distance > intersection->distance) 
					continue;
				
			}
		
		if (act_node.tri_counter != 0)
		    {
		        bool success = false;
		        //printf("%d\n",act_node.tri_counter);
		        for (int i = 0; i < act_node.tri_counter; i++) {
		        	struct Triangle tri;
		        	tri.a = (float3)(act_node.objects[i*9+0],act_node.objects[i*9+1],act_node.objects[i*9+2]);
		        	tri.b = (float3)(act_node.objects[i*9+3],act_node.objects[i*9+4],act_node.objects[i*9+5]);
		        	tri.c = (float3)(act_node.objects[i*9+6],act_node.objects[i*9+7],act_node.objects[i*9+8]);
		        	
		        	
		        	success |= intersectTriangle(ray, intersection, &tri);
		        }
		        continue;
		    }
		
		fringe[f] = act_node.leftIndex;
		//printf("%d\n",f);
		f=f+1;
		
		fringe[f] = act_node.rightIndex;
		f=f+1;
	}
	return true;
}



float ambientOcclusion(float3  point,
		float3 normal, int size_scene_node,__global struct cl_bvh_node*  nodes,int aoSamples,float aoMaxDistance, __global float* aoRandoms) {
	
	float3 basisY = normalize(normal);
	    float3 h = basisY;
	    if (fabs(h.x)<=fabs(h.y) && fabs(h.x)<=fabs(h.z))
	        h.x= 1.0;
	    else if (fabs(h.y)<=fabs(h.x) && fabs(h.y)<=fabs(h.z))
	        h.y= 1.0;
	    else
	        h.z= 1.0;
	    
	    float3 basisX = cross(h,basisY);
	    basisX = normalize(basisX);
	    float3 basisZ = cross(basisX, basisY);
	    basisZ = normalize(basisZ);
	    
	    //this->randGenrator.initialize(536870923u * pID);
	    
	    /* use better random generator */
	    
	

	unsigned int aoCounter = 0;
	float3 interPoint = point + (float3)(1e-6 * normal.x,1e-6 * normal.y,1e-6 * normal.z) ;
	// Cast ao rays from intersection point into different directions
	
	for (int i = 0; i < aoSamples; ++i) {
		struct Intersection intersection;
		intersection.distance= 1e6;
		struct Ray ray;
		ray.position = interPoint;
		
		// Generate "random" number
		//float xi1 = (float)(i)/(float)(aoSamples)*aoMaxDistance;
		//float xi2 = (float)(aoSamples)/(float)(i)*aoMaxDistance;
	    
	    //float theta = acos(sqrt(1.0-xi1));
	    //float phi = 2.0 * 2.71828182846 * xi2;
	    
	    float xs = aoRandoms[i*3+0];//sin(theta) * cos(phi);
	    float ys = aoRandoms[i*3+1];//cos(theta);
	    float zs = aoRandoms[i*3+2];//sin(theta) * sin(phi);
	    
	    float3 direction = basisX * xs + basisY * ys + basisZ * zs;
			        
			        
		
		
		ray.direction = normalize(direction);
		
		intersect(&ray, &intersection, nodes, size_scene_node, aoMaxDistance);
		if (intersection.distance == 1e6)
			continue;

		if (intersection.distance <= aoMaxDistance)
			aoCounter++;

	}
	//rintf("%i\n",aoCounter);
	return 1.0f
			- ((float)(aoCounter)
					/ (float)(aoSamples));
}
__kernel void tracer( __global unsigned char *image,
	uint width, uint height,
	float ax, float ay, int size_scene, int size_scene_nodes, __global struct cl_bvh_node* nodes, int nSamples, int aoSamples, float aoMaxDistance, int shade, __global float* aoRandoms)
{	


	
float x = (float)(get_global_id(0) % width);
float y = (float)(get_global_id(0) / width);

struct Intersection mainIntersection;
// *** Main ray trace from camera position if AO ***
if(aoSamples != 0) {
struct Ray ray;
ray.position = (float3)(0.0f, 0.0f, 2.0f);
ray.direction = (float3)((x + 0.5f) / ax - (float)(width) / (2.0f * ax),
		(-1.0f)*((y + 0.5f) / ay - (float)(height) / (2.0f * ay)), -1.0f);
//printf("%f+%f+%f\n",ray.direction.x,ray.direction.y,ray.direction.z);
//ray.direction[1] *= -1.0f;// Invert y-axis for image.
//ray.direction[2] *= -1.0f;// Look along negative z-axis.
ray.direction = normalize(ray.direction);
//ray.div= Vec3f(1.0f/ray.direction[0], 1.0f/ray.direction[1], 1.0f/ray.direction[2]);
// Compute main intersection through the center of the pixel

mainIntersection.distance = 1e6;
intersect(&ray, &mainIntersection, nodes, size_scene_nodes, 1e6);
//printf("%f",mainIntersection.bary.x);
}
 
int n = nSamples;
int n2 = nSamples*nSamples;
// *** Compute shading value for this pixel ***
float value = 0.0f;
if(shade == 1) {
// Supersampling (1 is equal to no super sampling, hence, just one ray with its shading)            float n = static_cast<float>(this->opts.nSuperSamples);
//printf("%f\n",value);
for(int s=0; s < n; s++) {
	// Compute direction
	float offset= (s + 0.5f)/n;
	struct Ray ray;
	ray.position = (float3)(0.0f, 0.0f, 2.0f);
	ray.direction= (float3)(
			(x + offset)/ax - (float)(width)/(2.0f * ax),
			((y + offset)/ay - (float)(height)/(2.0f * ay))*(-1.0f),
			-1.0f
			);
	//ray.div= Vec3f(1.0f/ray.direction[0], 1.0f/ray.direction[1], 1.0f/ray.direction[2]);
	// Compute intersection
	struct Intersection intersection;
	intersection.distance= 1e6;
	//ray.direction = ray.direction.normalized();
	intersect(&ray, &intersection, nodes, size_scene_nodes, 1e6);
	//printf("%d",sizeof(scene));
	// Compute and add shading value of this spuer sample ray
	//printf("%f\n",value);
	if (intersection.distance == 1e6)
		continue;
	float3 normal = getNormal(&intersection);
	//printf("%f+%f+%f\n",normal.x,normal.y,normal.z);
	//printf("%i\n", shade);
	
	value += shading(&ray, &normal);
		
	//printf("%f\n",value);
	///value += shade;
}
value/= (float)(n); // Average shading contributions from super sampling
//image[10] = (unsigned char)(0.2f);
//printf("testee: %f\n", value*255.0f);
// Set final shading for this pixel
} else {
	value = 1.0f;
}
float ao = 1.0f;
if(aoSamples != 0) {

	float3 normal = getNormal(&mainIntersection);	

	ao = ambientOcclusion(mainIntersection.position,normal, size_scene_nodes,nodes,aoSamples,aoMaxDistance, aoRandoms);
}
value *= ao;
image[get_global_id(0)] = (unsigned char)(value * 255.0f);

 
}

