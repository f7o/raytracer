
struct Ray {
	float3 position;
	float3 direction;
	float3 div;
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
    };
struct Intersection {
	float3 bary;
	float3 position;
	float distance;
	float objects[9];
};

bool intersectAABB(struct Ray* ray, struct Intersection* intersection, struct AABB* aabb, float* dis)
{
	// Smits ray-box intersection test using slabs
	    // http://www.cs.utah.edu/~awilliam/box/box.pdf
	    float tmin, tmax, tymin, tymax, tzmin, tzmax;
	    float div = ray->div.x;
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

	    div = ray->div.y;
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

	    div = ray->div.z;
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
		        intersection->bary = (float3)(1.0f - s - t, s, t);
		        intersection->position = intersectionPoint;
		        intersection->distance = inter_div;
		        
		        intersection->objects[0] = a.x;	
		        intersection->objects[1] = a.y;
		        intersection->objects[2] = a.z;
		        intersection->objects[3] = b.x;
		        intersection->objects[4] = b.y;
		        intersection->objects[5] = b.z;
		        intersection->objects[6] = c.x;
		        intersection->objects[7] = c.y;
		        intersection->objects[8] = c.z;
		        
		    }	    
	
	
	return true;

}
bool pointInside(struct AABB* aabb, float3* point) {
	
		if (point->x > aabb->max.x || point->x < aabb->min.x)
			return false;
		if (point->y > aabb->max.y || point->y < aabb->min.y)
					return false;
		if (point->z > aabb->max.z || point->z < aabb->min.z)
					return false;
	return true;
}

float3 getNormal(struct Intersection* inter) {
	
	float3 v0 = (float3)(inter->objects[0],inter->objects[1],inter->objects[2]);
	float3 v1 = (float3)(inter->objects[3],inter->objects[4],inter->objects[5]);
	float3 v2 = (float3)(inter->objects[6],inter->objects[7],inter->objects[8]);

return normalize(cross(v1-v0,v2-v0));

}

float shading(struct Ray* ray, float3* normal) {
	float doted = 0.5f;//
	doted = dot(*normal, ray->direction);
	return fmin(1.0f, fabs(doted));
	
}
__kernel void compute_rays(__global struct Ray* rays,
							 float positionx,
							 float positiony,
							 float positionz,
							 float ax, 
							 float ay)
{		
	//printf("%i\n", sizeof(struct Intersection));

	for(int s=0; s<SAMPLES; s++) {
		struct Ray ray;
		float3 pos = (float3)(positionx, positiony, positionz);
		float offset = (s+0.5f)/SAMPLES; 
		float3 direction = (float3)(
				(get_global_id(0)+offset)/ax-WIDTH/(2.0f*ax),
				-((get_global_id(1)+offset)/ay-HEIGHT/(2.0f*ay)),
				-1.0f);
		ray.position = pos;
		ray.direction = normalize(direction);
		float3 div = (float3)(1.0f/ray.direction.x, 1.0f/ray.direction.y, 1.0f/ray.direction.z);
		ray.div = div;
		rays[(get_global_id(0)+get_global_id(1)*WIDTH)*SAMPLES+s] = ray;
	}
}

__kernel void compute_ao_rays(__global struct Ray* rays,
							 float positionx,
							 float positiony,
							 float positionz,
							 float ax, 
							 float ay)
{		
	

	
		struct Ray ray;
		float3 pos = (float3)(positionx, positiony, positionz);
		float offset =  0.5f; 
		float3 direction = (float3)(
				(get_global_id(0)+offset)/ax-WIDTH/(2.0f*ax),
				-((get_global_id(1)+offset)/ay-HEIGHT/(2.0f*ay)),
				-1.0f);
		ray.position = pos;
		ray.direction = normalize(direction);
		float3 div = (float3)(1.0f/ray.direction.x, 1.0f/ray.direction.y, 1.0f/ray.direction.z);
		ray.div = div;
		//printf("%f:%f:%f\n",direction.x, direction.y,direction.z);
		rays[(get_global_id(0)+get_global_id(1)*WIDTH)] = ray;
	
}

__kernel void compute_intersections(__global struct Intersection* intersections, 
										int ray_count,
										__global struct Ray* rays,
										__global struct cl_bvh_node* nodes) {
	
		struct Intersection intersection;
		intersection.distance = 1e6;
		struct Ray ray = rays[get_global_id(0)];
		
		int fringe[NODE_COUNT];
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
			//if(!intersectAABB(&ray, &intersection, &aabb, &distance))
				//continue;
			if(!pointInside(&aabb, &ray.position)) {
						
							float distance= 0.0f;
							if(!intersectAABB(&ray, &intersection, &aabb, &distance)) {
								continue;
							}
								
							if(distance > intersection.distance) 
								continue;
							
						}
				
			
			if (act_node.tri_counter != 0)
			    {
			        
			        
			        for (int i = 0; i < act_node.tri_counter; i++) {
			        	struct Triangle tri;
			        	tri.a = (float3)(act_node.objects[i*9+0],act_node.objects[i*9+1],act_node.objects[i*9+2]);
			        	tri.b = (float3)(act_node.objects[i*9+3],act_node.objects[i*9+4],act_node.objects[i*9+5]);
			        	tri.c = (float3)(act_node.objects[i*9+6],act_node.objects[i*9+7],act_node.objects[i*9+8]);
			        	
			        	
			        	intersectTriangle(&ray, &intersection, &tri);
			        }
			        continue;
			    }
			
			fringe[f] = act_node.leftIndex;
			
			f=f+1;
			
			fringe[f] = act_node.rightIndex;
			f=f+1;
		}
		
		intersections[get_global_id(0)] = intersection;
}

__kernel void compute_shading(__global struct Intersection* intersections, 
								__global struct Ray* rays,
								__global unsigned char* image) 
{
	
	float value = 0.0f;
		for(int n=0; n<SAMPLES; n++) {
			struct Ray ray = rays[(get_global_id(0)+get_global_id(1)*WIDTH)*SAMPLES+n];
			struct Intersection inter = intersections[(get_global_id(0)+get_global_id(1)*WIDTH)*SAMPLES+n];
			float3 normal = getNormal(&inter);
			
			value +=  shading(&ray, &normal);
			
		}
		
		value = value /(float)(SAMPLES);
	
		
		image[get_global_id(0)+get_global_id(1)*WIDTH] = (unsigned char)(value * 255.0f);
}

__kernel void compute_ao_point_rays(__global struct Ray* ao_rays,
					__global struct cl_bvh_node* nodes,
					__global struct Ray* rays,
					__global struct Intersection* intersections,
					__global float* aoRandoms) {
	//printf("fdfd");
	struct Ray ray = rays[get_global_id(0)+get_global_id(1)*WIDTH];
	struct Intersection inter = intersections[get_global_id(0)+get_global_id(1)*WIDTH];
	//printf("%f\n",inter.distance);
	float3 normal = getNormal(&inter);
	//printf("%f:%f:%f\n",normal.x,normal.y,normal.z);
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
		    
		    
		    
		

	
	float3 interPoint = inter.position + (float3)(1e-6 * normal.x,1e-6 * normal.y,1e-6 * normal.z) ;
		// Cast ao rays from intersection point into different directions
		
		for (int i = 0; i < AO_SAMPLES; ++i) {
			struct Intersection intersection;
			intersection.distance= 1e6;
			struct Ray ray;
			ray.position = interPoint;
			//printf("%f\n", intersection.distance);
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
		
			ao_rays[(get_global_id(0)+get_global_id(1)*WIDTH)*AO_SAMPLES+i] = ray;

		}

}

__kernel void compute_ao_factor(__global struct Intersection* intersections, 
									__global unsigned char* image,
									__global unsigned char* image_out) 
	{
		int Counter = 0;
		for(int n=0; n<AO_SAMPLES; n++) {
				
				struct Intersection inter = intersections[(get_global_id(0)+get_global_id(1)*WIDTH)*AO_SAMPLES+n];
				
				if(inter.distance <= AO_MAX_DISTANCE) {
					Counter++;
					
				}	
			}		
			float old = (float)(image[get_global_id(0)+get_global_id(1)*WIDTH]);
			//printf("%i\n",old);
			float value = (1.0f-((float)(Counter)/(float)(AO_SAMPLES)))*(old/255.0f);
			//printf("%i\n",Counter);
			image_out[get_global_id(0)+get_global_id(1)*WIDTH] = (unsigned char)(value * 255.0f);
		
}



