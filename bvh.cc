#include <iostream>
#include <limits>
#include <stdexcept>
#include <list>

#include "visual.h"
#include "bvh.h"

void BVH::buildBVH(Mesh const& mesh)
{
    delete this->root;
    this->root = new BVHNode();
    std::vector<unsigned int> faceIDs(mesh.faces.size() / 3);
    for (std::size_t i = 0; i < faceIDs.size(); ++i)
        faceIDs[i] = i;
    root->insert(mesh, &faceIDs);
}

bool BVH::intersect(Ray const& ray, Intersection* intersection, float const maxDistance) const
{
    if (this->root == NULL)
        return false;

    intersection->distance = std::numeric_limits<float>::max();
    if (this->root->intersect(ray, intersection, maxDistance))
        return true;

    intersection->distance = 0.0f;
    return false;
}

void BVH::addAABBsToMesh(Mesh* mesh, int nrOfDots) const
{
    if (this->root == NULL)
        return;

    typedef std::list<BVHNode const*> NodeList;
    NodeList nodes;
    nodes.push_back(this->root);
    while (!nodes.empty())
    {
        BVHNode const* node = nodes.back();
        nodes.pop_back();
        sampleAABBToMesh(mesh, nrOfDots, node->aabb);
        if (node->left)
            nodes.push_back(node->left);
        if (node->right)
            nodes.push_back(node->right);
    }
}

int BVH::printStatistics(std::ostream& out)
{
    if (this->root == NULL)
        return 0;

    int numInnerNodes = 0;
    int minLeafDepth = std::numeric_limits<int>::max();
    int maxLeafDepth = 0;
    int numLeafs = 0;

    typedef std::list<std::pair<BVHNode const*, int> > NodeLevelList;
    NodeLevelList nodes;
    nodes.push_back(std::make_pair(this->root, 1));
    while (!nodes.empty())
    {
        BVHNode const* node = nodes.back().first;
        int const level = nodes.back().second;
        nodes.pop_back();
        if (node->left == NULL && node->right == NULL)
        {
            minLeafDepth = std::min(minLeafDepth, level);
            maxLeafDepth = std::max(maxLeafDepth, level);
            numLeafs += 1;
        }
        else
        {
            numInnerNodes += 1;
            if (node->left)
                nodes.push_back(std::make_pair(node->left, level + 1));
            if (node->right)
                nodes.push_back(std::make_pair(node->right, level + 1));
        }
    }

    out << "  Number of inner nodes: " << numInnerNodes << std::endl;
    out << "  Number of leaf nodes: " << numLeafs << std::endl;
    out << "  Minimum leaf depth: " << minLeafDepth << std::endl;
    out << "  Maximum leaf depth: " << maxLeafDepth << std::endl;
    return numInnerNodes + numLeafs;
}
void BVH::getBVH(cl_bvh_node* scene_array) const {
	int counter =0;
	parseBVH(this->root,scene_array, counter);
}

int BVH::parseBVH(BVHNode const* node, cl_bvh_node* scene_array, int &counter) const {
	// Create AABB as cl struct
	//AABB note_aabb;
	//note_aabb.max= node.aabb.max;
	//note_aabb.min= node.aabb.min;
	// Create bvh node as cl struct
	//std::cout << counter << std::endl;
	cl_bvh_node cl_node;

	int index = counter;
	counter++;
	//cl_node.leftIndex= counter+1;
	//cl_node.rightIndex= counter+2;
	cl_node.min[0] = node->aabb.getAABBMin()[0];
	cl_node.min[1] = node->aabb.getAABBMin()[1];
	cl_node.min[2] = node->aabb.getAABBMin()[2];
	cl_node.max[0] = node->aabb.getAABBMax()[0];
	cl_node.max[1] = node->aabb.getAABBMax()[1];
	cl_node.max[2] = node->aabb.getAABBMax()[2];
	cl_node.tri_counter = node->triangles.size();
	//cl_node.cl_aabb= cl_aabb;
	// Add cl bvh node to array and descent

	if(node->triangles.empty()) {
		cl_node.leftIndex = parseBVH(node->left, scene_array, counter);
		cl_node.rightIndex = parseBVH(node->right, scene_array, counter);
	} else {
		for(int i= 0; i<cl_node.tri_counter;i++) {
//std::cout <<"test"<<std::endl;
			cl_node.objects[i*9+0] = node->triangles.at(i)[0][0];
			cl_node.objects[i*9+1] = node->triangles.at(i)[0][1];
			cl_node.objects[i*9+2] = node->triangles.at(i)[0][2];
			cl_node.objects[i*9+3] = node->triangles.at(i)[1][0];
			cl_node.objects[i*9+4] = node->triangles.at(i)[1][1];
			cl_node.objects[i*9+5] = node->triangles.at(i)[1][2];
			cl_node.objects[i*9+6] = node->triangles.at(i)[2][0];
			cl_node.objects[i*9+7] = node->triangles.at(i)[2][1];
			cl_node.objects[i*9+8] = node->triangles.at(i)[2][2];

		}
	}
	scene_array[index]= cl_node;
	return index;

}
