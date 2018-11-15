#include "BVH.h"

BVH::BVH()
{
}

static bool sortX(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v[0].x + a.v[1].x + a.v[2].x) / 3 < (b.v[0].x + b.v[1].x + b.v[2].x) / 3) return 1;
	else return 0;
}
static bool sortY(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v[0].y + a.v[1].y + a.v[2].y) / 3 < (b.v[0].y + b.v[1].y + b.v[2].y) / 3) return 1;
	else return 0;
}
static bool sortZ(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v[0].z + a.v[1].z + a.v[2].z) / 3 < (b.v[0].z + b.v[1].z + b.v[2].z) / 3) return 1;
	else return 0;
}
BVH::BVHNode enclose(BVH::BVHNode node, BVH::Triangle data) {
	for (int j = 0; j < 3; j++) {
		if (data.v[j].x < node.bboxMin.x) node.bboxMin.x = data.v[j].x;
		if (data.v[j].y < node.bboxMin.y) node.bboxMin.y = data.v[j].y;
		if (data.v[j].z < node.bboxMin.z) node.bboxMin.z = data.v[j].z;

		if (data.v[j].x > node.bboxMax.x) node.bboxMax.x = data.v[j].x;
		if (data.v[j].y > node.bboxMax.y) node.bboxMax.y = data.v[j].y;
		if (data.v[j].z > node.bboxMax.z) node.bboxMax.z = data.v[j].z;
	}
	return node;
}

std::vector<BVH::Triangle>::iterator first;
void BVH::CreateBVH(const std::vector<BVH::Triangle>::iterator begin, const std::vector<BVH::Triangle>::iterator end) {
	if (br == 0) first = begin;
	br++;

	const int objectCount = end - begin;

	BVHNode node;
	node.bboxMax = glm::vec4(-100000, -1000000, -1000000, 1);
	node.bboxMin = glm::vec4(100000, 1000000, 1000000, 1);

	for (auto iter = begin; iter != end; ++iter) {
		node = enclose(node, iter->v);
	}
	
	int index = boundingVolume.size();
	boundingVolume.push_back(node);

	float diffX = node.bboxMax.x - node.bboxMin.x;
	float diffY = node.bboxMax.y - node.bboxMin.y;
	float diffZ = node.bboxMax.z - node.bboxMin.z;

	//std::cout << index << " " << node.bboxMin.x << " " << node.bboxMax.x << std::endl;
	if (objectCount <= 2) {
		boundingVolume[index].isLeaf = 1;
		boundingVolume[index].id = begin - first;
	} else {
		boundingVolume[index].isLeaf = 0;
		const std::vector<Triangle>::iterator medianIter = begin + (end - begin) / 2;

		if(diffX > diffY && diffX > diffZ)
			std::nth_element(begin, medianIter, end, sortX);
		else if(diffY > diffX && diffY > diffZ)
			std::nth_element(begin, medianIter, end, sortY);
		else
			std::nth_element(begin, medianIter, end, sortZ);

		boundingVolume[index].id = boundingVolume.size();
		const int firstChildIndex = boundingVolume.size();	
		CreateBVH(begin, medianIter);
		boundingVolume[index].secondChild = boundingVolume.size();
		CreateBVH(medianIter, end);
	}
}

void RecursiveSplit(int l, int r, std::vector<BVH::Triangle> data) {
	
}

