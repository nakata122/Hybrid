#include "BVH.h"

BVH::BVH()
{
}


void BVH::CreateBVH(std::vector<Triangle> data) {
	
}

void RecursiveSplit(int l, int r) {
	
}

static bool sortX(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v0.x + a.v1.x + a.v2.x) / 3 < (b.v0.x + b.v1.x + b.v2.x) / 3) return 1;
	else return 0;
}
static bool sortY(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v0.y + a.v1.y + a.v2.y) / 3 < (b.v0.y + b.v1.y + b.v2.y) / 3) return 1;
	else return 0;
}
static bool sortZ(const BVH::Triangle &a, const BVH::Triangle &b) {
	if ((a.v0.z + a.v1.z + a.v2.z) / 3 < (b.v0.z + b.v1.z + b.v2.z) / 3) return 1;
	else return 0;
}