#pragma once
#include <GL/glew.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <glm/glm.hpp>

class BVH
{
public:
	int br = 0;
	struct Triangle {
		glm::vec4 v[3];
		Triangle(glm::vec4 a, glm::vec4 b, glm::vec4 c) { v[0] = a; v[1] = b; v[2] = c; };
		Triangle(glm::vec4 temp[3]) { v[0] = temp[0]; v[1] = temp[1]; v[2] = temp[2]; };
	};
	struct BVHNode
	{
		glm::vec4 bboxMin;
		glm::vec4 bboxMax;
		int secondChild;
		int id;
		int isLeaf;
		float padding;
	};
	std::vector<BVHNode> boundingVolume;
	BVH();
	void CreateBVH(const std::vector<BVH::Triangle>::iterator begin, const std::vector<BVH::Triangle>::iterator end);

};

