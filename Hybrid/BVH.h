#pragma once
#include <GL/glew.h>
#include <vector>
#include <glm/glm.hpp>

class BVH
{
public:
	struct Triangle {
		glm::vec4 v0, v1, v2;
		Triangle(glm::vec4 a, glm::vec4 b, glm::vec4 c) { v0 = a; v1 = b; v2 = c; };
	};
	struct BVHNode
	{
		glm::vec4 bboxMin;
		glm::vec4 bboxMax;
	};
	std::vector<BVHNode> boundingVolume;
	BVH();
	void CreateBVH(std::vector<Triangle> data);
};

