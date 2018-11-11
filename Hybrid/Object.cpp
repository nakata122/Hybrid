#include <GL/glew.h>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <ios>
#include "BVH.h"
#define BUFFER_OFFSET(i) ((void*)(i))


class Object {
	std::vector<float> normals;
	std::vector<glm::vec4> vertices;
	std::vector<BVH::Triangle> triangles;
	std::vector<GLuint> elements;
	glm::mat4 Model = glm::mat4(1.0f);
	glm::mat4 Scale = glm::mat4(1.0f);
	glm::mat4 Translate = glm::mat4(1.0f);
	glm::mat4 Rotate = glm::mat4(1.0f);
	glm::mat4 mvp = glm::mat4(1.0f);

	GLuint vbo, ibo, vao, ssbo;
public:
	void LoadObj(const char* filename)
	{
		std::ifstream in(filename, std::ios::in);
		if (!in)
		{
			std::cout << "Cannot open " << filename << std::endl; exit(1);
		}

		std::string line;
		while (std::getline(in, line))
		{
			if (line.substr(0, 2) == "v ")
			{
				std::istringstream s(line.substr(2));
				glm::vec4 v; s >> v.x; s >> v.y; s >> v.z; v.w = 1.0f;
				vertices.push_back(v);
			}
			else if (line.substr(0, 2) == "f ")
			{
				std::istringstream s(line.substr(2));
				char dummy;
				GLuint a, b, c, d, e, f, t;
				//s >> a >> dummy >> d >> dummy >> t;
				//s >> b >> dummy >> e >> dummy >> t;
				//s >> c >> dummy >> f >> dummy >> t;
				s >> a;
				s >> b;
				s >> c;
				a--; b--; c--;
				triangles.push_back(BVH::Triangle(vertices[a], vertices[b], vertices[c]));
				elements.push_back(a); elements.push_back(b); elements.push_back(c);
				//normals.push_back(d); normals.push_back(e); normals.push_back(f);

				//std::cout << a << " " << b << " " << c << std::endl;
			}
			else if (line[0] == '#')
			{
				/* ignoring this line */
			}
			else
			{
				/* ignoring this line */
			}
		}
		
		/*normals.resize(vertices.size(), glm::vec3(0.0, 0.0, 0.0));
		for (int i = 0; i < elements.size(); i += 3)
		{
			GLuint ia = elements[i];
			GLuint ib = elements[i + 1];
			GLuint ic = elements[i + 2];
			glm::vec3 normal = glm::normalize(glm::cross(
				glm::vec3(vertices[ib]) - glm::vec3(vertices[ia]),
				glm::vec3(vertices[ic]) - glm::vec3(vertices[ia])));
			normals[ia] = normals[ib] = normals[ic] = normal;
		}*/
	}
	
	void CreateVBO() {

		std::cout << vertices.size() << std::endl;

		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), &vertices[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);


		glGenBuffers(1, &ibo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, elements.size() * sizeof(GLuint), &elements[0], GL_STATIC_DRAW);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);

		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, (verts.size()/3) * sizeof(float), 0);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glEnableVertexAttribArray(0);
		glBindVertexArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		
	}
	static bool sortX(const glm::vec4 &a, const glm::vec4 &b) {
		if (a.x < b.x) return 1;
		else if (a.x > b.x) return 0;
	}
	static bool sortXTri(const BVH::Triangle &a, const BVH::Triangle &b) {
		if ((a.v0.x + a.v1.x + a.v2.x)/3 < (b.v0.x + b.v1.x + b.v2.x) / 3) return 1;
		else return 0;
	}
	void createSSBO() {

		BVH bvh_builder;

		bvh_builder.CreateBVH(triangles);
		std::sort(triangles.begin(), triangles.end(), sortXTri);
		glGenBuffers(1, &ssbo);
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
		glBufferData(GL_SHADER_STORAGE_BUFFER, triangles.size() * sizeof(BVH::Triangle), &triangles[0], GL_STATIC_DRAW);
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo); // Buffer Binding 1
		glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
	}
	void ScaleObject(glm::mat4 Scaling){
		Scale = Scaling;
	}
	void TranslateObject(glm::mat4 Translation) {
		Translate = Translation;
	}
	void RotateObject(glm::mat4 Rotation) {
		Rotate = Rotation;
	}
	void ProjectObject(glm::mat4 Camera, glm::mat4 Projection, GLuint id) {
		Model = Translate * Rotate * Scale;
		mvp = Projection * Camera * Model;
		// Send our transformation to the currently bound shader, in the "MVP" uniform
		// This is done in the main loop since each model will have a different MVP matrix (At least for the M part)
		glUniformMatrix4fv(id, 1, GL_FALSE, &mvp[0][0]);
	}
	void DrawObject() {
		glBindVertexArray(vao);
		//glDrawArrays(GL_POINTS, 0, elements.size());
		//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
		glDrawElements(GL_POINTS, elements.size(), GL_UNSIGNED_INT, NULL);
		glBindVertexArray(0);

	}
	void SendData() {
		glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, ssbo); // Buffer Binding 1
	}
};