#pragma once
#include <GL/glew.h>
#include <vector>
#include <glm/glm.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <ios>
class Utils
{
public:
	static GLuint CreateQuad() {
		GLuint vao = 0, vbo = 0;
		float verts[] = { 
			-1.0f, -1.0f, 0.0f, 0.0f,
			-1.0f, 1.0f, 0.0f, 1.0f, 
			1.0f, -1.0f, 1.0f, 0.0f, 
			1.0f, 1.0f, 1.0f, 1.0f 
		};
		glGenBuffers(1, &vbo);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBufferData(GL_ARRAY_BUFFER, 16 * sizeof(float), verts, GL_STATIC_DRAW);
		glGenVertexArrays(1, &vao);
		glBindVertexArray(vao);
		glEnableVertexAttribArray(0);
		GLintptr stride = 4 * sizeof(float);
		glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, NULL);
		glEnableVertexAttribArray(1);
		GLintptr offset = 2 * sizeof(float);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (GLvoid *)offset);
		return vao;
	}
	static GLuint CreateQuadProgram() {

		std::string vertexShader =
			"#version 430 \n"
			"layout (location = 0) in vec2 vp; \n"
			"layout (location = 1) in vec2 vt; \n"
			"out vec2 st; \n"
			"void main () { \n"
			"	st = vt; \n"
			"	gl_Position = vec4 (vp, 0.0, 1.0); \n"
			"} \n";

		std::string fragmentShader =
			"#version 430 \n"
			"in vec2 st; \n "
			"uniform sampler2D img; \n"
			"out vec4 fc; \n"
			"void main () { \n"
			"	fc = texture (img, st); \n"
			"} \n";

		GLuint program =  CreateShader(vertexShader, fragmentShader);
		return program;
	}
	static GLuint GenerateTexture() {
		GLuint texture;
		glGenTextures(1, &texture);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texture);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, 640, 480, 0, GL_RGBA, GL_FLOAT,
			NULL);
		glBindImageTexture(0, texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
		return texture;
	}
	static int CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
		unsigned int program = glCreateProgram();
		unsigned int vs = CompileShader(vertexShader, GL_VERTEX_SHADER);
		unsigned int fs = CompileShader(fragmentShader, GL_FRAGMENT_SHADER);

		glAttachShader(program, vs);
		glAttachShader(program, fs);
		glLinkProgram(program);
		glValidateProgram(program);

		glDeleteShader(vs);
		glDeleteShader(fs);

		return program;
	}
	static int CreateComputeShader(const std::string& computeShader) {
		unsigned int program = glCreateProgram();
		unsigned int shader = CompileShader(computeShader, GL_COMPUTE_SHADER);

		glAttachShader(program, shader);
		glLinkProgram(program);
		glValidateProgram(program);

		glDeleteShader(shader);

		return program;
	}
	static unsigned int CompileShader(const std::string& source, unsigned int type) {
		unsigned int id = glCreateShader(type);
		const char* src = source.c_str();
		glShaderSource(id, 1, &src, nullptr);
		glCompileShader(id);

		int err;
		glGetShaderiv(id, GL_COMPILE_STATUS, &err);
		if (err == GL_FALSE) {
			int maxLength = 1000;
			std::vector<GLchar> errorLog(maxLength);
			glGetShaderInfoLog(id, maxLength, &maxLength, &errorLog[0]);
			for (auto c : errorLog) std::cout << c;
			std::cout << std::endl;
		}

		return id;
	}
};

