#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include "Object.cpp"
#include "Utils.h"

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
double posX=0, posY=0, oldX=0, oldY=0;


std::string loadshader(const char* filename)
{
	std::ifstream file;
	file.open(filename, std::ios::in); // opens as ASCII!
	if (!file) return "";
	std::string str((std::istreambuf_iterator<char>(file)),
		std::istreambuf_iterator<char>());

	file.close();

	return str;
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	glfwGetCursorPos(window, &posX, &posY);
	//std::cout << posX << " " << oldX << " " << posY << " " << oldY << std::endl;
}

int main(void)
{
	GLFWwindow* window;

	/* Initialize the library */
	if (!glfwInit())
		return -1;

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
	if (!window)
	{
		glfwTerminate();
		return -1;
	}

	/* Make the window's context current */
	glfwMakeContextCurrent(window);
	glfwSetCursorPosCallback(window, cursor_position_callback);

	glewInit();

	std::cout << glGetString(GL_VERSION) << std::endl;
	

	std::string vertexShader = 
		"#version 330 core \n"
		"layout(location = 0) in vec3 position; \n"
		"uniform mat4 MVP; \n"
		"void main(void){ \n"
		"	gl_Position = MVP * vec4(position, 1.0); \n"
		"} \n";

	std::string fragmentShader =
		"#version 330 core \n"
		"layout(location = 0) out vec4 color; \n"
		"void main(){ \n"
		"	color = vec4(1.0, 0.0, 0.0, 1.0); \n"
		"} \n";

	std::string computeShader = loadshader("shader.comp");

	glm::mat4 ProjectionMatrix = glm::perspective(glm::radians(45.0f), (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT, 0.1f, 1000.0f);

	glm::mat4 CameraMatrix = glm::lookAt(
		glm::vec3(4, 3, 3), // the position of your camera, in world space
		glm::vec3(0,0,0),   // where you want to look at, in world space
		glm::vec3(0,1,0)    // probably glm::vec3(0,1,0), but (0,-1,0) would make you looking upside-down, which can be great too
	);

	Object box;
	box.LoadObj("dragon.obj");
	box.createSSBO();
	
	
	unsigned int shader = Utils::CreateShader(vertexShader, fragmentShader);
	glUseProgram(shader);

	GLuint quadVao = Utils::CreateQuad();
	GLuint quadProgram = Utils::CreateQuadProgram();
	GLuint rayProgram = Utils::CreateComputeShader(computeShader);
	GLuint tex = Utils::GenerateTexture();

	GLint loc = glGetUniformLocation(rayProgram, "pos");
	GLint direction = glGetUniformLocation(rayProgram, "direction");
	

	double phi = 0;
	double theta = 0;

	glDisable(GL_CULL_FACE);
	glfwSwapInterval(0);

	double lastTime = glfwGetTime();
	int nbFrames = 0;

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		/* Render here */
		glClear(GL_COLOR_BUFFER_BIT);


		//Camera movement
		CameraMatrix = glm::lookAt(
			glm::vec3(20*cos(phi)*sin(theta),20*sin(phi)*sin(theta), 20*cos(theta)), // the position of your camera, in world space
			glm::vec3(0, 0, 0),   // where you want to look at, in world space
			glm::vec3(0, cos(phi), 0)        // probably glm::vec3(0,1,0)
		);

		glUseProgram(rayProgram);
		if (loc != -1)
		{
			glUniform2f(loc, (SCREEN_WIDTH/2 - posX)/100, (SCREEN_HEIGHT/2 - posY)/100);
		}
		if (direction != -1)
		{
			glUniform3f(direction, (rand() % 100) / 100, (rand() % 100) / 100, (rand() % 100) / 100);
		}
		//box.createSSBO();
		glDispatchCompute(SCREEN_WIDTH/8, SCREEN_HEIGHT/8, 1);

		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		glUseProgram(quadProgram);
		glBindVertexArray(quadVao);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, tex);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

		/*glUseProgram(shader);
		GLuint mvpId = glGetUniformLocation(shader, "MVP");
		box.ProjectObject(CameraMatrix, ProjectionMatrix, mvpId);
		box.DrawObject();*/

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		theta += (posX - oldX)*0.01f;
		phi += (posY - oldY)*0.01f;

		oldX = posX;
		oldY = posY;

		// Measure speed
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0) { // If last prinf() was more than 1 sec ago
			// printf and reset timer
			printf("%f ms/frame\n", 1000.0 / double(nbFrames));
			nbFrames = 0;
			lastTime += 1.0;
		}
	}

	glfwTerminate();
	return 0;
}