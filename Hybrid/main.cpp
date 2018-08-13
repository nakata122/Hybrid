#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include "Object.cpp"
#include "Utils.h"

const int SCREEN_WIDTH = 640;
const int SCREEN_HEIGHT = 480;
double posX=0, posY=0, oldX=0, oldY=0;


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

	std::string computeShader =
		"#version 430 \n"
		"layout (local_size_x = 1, local_size_y = 1) in; \n"
		"layout (rgba32f, binding = 0) uniform image2D img_output; \n"
		"const float PI = 3.1415926535897932384626433832795;"
		"struct Sphere \n"
		"{ \n"
		"	vec3 center; \n"
		"	float radius; \n"
		"	vec3 surfaceColor; \n"
		"	float reflection; \n"
		"}; \n"
		"struct Ray \n"
		"{ \n"
		"	vec3 origin; \n"
		"	vec3 dir; \n"
		"	float power; \n"
		"}; \n"
		"void main () { \n"
		"	vec4 pixel = vec4 (0.0, 0.0, 0.0, 1.0); \n"
		"	ivec2 pixel_coords = ivec2 (gl_GlobalInvocationID.xy); \n"
		"	float max_x = 5.0; \n"
		"	float max_y = 5.0; \n"
		"	ivec2 dims = imageSize (img_output); \n"
		"	float fov = 30.0; \n"
		"	float aspectratio = dims.x / float(dims.y); \n"
		"	float angle = tan(PI * 0.5 * fov / 180.);"
		"	float x = ((2.0 * (pixel_coords.x + 0.5) / dims.x) - 1.0) * aspectratio; \n"
		"	float y = 1.0 - (2.0 * (pixel_coords.y + 0.5f) / dims.y); \n"
		"	Ray camera = Ray(vec3 (x * max_x, y * max_y, 0.0), vec3 (0.0, 0.0, -1.0), 1.0);"
		"	Sphere sphere1 = Sphere(vec3(0.0, 0.0, -10.0), 2.0, vec3(0.4, 0.4, 1.0), 0.5); \n"
		"	vec3 omc = sphere1.center - camera.origin; \n"
		"	float b = dot (camera.dir, omc); \n"
		"	float c = dot (omc, omc) - sphere1.radius * sphere1.radius; \n"
		"	float bsqmc = b * b - c; \n"

		"	float t = b + sqrt(bsqmc); \n"
		"	vec3 hitpoint = camera.origin + camera.dir * t; \n"
		"	vec3 normal = normalize(hitpoint - sphere1.center); \n"
		"	float cosine_factor = dot(normal, camera.dir); \n"

		"	if (bsqmc >= 0.0) { \n"
		"		pixel = vec4 (vec3(0.4,0.4,1.0) * cosine_factor, 1.0); \n"
		"	} \n"
		"	imageStore (img_output, pixel_coords, pixel); \n"
		"} \n";

	glm::mat4 ProjectionMatrix = glm::perspective(glm::radians(45.0f), (float)SCREEN_WIDTH / (float)SCREEN_HEIGHT, 0.1f, 1000.0f);

	glm::mat4 CameraMatrix = glm::lookAt(
		glm::vec3(4, 3, 3), // the position of your camera, in world space
		glm::vec3(0,0,0),   // where you want to look at, in world space
		glm::vec3(0,1,0)    // probably glm::vec3(0,1,0), but (0,-1,0) would make you looking upside-down, which can be great too
	);

	Object box;
	box.LoadObj("dragon.obj");
	box.CreateVBO();
	
	
	unsigned int shader = Utils::CreateShader(vertexShader, fragmentShader);
	glUseProgram(shader);

	GLuint quadVao = Utils::CreateQuad();
	GLuint quadProgram = Utils::CreateQuadProgram();
	GLuint rayProgram = Utils::CreateComputeShader(computeShader);
	GLuint tex = Utils::GenerateTexture();


	double phi = 0;
	double theta = 0;

	glDisable(GL_CULL_FACE);

	/* Loop until the user closes the window */
	while (!glfwWindowShouldClose(window))
	{
		/* Render here */
		glClear(GL_COLOR_BUFFER_BIT);


		//Camera movement
		CameraMatrix = glm::lookAt(
			glm::vec3(20*cos(phi)*sin(theta),20*sin(phi)*sin(theta), 20*cos(theta)), // the position of your camera, in world space
			glm::vec3(0, 0, 0),   // where you want to look at, in world space
			glm::vec3(0, cos(phi), 0)        // probably glm::vec3(0,1,0), but (0,-1,0) would make you looking upside-down, which can be great too
		);

		/*glUseProgram(rayProgram);
		glDispatchCompute((GLuint)640, (GLuint)480, 1);

		glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

		glUseProgram(quadProgram);
		glBindVertexArray(quadVao);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, tex);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);*/

		glUseProgram(shader);
		GLuint mvpId = glGetUniformLocation(shader, "MVP");
		box.ProjectObject(CameraMatrix, ProjectionMatrix, mvpId);
		box.DrawObject();

		/* Swap front and back buffers */
		glfwSwapBuffers(window);

		/* Poll for and process events */
		glfwPollEvents();

		theta += (posX - oldX)*0.01f;
		phi += (posY - oldY)*0.01f;

		oldX = posX;
		oldY = posY;
	}

	glfwTerminate();
	return 0;
}