#include "LFrame.h"

unsigned int LFrame::frameCount = 0;

LFrame::LFrame(void)
	: window(NULL)
	, width(1024)
	, height(768)
	, title("JFrame")
{
	if(init())
	{
		window = glfwCreateWindow(width, height, title.c_str(), NULL, NULL);
		if (!window)
		{
			glfwTerminate();
			//exit(EXIT_FAILURE);
		}

		glfwSetErrorCallback(errorCallback);

		//binding context
		glfwMakeContextCurrent(window);
		glfwSetKeyCallback(window, keyCallback);
		frameCount++;
	}
}


LFrame::LFrame(unsigned int width, unsigned int height, std::string title)
	: window(NULL)
	, width(width)
	, height(height)
	, title(title)
{
	if(init())
	{
		window = glfwCreateWindow(width, height, title.c_str(), NULL, NULL);
		if (!window)
		{
			glfwTerminate();
			//exit(EXIT_FAILURE);
		}

		glfwSetErrorCallback(errorCallback);

		//binding context
		glfwMakeContextCurrent(window);
		glfwSetKeyCallback(window, keyCallback);
		frameCount++;
	}
}


LFrame::~LFrame(void)
{
	destory();
	if (frameCount < 1)
	{
		glfwTerminate(); 
	}
}


void LFrame::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
		glfwSetWindowShouldClose(window, GL_TRUE);
}


void LFrame::errorCallback(int error, const char* description)
{
	fputs(description, stderr);
}


bool LFrame::init(void)
{
	if (frameCount < 1)
	{
		if (!glfwInit())
		{
			//exit(EXIT_FAILURE);
			return false;
		}
	}

	return true;
}


void LFrame::draw(void)
{
	// Enable Z-buffer read and write
	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glClearDepth(1.f);

	GLuint texture;
	Image img = OpenGLTools::Image2GL("img/imR.png");

	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D, texture);
	gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, img.width, img.height, GL_RGB, GL_UNSIGNED_BYTE, img.getPtr<BYTE>());
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);

	// Bind the texture
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture);

	while (!glfwWindowShouldClose(window))
	{
		float ratio;
		int width, height;
		glfwGetFramebufferSize(window, &width, &height);
		ratio = width / (float) height;
		glViewport(0, 0, width, height);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		//正交投影 物体大小不变
		//glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
		//透视投影
		glFrustum(-ratio*2, ratio, -1.f, 1.f, 1.f,-10.5f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glPushMatrix();
		glRotatef((float) glfwGetTime() * 50.f, 0.f, 0.f, 1.f);
		glRotatef((float) glfwGetTime() * 50.f, 0.f, 1.f, 0.f);

		glBegin(GL_TRIANGLES);
			//glColor3f(1.f, 0.f, 0.f);
			glTexCoord2f(1.f,0.f);
			glVertex3f(-0.6f, -0.4f, 0.5f);
			//glColor3f(0.f, 1.f, 0.f);
			glTexCoord2f(0.f,1.f);
			glVertex3f(0.6f, -0.4f, 0.7f);
			//glColor3f(0.f, 0.f, 1.f);
			glTexCoord2f(1.f,1.f);
			glVertex3f(0.f, 0.6f, -0.7f);
		glEnd();
		glPopMatrix();

		glfwSwapBuffers(window);
		glfwPollEvents();
	}
	glDeleteTextures(1, &texture);
	destory();
}


void LFrame::destory(void)
{
	if (window)
	{
		glfwDestroyWindow(window);
		window = NULL;
	    frameCount--;
	}
}
