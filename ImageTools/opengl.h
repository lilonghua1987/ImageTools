#ifndef OPENGL_GLFW_H
#define OPENGL_GLFW_H

#pragma once

//#include <GL/glew.h>
//
//#if defined(_WIN32)
//#include <GL/wglew.h>
//#elif !defined(__APPLE__) || defined(GLEW_APPLE_GLX)
//#include <GL/glxew.h>
//#endif

#define GLFW_INCLUDE_GLU

#include "GLFW/glfw3.h"


#ifdef _DEBUG
#  ifdef _WIN32
#     pragma comment(lib, "lib/x86/glew32d.lib")
#     pragma comment(lib, "lib/x86/glfw3dll.lib")
#  else
#     pragam comment(lib, "lib/x64/glew32d.lib")
#     pragam comment(lib, "lib/x64/glfw3dll.lib")
#  endif
#else
#  ifdef _WIN32
#     pragma comment(lib, "lib/x86/glew32.lib")
#     pragma comment(lib, "lib/x86/glfw3.lib")
#  else
#     pragam comment(lib, "lib/x64/glew32.lib")
#     pragam comment(lib, "lib/x64/glfw3.lib")
#  endif
#endif

#pragma comment( lib, "OpenGL32.lib")

#ifdef GLFW_INCLUDE_GLU
#   pragma comment( lib, "GlU32.lib")
#endif

//#define FREEGLUT_LIB_PRAGMAS 0
//#include "freeglut/freeglut.h"
//#pragma comment(lib, "opengl32.lib")
//#pragma comment(lib, "lib/freeglut.lib")


#include "Image.h"

class OpenGLTools
{
public:
	OpenGLTools();
	~OpenGLTools();

	//OpenGL
	static Image Image2GL(const Image& img);
	//SFML
	static Image ImageToSFML(const Image& src);
};

#endif