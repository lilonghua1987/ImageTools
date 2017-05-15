#pragma once
#include "opengl.h"
#include<iostream>


class LFrame
{
public:
	LFrame();
	LFrame(unsigned int width, unsigned int height, std::string title);
	~LFrame(void);
	
private:

	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void errorCallback(int error, const char* description);
	bool init(void);

public:
	GLFWwindow* window;
	unsigned int width;
	unsigned int height;
	std::string title;

	void draw(void);
private:
	static unsigned int frameCount;
public:
	void destory(void);
};

