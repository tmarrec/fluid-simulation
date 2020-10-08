#pragma once

#include <GL/gl.h>
#include <QOpenGLContext>
#include <QOpenGLWidget>
#include <cstdint>
#include <memory>

#include "Shader.h"
#include "src/ECS.h"


class DrawableComponent;
class CameraComponent;

class Renderer
{
public:
	void initGl();
	void resizeGl(int w, int h) const;

	void initDrawable(DrawableComponent* __drawableComponent);
	void freeDrawable();
	void draw();
	void useShader(DrawableComponent* __drawableComponent);
	void setActiveCamera(CameraComponent* __cameraComponent);

private:
	CameraComponent* _activeCamera = nullptr;
	std::vector<GLuint> _VAOs;
	GLuint _VAO;
	std::uint_fast64_t arrayObjetsNB = 0;
	GLsizei _sceneIndices = 0;
};

