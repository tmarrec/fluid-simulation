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

	void drawDrawable(DrawableComponent* __drawableComponent);
	void initDrawable(DrawableComponent* __drawableComponent);
	void freeDrawable(DrawableComponent* __drawableComponent);
	void setActiveCamera(CameraComponent* __cameraComponent);
	void clear() const;


private:
	CameraComponent* _activeCamera = nullptr;
	GLsizei _sceneIndices = 0;

	void _useShader(DrawableComponent* __drawableComponent);
};

