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
class LightComponent;
class Shader;

class Renderer
{
public:
	void initGl();
	void resizeGl(int __w, int __h);

	void drawDrawable(DrawableComponent* __drawableComponent);
	void initDrawable(DrawableComponent* __drawableComponent);
	void freeDrawable(DrawableComponent* __drawableComponent);

	void setActiveCamera(CameraComponent* __cameraComponent);
	void addLight(LightComponent* __lightComponent);
	void startFrame();
	void endFrame(GLuint qtFramebuffer);


private:
	void _screenbufferInit(int __w, int __h);

	CameraComponent* _activeCamera = nullptr;
	GLsizei _sceneIndices = 0;
	std::vector<LightComponent*> _lights;

	void _useShader(DrawableComponent* __drawableComponent);
	Shader* _screenquadShader = nullptr;
	GLuint _fbo;
	GLuint _rbo;
	GLuint _screenquadVAO;
	GLuint _screenquadVBO;
	GLuint _textureColorbuffer;
	bool _init = false;
};

