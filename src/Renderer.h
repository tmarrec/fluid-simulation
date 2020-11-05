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
class ECS_Manager;

class Renderer
{
public:
	Renderer(std::shared_ptr<ECS_Manager> __ECS_manager);
	void initGl();
	void resizeGl(int __w, int __h);
	void draw(double __deltaTime, int __qtFramebuffer);

	void drawDrawable(DrawableComponent* __drawableComponent);
	void initDrawable(DrawableComponent* __drawableComponent);
	void freeDrawable(DrawableComponent* __drawableComponent);

	void setActiveCamera(CameraComponent* __cameraComponent);
	void addLight(LightComponent* __lightComponent);

	void switchPolygonmode();


private:
	void _depthMapPass();
	void _colorPass(int __qtFramebuffer);
	void _screenbufferInit(int __w, int __h);

	double _deltaTime = 0;

	std::shared_ptr<ECS_Manager> _ECS_manager = nullptr;
	CameraComponent* _activeCamera = nullptr;
	GLsizei _sceneIndices = 0;
	std::vector<LightComponent*> _lights;

	void _useShader(DrawableComponent* __drawableComponent);
	bool _init = false;
	GLenum _polygonMode = GL_FILL;
	int _glWidth = 0;
	int _glHeight = 0;

	// Depth Map framebuffer
	GLuint _depthMapFBO;
	GLuint _depthMapTexture;
	std::uint64_t _depthShadowWidth = 1024;
	std::uint64_t _depthShadowHeight = 1024;

	// Screen quad framebuffer
	Shader* _screenquadShader = nullptr;
	GLuint _screenquadFBO;
	GLuint _screenquadRBO;
	GLuint _screenquadVAO;
	GLuint _screenquadVBO;
	GLuint _textureColorbuffer;
};

