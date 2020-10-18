#include "Renderer.h"

#include "src/DrawableComponent.h"
#include "src/CameraComponent.h"
#include "src/LightComponent.h"

#include <GL/gl.h>
#include <memory>
#include <thread>

void Renderer::initGl() const
{
	resizeGl(0, 0);
}

void Renderer::resizeGl(int __w, int __h) const
{
	if (_activeCamera)
	{
		_activeCamera->setProjection(__w, __h);
	}
	glViewport(0, 0, __w, __h);
}

void Renderer::initDrawable(DrawableComponent* __drawableComponent)
{
	glGenVertexArrays(1, &__drawableComponent->VAO());

	glGenBuffers(1, &__drawableComponent->VBO());
	glGenBuffers(1, &__drawableComponent->NBO());
	glGenBuffers(1, &__drawableComponent->EBO());

	_sceneIndices += __drawableComponent->indices()->size();

	glBindVertexArray(__drawableComponent->VAO());
		glBindBuffer(GL_ARRAY_BUFFER, __drawableComponent->VBO());
		glBufferData(
			GL_ARRAY_BUFFER,
			__drawableComponent->vertices()->size()*sizeof(GLfloat),
			__drawableComponent->vertices()->data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, __drawableComponent->NBO());
		glBufferData(
			GL_ARRAY_BUFFER,
			__drawableComponent->normals()->size()*sizeof(GLfloat),
			__drawableComponent->normals()->data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, __drawableComponent->EBO());
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			__drawableComponent->indices()->size()*sizeof(GLuint),
			__drawableComponent->indices()->data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLuint), (GLvoid*)nullptr);
		glEnableVertexAttribArray(2);

	glBindVertexArray(0);
}

void Renderer::freeDrawable(DrawableComponent* __drawableComponent)
{
	glDeleteVertexArrays(1, &__drawableComponent->VAO());	
	glDeleteBuffers(1, &__drawableComponent->VBO());
	glDeleteBuffers(1, &__drawableComponent->NBO());
	glDeleteBuffers(1, &__drawableComponent->EBO());
}

void Renderer::_useShader(DrawableComponent* __drawableComponent)
{
	ASSERT(_activeCamera, "_activeCamera should not be nullptr");

	auto shader = __drawableComponent->shader();
	auto color = __drawableComponent->color();
	shader->use();
	shader->set_3f("_object_color", color);

	shader->set_mat4("model", __drawableComponent->getModel());
	shader->set_mat4("view", _activeCamera->getView());
	shader->set_mat4("projection", _activeCamera->projection());
	uint64_t time = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	float time_sin = sin(time/50000);
	shader->set_1f("time", time_sin);

	shader->set_3f("_view_pos", _activeCamera->getView()[3]);

	shader->set_1i("_light_nb", _lights.size());
	// Envoie les uniforms pour toutes les lumieres
	for (size_t i = 0; i < _lights.size(); ++i) {
		auto light = static_cast<LightComponent*>(_lights[i]);
		auto temp = std::string("_point_lights[") + std::to_string(i) + "].position";
		shader->set_3f(temp.c_str(), light->getPosition());
		temp = std::string("_point_lights[") + std::to_string(i) + "].color";
		shader->set_3f(temp.c_str(), light->color());
		temp = std::string("_point_lights[") + std::to_string(i) + "].intensity";
		shader->set_1f(temp.c_str(), light->intensity());
		temp = std::string("_point_lights[") + std::to_string(i) + "].constant";
		shader->set_1f(temp.c_str(), 1.0f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].linear";
		shader->set_1f(temp.c_str(), 0.0014f);
		temp = std::string("_point_lights[") + std::to_string(i) + "].quadratic";
		shader->set_1f(temp.c_str(), 0.000007f);
	}
}

void Renderer::clear() const
{
	glClearColor(0.85f, 0.9f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Renderer::drawDrawable(DrawableComponent* __drawableComponent)
{
	_useShader(__drawableComponent);
	glBindVertexArray(__drawableComponent->VAO());
	glDrawElements(__drawableComponent->drawMode(), __drawableComponent->indices()->size(), GL_UNSIGNED_INT, nullptr);
	glBindVertexArray(0);
}

void Renderer::setActiveCamera(CameraComponent* __cameraComponent)
{
	ASSERT(__cameraComponent, "__cameraComponent should not be nullptr");
	_activeCamera = __cameraComponent;	
}

void Renderer::addLight(LightComponent* __lightComponent)
{
	ASSERT(__lightComponent, "__lightComponent should not be nullptr");
	_lights.emplace_back(__lightComponent);
}
