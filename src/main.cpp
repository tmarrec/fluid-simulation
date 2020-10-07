#include <iostream>

#include <GL/gl.h>
#include <QApplication>

#include "config.h"
#include "Renderer.h"
#include "ui/MainWindow.h"
#include "ECS.h"
#include "DrawableComponent.h"
#include "CameraComponent.h"
#include "TransformComponent.h"

using Renderer__ = std::shared_ptr<Renderer>; 

int main(int argc, char *argv[])
{
	PRINT_TITLE()

	QApplication gui{argc, argv};

	Renderer renderer {};
	Renderer__ renderer_ = std::make_shared<Renderer>(renderer);

	Manager manager {};

	MainWindow mainWindow {renderer_};
	mainWindow.show();

	std::vector<GLfloat> v = std::vector<GLfloat>{
					-0.5f, -0.5f, 0.5f,
					0.5f, -0.5f, 0.5f,
					0.5f, 0.5f, 0.5f,
					-0.5f, 0.5f, 0.5f,

					0.5f, 0.5f, 0.5f,
					0.5f, 0.5f, -0.5f,
					0.5f, -0.5f, -0.5f,
					0.5f, -0.5f, 0.5f,

					-0.5f, -0.5f, -0.5f,
					0.5f, -0.5f, -0.5f,
					0.5f, 0.5f, -0.5f,
					-0.5f, 0.5f, -0.5f,

					-0.5f, -0.5f, -0.5f,
					-0.5f, -0.5f, 0.5f,
					-0.5f, 0.5f, 0.5f,
					-0.5f, 0.5f, -0.5f,
						
					0.5f, 0.5f, 0.5f,
					-0.5f, 0.5f, 0.5f,
					-0.5f, 0.5f, -0.5f,
					0.5f, 0.5f, -0.5f,

					-0.5f, -0.5f, -0.5f,
					0.5f, -0.5f, -0.5f,
					0.5f, -0.5f, 0.5f,
					-0.5f, -0.5f, 0.5f,
	};
	std::vector<GLfloat> n = std::vector<GLfloat>{
					0.0f, 0.0f, 1.0f,
					0.0f, 0.0f, 1.0f,
					0.0f, 0.0f, 1.0f,
					0.0f, 0.0f, 1.0f,

					1.0f, 0.0f, 0.0f,
					1.0f, 0.0f, 0.0f,
					1.0f, 0.0f, 0.0f,
					1.0f, 0.0f, 0.0f,

					0.0f, 0.0f, -1.0f,
					0.0f, 0.0f, -1.0f,
					0.0f, 0.0f, -1.0f,
					0.0f, 0.0f, -1.0f,

					-1.0f, 0.0f, 0.0f,
					-1.0f, 0.0f, 0.0f,
					-1.0f, 0.0f, 0.0f,
					-1.0f, 0.0f, 0.0f,

					0.0f, 1.0f, 0.0f,
					0.0f, 1.0f, 0.0f,
					0.0f, 1.0f, 0.0f,
					0.0f, 1.0f, 0.0f,

					0.0f, -1.0f, 0.0f,
					0.0f, -1.0f, 0.0f,
					0.0f, -1.0f, 0.0f,
					0.0f, -1.0f, 0.0f,
				};
	std::vector<GLuint> i = std::vector<GLuint>{
					0,  1,  2,  0,  2,  3,   //front
					4,  5,  6,  4,  6,  7,   //right
					8,  9,  10, 8,  10, 11,  //back
					12, 13, 14, 12, 14, 15,  //left
					16, 17, 18, 16, 18, 19,  //upper
					20, 21, 22, 20, 22, 23 	 //bottom
				};

	auto & cube(manager.addEntity());
	cube.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
		glm::vec3{50.0f, 50.0f, 50.0f});
	cube.addComponent<DrawableComponent>(renderer_, "shaders/vert.vert", "shaders/frag.frag", std::make_shared<std::vector<GLfloat>>(v), std::make_shared<std::vector<GLfloat>>(n), std::make_shared<std::vector<GLuint>>(i));


	/*
	auto & cube2(manager.addEntity());
	cube2.addComponent<TransformComponent>(glm::vec3{5.0f, 0.0f, 0.0f}, glm::vec3{3.0f, 0.0f, 0.0f},
		glm::vec3{4.0f, 5.0f, 5.0f});
	cube2.addComponent<DrawableComponent>(renderer_, "shaders/vert.vert", "shaders/frag.frag", std::make_shared<std::vector<GLfloat>>(v), std::make_shared<std::vector<GLfloat>>(n), std::make_shared<std::vector<GLuint>>(i));
	*/

	auto & camera(manager.addEntity());
	camera.addComponent<CameraComponent>(0.0f, 0.0f, 15.0f, 90.0f);
	camera.addComponent<TransformComponent>(glm::vec3{-250.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
		glm::vec3{1.0f, 1.0f, 1.0f});
	
	std::cout << "start draw" << std::endl;

	// Update in the game loop
	for (int i = 0; i < 1; ++i)
	{
		renderer_->setActiveCamera(&camera.getComponent<CameraComponent>());
		manager.update();
		renderer_->draw();
	}
	
	QApplication::exec();

	return EXIT_SUCCESS;
}
