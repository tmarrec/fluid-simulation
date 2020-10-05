#include "config.h"

#include "MessageBus.h"
#include "Renderer.h"
#include "ui/MainWindow.h"

#include "ECS.h"
#include "DrawableComponent.h"
#include "CameraComponent.h"
#include "TransformComponent.h"

#include <iostream>

#include <QApplication>

void printTitle()
{
	std::cout << "\033[96m\033[1m";
	std::cout << "                          __                                                                            \n                         /\\ \\                                                      __                   \n  ___    ___   __  __  __\\ \\ \\____    ___   __  __              __    ___      __ /\\_\\    ___      __   \n /'___\\ / __`\\/\\ \\/\\ \\/\\ \\\\ \\ '__`\\  / __`\\/\\ \\/\\ \\  _______  /'__`\\/' _ `\\  /'_ `\\/\\ \\ /' _ `\\  /'__`\\ \n/\\ \\__//\\ \\L\\ \\ \\ \\_/ \\_/ \\\\ \\ \\L\\ \\/\\ \\L\\ \\ \\ \\_\\ \\/\\______\\/\\  __//\\ \\/\\ \\/\\ \\L\\ \\ \\ \\/\\ \\/\\ \\/\\  __/ \n\\ \\____\\ \\____/\\ \\___x___/' \\ \\_,__/\\ \\____/\\/`____ \\/______/\\ \\____\\ \\_\\ \\_\\ \\____ \\ \\_\\ \\_\\ \\_\\ \\____\\\n \\/____/\\/___/  \\/__//__/    \\/___/  \\/___/  `/___/> \\        \\/____/\\/_/\\/_/\\/___L\\ \\/_/\\/_/\\/_/\\/____/\n                                                /\\___/                         /\\____/                  \n                                                \\/__/                          \\_/__/                   \n" << std::endl;
	std::cout << "\033[39m\033[0m";
}

int main(int argc, char *argv[])
{
	printTitle();
	MessageBus * msgBus_ptr = new MessageBus();
	MsgBus_ptr messageBus {msgBus_ptr};
	Renderer renderer {messageBus};

	QApplication gui{argc, argv};
	MainWindow mainWindow {messageBus};
	mainWindow.show();

	// Testings
	// ECS Manager
	

	// TODO Maybe link Entity to mainbus instead of each component...
	Manager manager {messageBus};

	auto & cube(manager.addEntity(messageBus));
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

	cube.addComponent<TransformComponent>(glm::vec3{0.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
		glm::vec3{1.0f, 1.0f, 1.0f});
	cube.addComponent<DrawableComponent>("shaders/vert.vert", "shaders/frag.frag", v, n, i);

	auto & camera(manager.addEntity(messageBus));
	camera.addComponent<CameraComponent>(0.0f, 0.0f, 15.0f, 90.0f);
	camera.addComponent<TransformComponent>(glm::vec3{-250.0f, 0.0f, 0.0f}, glm::vec3{0.0f, 0.0f, 0.0f},
		glm::vec3{1.0f, 1.0f, 1.0f});

	// Update in the game loop
	manager.update();
	
	QApplication::exec();

	return EXIT_SUCCESS;
}
