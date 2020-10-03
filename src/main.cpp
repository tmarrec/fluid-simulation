#include "config.h"

#include "MessageBus.h"
#include "Renderer.h"
#include "ui/MainWindow.h"

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
	MessageBus messageBus = MessageBus();
	Renderer renderer {messageBus};

	QApplication gui{argc, argv};
	MainWindow mainWindow {messageBus};
	mainWindow.show();
	
	QApplication::exec();

	return EXIT_SUCCESS;
}
