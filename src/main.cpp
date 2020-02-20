#include "config.h"
#include "MainWindow.h"

#include <iostream>

#include <QApplication>
#include <QSurfaceFormat>

int main(int argc, char *argv[])
{
	std::cout << "#- " << PROJECT_NAME << " " << PROJECT_VER << " -#" << std::endl << std::endl;

	QApplication a(argc, argv);

	MainWindow w;
	w.show();

	return a.exec();
}
