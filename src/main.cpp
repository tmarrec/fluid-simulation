#include "config.h"
#include "MainWindow.h"

#include <iostream>

#include <QApplication>
#include <QSurfaceFormat>

int main(int argc, char *argv[])
{
	std::cout << "#- " << PROJECT_NAME << " " << PROJECT_VER << " -#" << std::endl << std::endl;

	QApplication a(argc, argv);

	QSurfaceFormat fmt;
	fmt.setSamples(24);
	QSurfaceFormat::setDefaultFormat(fmt);

	MainWindow w;
	w.show();

	return a.exec();
}
