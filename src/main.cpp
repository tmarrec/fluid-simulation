#include "MainWindow.h"

#include <iostream>

#include <QApplication>
#include <QSurfaceFormat>

int main(int argc, char *argv[])
{
	std::cout << "#- Cowboy Engine 0.1 -#" << std::endl;

	QApplication a(argc, argv);

	QSurfaceFormat fmt;
	fmt.setSamples(24);
	QSurfaceFormat::setDefaultFormat(fmt);

	MainWindow w;
	w.show();

	return a.exec();
}
