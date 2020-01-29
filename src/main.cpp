#include "mainwindow.h"

#include <iostream>

#include <QApplication>

int main(int argc, char *argv[])
{
	std::cout << "#- Cowboy Engine 0.1 -#" << std::endl;

	QApplication a(argc, argv);
	MainWindow w;
	w.show();
	return a.exec();
}
