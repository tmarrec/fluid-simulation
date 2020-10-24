#include <iostream>

#include <QApplication>

#include "config.h"
#include "ui/MainWindow.h"


int main(int argc, char *argv[])
{
	PRINT_TITLE()
	QApplication gui{argc, argv};
	MainWindow mainWindow {};
	mainWindow.show();
	return QApplication::exec();
}
