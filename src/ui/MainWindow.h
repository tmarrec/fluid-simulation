#pragma once

#include <QtWidgets/QMainWindow>
#include <QOpenGLWidget>

#include "ui/ui_MainWindow.h"

namespace Ui {
	class MainWindow;
}

class MainWindow : public QMainWindow
{
Q_OBJECT
public:
	explicit MainWindow();
	~MainWindow();

private:
	Ui::MainWindow* _ui;

private slots:
};
