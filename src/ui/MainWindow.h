#pragma once

#include <memory>

#include <QtWidgets/QMainWindow>
#include <QOpenGLWidget>

#include "../utils.h"
#include "../Renderer.h"
#include "GlWidget.h"
#include "ui/ui_MainWindow.h"

namespace Ui {
	class MainWindow;
}

class MainWindow : public QMainWindow
{
Q_OBJECT

public:
	explicit MainWindow(Renderer__ __renderer);
	void paint();

	~MainWindow();

private:
	Ui::MainWindow* _ui;
	GlWidget* _QOpenGLWidget;

private slots:

};
