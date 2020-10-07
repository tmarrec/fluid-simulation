#include "MainWindow.h"


MainWindow::MainWindow(Renderer__ __renderer)
: QMainWindow{nullptr}
, _ui{new Ui::MainWindow}
{
	QSurfaceFormat format;
	format.setVersion(4, 6);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSamples(4);
	QSurfaceFormat::setDefaultFormat(format);

	_ui->setupUi(this);
	_QOpenGLWidget = _ui->openGLWidget;
	_QOpenGLWidget->setFocus();
	_QOpenGLWidget->setRenderer(__renderer);
}

MainWindow::~MainWindow()
{
	delete _ui;
}

void MainWindow::paint()
{
	_QOpenGLWidget->update();
}
