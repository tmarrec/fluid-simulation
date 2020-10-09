#include "MainWindow.h"

MainWindow::MainWindow()
: QMainWindow{nullptr}
, _ui{new Ui::MainWindow}
{
	QSurfaceFormat format;
	format.setVersion(4, 6);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSamples(4);
	format.setSwapInterval(0);
	QSurfaceFormat::setDefaultFormat(format);

	_ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete _ui;
}
