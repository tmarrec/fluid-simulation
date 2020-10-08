#include "MainWindow.h"

MainWindow::MainWindow()
: QMainWindow{nullptr}
, _ui{new Ui::MainWindow}
{
	std::cout << ":)" << std::endl;
	QSurfaceFormat format;
	format.setVersion(4, 6);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSamples(4);
	QSurfaceFormat::setDefaultFormat(format);

	_ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete _ui;
}
