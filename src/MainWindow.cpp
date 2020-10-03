#include "MainWindow.h"
#include "ui/ui_MainWindow.h"

MainWindow::MainWindow(MessageBus & messageBus)
: QMainWindow{nullptr}
, System{messageBus}
, ui{new Ui::MainWindow}
{

	QSurfaceFormat format;
	format.setVersion(4, 6);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSamples(8);
	QSurfaceFormat::setDefaultFormat(format);

	ui->setupUi(this);
	ui->openGLWidget->setFocus();

	Message helloMsg {HELLO, this};
	postMessage(helloMsg);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::cout(std::string string) const
{
	std::cout << "0x" << std::hex << std::this_thread::get_id() << " ";
	std::cout << "\033[43m\033[1m";
	std::cout << "[MainWindow]";
	std::cout << "\033[49m\033[0m";
	std::cout << " " << string << std::endl;
}

