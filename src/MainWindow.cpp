#include "MainWindow.h"
#include "GLWidget.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>


#include <iostream>

MainWindow::MainWindow() {
	setWindowTitle(tr("cowboy-engine 0.1"));

	// Parametres OpenGL 
	QSurfaceFormat format;
	format.setVersion(4, 5);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSwapInterval(0);
	QSurfaceFormat::setDefaultFormat(format);

    // Layout QT
	QVBoxLayout *main_layout = new QVBoxLayout;
    QHBoxLayout *container = new QHBoxLayout;


	// Widget OpenGL
	GLWidget *glw = new GLWidget(this);
    container->addWidget(glw);

	// Side Panel
	QVBoxLayout *side_panel_l = new QVBoxLayout;	
	QPushButton *test = new QPushButton;
    side_panel_l->addWidget(test);
	QPushButton *test2 = new QPushButton;
    side_panel_l->addWidget(test2);

    QWidget *side_panel_w = new QWidget;
    side_panel_w->setLayout(side_panel_l);
	container->addWidget(side_panel_w);

	// Main Widget
    QWidget *w = new QWidget;
    w->setLayout(container);

	
    main_layout->addWidget(w);
    setLayout(main_layout);

	setCentralWidget(w);
}

MainWindow::~MainWindow() {
}


void MainWindow::on_triangleButton_clicked()
{
	std::cout << "triangle" << std::endl;
}

void MainWindow::on_rectangleButton_clicked()
{
	std::cout << "rect" << std::endl;
}

