#include "MainWindow.h"
#include "GLWidget.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QMenu>
#include <QMenuBar>


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
	
	// Menu
	QWidget *topFiller = new QWidget;
    topFiller->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	container->addWidget(topFiller);
	QMenu *add_menu = menuBar()->addMenu("Add");
	
	// Menu Actions
	QAction *add_triangle = new QAction("Triangle", this);
    connect(add_triangle, &QAction::triggered, this, &MainWindow::add_triangle);
	add_menu->addAction(add_triangle);

	QAction *add_cube = new QAction("Cube", this);
    connect(add_cube, &QAction::triggered, this, &MainWindow::add_cube);
	add_menu->addAction(add_cube);

	QAction *add_sphere = new QAction("Sphere", this);
    connect(add_sphere, &QAction::triggered, this, &MainWindow::add_sphere);
	add_menu->addAction(add_sphere);

	// Left Panel, with TreeView
	QVBoxLayout *left_side_panel_l = new QVBoxLayout;	
	_list = new QListWidget;
    left_side_panel_l->addWidget(_list);
	QWidget *left_side_panel_w = new QWidget;
    left_side_panel_w->setLayout(left_side_panel_l);
	container->addWidget(left_side_panel_w);

	connect(_list, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(on_item_clicked(QListWidgetItem*)));

	// Widget OpenGL
	_glw = new GLWidget(this);
    container->addWidget(_glw);

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

void MainWindow::on_item_clicked(QListWidgetItem *item) {
	uint selected_id = item->data(100).value<uint>();
	std::cout << selected_id << std::endl;
}

void MainWindow::on_triangleButton_clicked()
{
	std::cout << "triangle" << std::endl;
}

void MainWindow::on_rectangleButton_clicked()
{
	std::cout << "rect" << std::endl;
}

void MainWindow::add_triangle() {
	_glw->add_shape("triangle");
}

void MainWindow::add_cube() {
	_glw->add_shape("cube");
}

void MainWindow::add_sphere() {
	_glw->add_shape("sphere");
}

void MainWindow::add_item_to_QListW(uint id, std::string name) {
	QListWidgetItem *item = new QListWidgetItem(name.c_str());
	item->setData(100, QVariant(id));
	_list->addItem(item);
}
