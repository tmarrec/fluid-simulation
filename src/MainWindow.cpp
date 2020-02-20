#include "config.h"
#include "MainWindow.h"
#include "GLWidget.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QMenu>
#include <QMenuBar>
#include <QSlider>
#include <QGroupBox>


#include <iostream>

MainWindow::MainWindow()
	: _selected_entity{0}
{
	setWindowTitle(TITLE);
	
	// Widget OpenGL
	_glw = new GLWidget(this);
	_openGL = _glw->openGL();

	// Parametres OpenGL 
	QSurfaceFormat format;
	format.setVersion(4, 5);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSwapInterval(0);
	format.setSamples(8);
	QSurfaceFormat::setDefaultFormat(format);

    // Layout QT
    QHBoxLayout *container = new QHBoxLayout;

	QGroupBox *properties_box = new QGroupBox("Properties");
	properties_box->setMinimumWidth(300);
	QGroupBox *tree_box = new QGroupBox("Tree");
	tree_box->setMinimumWidth(200);
	
	// Menu
	QWidget *topFiller = new QWidget;
    topFiller->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
	container->addWidget(topFiller);
	QMenu *add_menu = menuBar()->addMenu("Add");
	
	// Menu Actions
	QAction *add_triangle = new QAction("Triangle", this);
	connect(add_triangle, &QAction::triggered, _openGL, &OpenGL::add_triangle);
	add_menu->addAction(add_triangle);

	QAction *add_cube = new QAction("Cube", this);
	connect(add_cube, &QAction::triggered, _openGL, &OpenGL::add_cube);
	add_menu->addAction(add_cube);

	QAction *add_sphere = new QAction("Sphere", this);
	connect(add_sphere, &QAction::triggered, _openGL, &OpenGL::add_sphere);
	add_menu->addAction(add_sphere);

	// Left Panel, with TreeView
	QVBoxLayout *left_side_panel_l = new QVBoxLayout;	
	_list = new QListWidget;
    left_side_panel_l->addWidget(_list);
    tree_box->setLayout(left_side_panel_l);
	container->addWidget(tree_box);
	connect(_list, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(on_item_clicked(QListWidgetItem*)));

    container->addWidget(_glw);

	// Side Panel
	QVBoxLayout *side_panel_l = new QVBoxLayout;	
	QPushButton *test = new QPushButton;
    side_panel_l->addWidget(test);
	QPushButton *test2 = new QPushButton;
    side_panel_l->addWidget(test2);


	QSlider *slide_x = new QSlider(Qt::Orientation::Horizontal);
	slide_x->setMinimum(-100);
	slide_x->setMaximum(100);
    side_panel_l->addWidget(slide_x);
	connect(slide_x, &QSlider::valueChanged, this, &MainWindow::change_slide_x);
	QSlider *slide_y = new QSlider(Qt::Orientation::Horizontal);
	slide_y->setMinimum(-100);
	slide_y->setMaximum(100);
    side_panel_l->addWidget(slide_y);
	connect(slide_y, &QSlider::valueChanged, this, &MainWindow::change_slide_y);
	QSlider *slide_z = new QSlider(Qt::Orientation::Horizontal);
	slide_z->setMinimum(-100);
	slide_z->setMaximum(100);
    side_panel_l->addWidget(slide_z);
	connect(slide_z, &QSlider::valueChanged, this, &MainWindow::change_slide_z);

    properties_box->setLayout(side_panel_l);
	container->addWidget(properties_box);

	// Main Widget
    QWidget *w = new QWidget;
    w->setLayout(container);

	setCentralWidget(w);
}

MainWindow::~MainWindow() {
}

void MainWindow::change_slide_x(int value) {
	_glw->move(_selected_entity, 'x', value);
}

void MainWindow::change_slide_y(int value) {
	_glw->move(_selected_entity, 'y', value);
}

void MainWindow::change_slide_z(int value) {
	_glw->move(_selected_entity, 'z', value);
}

void MainWindow::on_item_clicked(QListWidgetItem *item) {
	uint selected_id = item->data(100).value<uint>();
	_selected_entity = selected_id;
}

void MainWindow::on_triangleButton_clicked()
{
	std::cout << "triangle" << std::endl;
}

void MainWindow::on_rectangleButton_clicked()
{
	std::cout << "rect" << std::endl;
}

void MainWindow::add_item_to_QListW(uint id, std::string name) {
	QListWidgetItem *item = new QListWidgetItem(name.c_str());
	item->setData(100, QVariant(id));
	_list->addItem(item);
}
