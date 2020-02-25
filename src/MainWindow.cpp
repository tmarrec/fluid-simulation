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
	
	// Parametres OpenGL 
	QSurfaceFormat format;
	format.setVersion(4, 5);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSwapInterval(0);
	format.setSamples(8);
	QSurfaceFormat::setDefaultFormat(format);
	
	// Widget OpenGL
	_glw = new GLWidget(this);
	_openGL = _glw->openGL();


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
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	connect(add_triangle, &QAction::triggered, _openGL, &OpenGL::add_triangle);
	add_menu->addAction(add_triangle);

	QAction *add_cube = new QAction("Cube", this);
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	connect(add_cube, &QAction::triggered, _openGL, &OpenGL::add_cube);
	add_menu->addAction(add_cube);

	QAction *add_sphere = new QAction("Sphere", this);
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
	//TODO connect a ECS directement putot que openGL
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
	auto s = _selected_entity->shape_ptr();
	auto p = s->position();
	s->set_position(glm::vec3{0.05f*value, p.y, p.z});
}

void MainWindow::change_slide_y(int value) {
	auto s = _selected_entity->shape_ptr();
	auto p = s->position();
	s->set_position(glm::vec3{p.x, 0.05f*value, p.z});
}

void MainWindow::change_slide_z(int value) {
	auto s = _selected_entity->shape_ptr();
	auto p = s->position();
	s->set_position(glm::vec3{p.x, p.y, 0.05f*value});
}

void MainWindow::on_item_clicked(QListWidgetItem *item) {
	auto s = item->data(100).value<Entity_Item*>();
	change_selected_entity(s);
}

void MainWindow::change_selected_entity(Entity_Item* e) {
	_selected_entity = e;
}



Q_DECLARE_METATYPE(Entity_Item*)

void MainWindow::add_item_to_QListW(std::shared_ptr<Entity> shape_ptr) {
	auto s = new Entity_Item(shape_ptr);
	change_selected_entity(s);
	QListWidgetItem *item = new QListWidgetItem(s->name().c_str());
	item->setData(100, QVariant::fromValue(s));
	_list->addItem(item);
}
