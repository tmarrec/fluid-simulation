#include "config.h"
#include "MainWindow.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QMenu>
#include <QMenuBar>
#include <QSlider>

#include <iostream>

MainWindow::MainWindow()
	: _selected_entity{NULL}
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
	QGroupBox *tree_box = new QGroupBox("Entities");
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

	QAction *add_light = new QAction("Light", this);
	connect(add_light, &QAction::triggered, _openGL, &OpenGL::add_light);
	add_menu->addAction(add_light);
	
	QAction *add_model = new QAction("Model", this);
	connect(add_model, &QAction::triggered, _openGL, &OpenGL::add_model);
	add_menu->addAction(add_model);



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

	// Position
	side_panel_l->addWidget(position_box());

	// Rotation
	side_panel_l->addWidget(rotation_box());

	// Scale
	side_panel_l->addWidget(scale_box());

	// Delete button
	QPushButton *delete_button = new QPushButton("Delete");
	connect(delete_button, &QPushButton::clicked, this, &MainWindow::delete_item_entities_tree_view);
    side_panel_l->addWidget(delete_button);

    properties_box->setLayout(side_panel_l);
	container->addWidget(properties_box);

	// Main Widget
    QWidget *w = new QWidget;
    w->setLayout(container);

	setCentralWidget(w);
}

MainWindow::~MainWindow() {

}

void MainWindow::delete_item_entities_tree_view() {
	std::cout << "kek" << std::endl;
}	

// TRANSLATION
void MainWindow::change_slide_x_position(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{value, p.y, p.z});
		_slide_x_position_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_y_position(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, value, p.z});
		_slide_y_position_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_z_position(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, p.y, value});
		_slide_z_position_label->setText(std::to_string(value).c_str());
	}
}

// ROTATIONS
void MainWindow::change_slide_x_rotation(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{value, r.y, r.z});
		_slide_x_rotation_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_y_rotation(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, value, r.z});
		_slide_y_rotation_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_z_rotation(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, r.y, value});
		_slide_z_rotation_label->setText(std::to_string(value).c_str());
	}
}

// SCALE
void MainWindow::change_slide_x_scale(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{value, r.y, r.z});
		_slide_x_scale_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_y_scale(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, value, r.z});
		_slide_y_scale_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::change_slide_z_scale(int value) {
	if (_selected_entity != NULL) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, r.y, value});
		_slide_z_scale_label->setText(std::to_string(value).c_str());
	}
}

void MainWindow::update_slide_position(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != NULL && _selected_entity->shape_ptr()->id() == id) {
		_slide_x_position->setValue(pos.x);
		_slide_x_position_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_position->setValue(pos.y);
		_slide_y_position_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_position->setValue(pos.z);
		_slide_z_position_label->setText(std::to_string(int(pos.z)).c_str());
	}
}

void MainWindow::update_slide_rotation(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != NULL && _selected_entity->shape_ptr()->id() == id) {
		_slide_x_rotation->setValue(pos.x);
		_slide_x_rotation_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_rotation->setValue(pos.y);
		_slide_y_rotation_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_rotation->setValue(pos.z);
		_slide_z_rotation_label->setText(std::to_string(int(pos.z)).c_str());
	}
}

void MainWindow::update_slide_scale(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != NULL && _selected_entity->shape_ptr()->id() == id) {
		_slide_x_scale->setValue(pos.x);
		_slide_x_scale_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_scale->setValue(pos.y);
		_slide_y_scale_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_scale->setValue(pos.z);
		_slide_z_scale_label->setText(std::to_string(int(pos.z)).c_str());
	}
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

QGroupBox* MainWindow::position_box() {
	QGroupBox *box = new QGroupBox("Position");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *slide_x_box = new QGroupBox();
	QHBoxLayout *slide_x_layout = new QHBoxLayout;
	_slide_x_position_label = new QLabel("x");
	_slide_x_position = new QSlider(Qt::Orientation::Horizontal);
	_slide_x_position->setMinimumWidth(180);
	_slide_x_position->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_x_position->setMinimum(-100);
	_slide_x_position->setMaximum(100);
	slide_x_layout->addWidget(_slide_x_position);
	slide_x_layout->addWidget(_slide_x_position_label);
	slide_x_box->setLayout(slide_x_layout);
	box_layout->addWidget(slide_x_box);
	connect(_slide_x_position, &QSlider::sliderMoved, this, &MainWindow::change_slide_x_position);

	QGroupBox *slide_y_box = new QGroupBox();
	QHBoxLayout *slide_y_layout = new QHBoxLayout;
	_slide_y_position_label = new QLabel("y");
	_slide_y_position = new QSlider(Qt::Orientation::Horizontal);
	_slide_y_position->setMinimumWidth(180);
	_slide_y_position->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_y_position->setMinimum(-100);
	_slide_y_position->setMaximum(100);
	slide_y_layout->addWidget(_slide_y_position);
	slide_y_layout->addWidget(_slide_y_position_label);
	slide_y_box->setLayout(slide_y_layout);
	box_layout->addWidget(slide_y_box);
	connect(_slide_y_position, &QSlider::sliderMoved, this, &MainWindow::change_slide_y_position);

	QGroupBox *slide_z_box = new QGroupBox();
	QHBoxLayout *slide_z_layout = new QHBoxLayout;
	_slide_z_position_label = new QLabel("z");
	_slide_z_position = new QSlider(Qt::Orientation::Horizontal);
	_slide_z_position->setMinimumWidth(180);
	_slide_z_position->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_z_position->setMinimum(-100);
	_slide_z_position->setMaximum(100);
	slide_z_layout->addWidget(_slide_z_position);
	slide_z_layout->addWidget(_slide_z_position_label);
	slide_z_box->setLayout(slide_z_layout);
	box_layout->addWidget(slide_z_box);
	connect(_slide_z_position, &QSlider::sliderMoved, this, &MainWindow::change_slide_z_position);
	
	box->setLayout(box_layout);
	return box;
}

QGroupBox* MainWindow::rotation_box() {
	QGroupBox *box = new QGroupBox("Rotation");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *slide_x_box = new QGroupBox();
	QHBoxLayout *slide_x_layout = new QHBoxLayout;
	_slide_x_rotation_label = new QLabel("x");
	_slide_x_rotation = new QSlider(Qt::Orientation::Horizontal);
	_slide_x_rotation->setMinimumWidth(180);
	_slide_x_rotation->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_x_rotation->setMinimum(0);
	_slide_x_rotation->setMaximum(360);
	slide_x_layout->addWidget(_slide_x_rotation);
	slide_x_layout->addWidget(_slide_x_rotation_label);
	slide_x_box->setLayout(slide_x_layout);
	box_layout->addWidget(slide_x_box);
	connect(_slide_x_rotation, &QSlider::sliderMoved, this, &MainWindow::change_slide_x_rotation);

	QGroupBox *slide_y_box = new QGroupBox();
	QHBoxLayout *slide_y_layout = new QHBoxLayout;
	_slide_y_rotation_label = new QLabel("y");
	_slide_y_rotation = new QSlider(Qt::Orientation::Horizontal);
	_slide_y_rotation->setMinimumWidth(180);
	_slide_y_rotation->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_y_rotation->setMinimum(0);
	_slide_y_rotation->setMaximum(360);
	slide_y_layout->addWidget(_slide_y_rotation);
	slide_y_layout->addWidget(_slide_y_rotation_label);
	slide_y_box->setLayout(slide_y_layout);
	box_layout->addWidget(slide_y_box);
	connect(_slide_y_rotation, &QSlider::sliderMoved, this, &MainWindow::change_slide_y_rotation);

	QGroupBox *slide_z_box = new QGroupBox();
	QHBoxLayout *slide_z_layout = new QHBoxLayout;
	_slide_z_rotation_label = new QLabel("z");
	_slide_z_rotation = new QSlider(Qt::Orientation::Horizontal);
	_slide_z_rotation->setMinimumWidth(180);
	_slide_z_rotation->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_z_rotation->setMinimum(0);
	_slide_z_rotation->setMaximum(360);
	slide_z_layout->addWidget(_slide_z_rotation);
	slide_z_layout->addWidget(_slide_z_rotation_label);
	slide_z_box->setLayout(slide_z_layout);
	box_layout->addWidget(slide_z_box);
	connect(_slide_z_rotation, &QSlider::sliderMoved, this, &MainWindow::change_slide_z_rotation);
	
	box->setLayout(box_layout);
	return box;
}

QGroupBox* MainWindow::scale_box() {
	QGroupBox *box = new QGroupBox("Scale");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *slide_x_box = new QGroupBox();
	QHBoxLayout *slide_x_layout = new QHBoxLayout;
	_slide_x_scale_label = new QLabel("x");
	_slide_x_scale = new QSlider(Qt::Orientation::Horizontal);
	_slide_x_scale->setMinimumWidth(180);
	_slide_x_scale->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_x_scale->setMinimum(0);
	_slide_x_scale->setMaximum(50);
	slide_x_layout->addWidget(_slide_x_scale);
	slide_x_layout->addWidget(_slide_x_scale_label);
	slide_x_box->setLayout(slide_x_layout);
	box_layout->addWidget(slide_x_box);
	connect(_slide_x_scale, &QSlider::sliderMoved, this, &MainWindow::change_slide_x_scale);

	QGroupBox *slide_y_box = new QGroupBox();
	QHBoxLayout *slide_y_layout = new QHBoxLayout;
	_slide_y_scale_label = new QLabel("y");
	_slide_y_scale = new QSlider(Qt::Orientation::Horizontal);
	_slide_y_scale->setMinimumWidth(180);
	_slide_y_scale->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_y_scale->setMinimum(0);
	_slide_y_scale->setMaximum(50);
	slide_y_layout->addWidget(_slide_y_scale);
	slide_y_layout->addWidget(_slide_y_scale_label);
	slide_y_box->setLayout(slide_y_layout);
	box_layout->addWidget(slide_y_box);
	connect(_slide_y_scale, &QSlider::sliderMoved, this, &MainWindow::change_slide_y_scale);

	QGroupBox *slide_z_box = new QGroupBox();
	QHBoxLayout *slide_z_layout = new QHBoxLayout;
	_slide_z_scale_label = new QLabel("z");
	_slide_z_scale = new QSlider(Qt::Orientation::Horizontal);
	_slide_z_scale->setMinimumWidth(180);
	_slide_z_scale->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_z_scale->setMinimum(0);
	_slide_z_scale->setMaximum(50);
	slide_z_layout->addWidget(_slide_z_scale);
	slide_z_layout->addWidget(_slide_z_scale_label);
	slide_z_box->setLayout(slide_z_layout);
	box_layout->addWidget(slide_z_box);
	connect(_slide_z_scale, &QSlider::sliderMoved, this, &MainWindow::change_slide_z_scale);
	
	box->setLayout(box_layout);
	return box;
}
