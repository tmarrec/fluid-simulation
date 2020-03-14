#include "config.h"
#include "MainWindow.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenu>
#include <QMenuBar>
#include <QKeyEvent>
#include <QMouseEvent>

#include <iostream>

MainWindow::MainWindow()
	: _selected_entity{nullptr}
	, _last_mouse_x{0}
	, _last_mouse_y{0}
{
	setWindowTitle(TITLE);

	setFocusPolicy(Qt::ClickFocus);
	
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
	add_menu->addAction(add_triangle);

	QAction *add_cube = new QAction("Cube", this);
	add_menu->addAction(add_cube);

	QAction *add_sphere = new QAction("Sphere", this);
	add_menu->addAction(add_sphere);

	QAction *add_light = new QAction("Light", this);
	add_menu->addAction(add_light);
	
	QAction *add_model = new QAction("Model", this);
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
	_delete_button = new QPushButton("Delete");
    side_panel_l->addWidget(_delete_button);

    properties_box->setLayout(side_panel_l);
	container->addWidget(properties_box);

	// Main Widget
    QWidget *w = new QWidget;
    w->setLayout(container);

	setCentralWidget(w);

	_glw->init();
	_openGL = _glw->openGL();

	connect(_delete_button, &QPushButton::clicked, this, &MainWindow::delete_item_entities_tree_view);

	connect(add_triangle, &QAction::triggered, _openGL, &OpenGL::add_triangle);
	connect(add_cube, &QAction::triggered, _openGL, &OpenGL::add_cube);
	connect(add_sphere, &QAction::triggered, _openGL, &OpenGL::add_sphere);
	connect(add_light, &QAction::triggered, _openGL, &OpenGL::add_light);
	connect(add_model, &QAction::triggered, _openGL, &OpenGL::add_model);

}

MainWindow::~MainWindow() {

}

void MainWindow::mousePressEvent(QMouseEvent *event) {
	_last_mouse_x = event->x();
	_last_mouse_y = event->y();
}

void MainWindow::mouseMoveEvent(QMouseEvent *event) {
	int x = event->x();
	int y = event->y();

	auto camera = std::static_pointer_cast<Camera>(_camera);
	float yaw = camera->yaw();
	float pitch = camera->pitch();

	float sensitivity = 0.3f;
	float offset_x = (x-_last_mouse_x)*sensitivity;
	float offset_y = (y-_last_mouse_y)*sensitivity;

	yaw += offset_x;
	pitch -= offset_y;
	pitch = glm::clamp(pitch, -89.9f, 89.9f);

	camera->set_yaw(yaw);
	camera->set_pitch(pitch);

	glm::vec3 dir;
	dir.x = cos(glm::radians(yaw))*cos(glm::radians(pitch));
	dir.y = sin(glm::radians(pitch));
	dir.z = sin(glm::radians(yaw))*cos(glm::radians(pitch));
	dir = glm::normalize(dir);

	camera->set_front(dir);

	_last_mouse_x = x;
	_last_mouse_y = y;
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
	auto key = event->key();
	auto p = _camera->position();
	auto camera = std::static_pointer_cast<Camera>(_camera);
	auto front = camera->front();
	auto up = camera->up();
	float speed = 5.0f;
	switch(key) {
		case Qt::Key_Z:
			_camera->set_position(p+(front*speed));
			break;
		case Qt::Key_Q:
			_camera->set_position(p-(glm::normalize(glm::cross(front, up))*speed));
			break;
		case Qt::Key_S:
			_camera->set_position(p-(front*speed));
			break;
		case Qt::Key_D:
			_camera->set_position(p+(glm::normalize(glm::cross(front, up))*speed));
			break;
		case Qt::Key_Shift:
			_camera->set_position(glm::vec3{p.x, p.y-speed, p.z});
			break;
		case Qt::Key_Control:
			_camera->set_position(glm::vec3{p.x, p.y+speed, p.z});
			break;
	}
}

void MainWindow::delete_item_entities_tree_view() {
	_openGL->remove_entity(_selected_entity->shape_ptr());
	_list->takeItem(_list->currentRow());

	_list->setCurrentRow(_list->count()-1);
	auto new_item = _list->item(_list->count()-1);
	auto s = new_item->data(100).value<Entity_Item*>();
	change_selected_entity(s);
}	

// TRANSLATION
void MainWindow::change_slide_x_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{value, p.y, p.z});
	}
}

void MainWindow::change_slide_y_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, value, p.z});
	}
}

void MainWindow::change_slide_z_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, p.y, value});
	}
}

// ROTATIONS
void MainWindow::change_slide_x_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{value, r.y, r.z});
	}
}

void MainWindow::change_slide_y_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, value, r.z});
	}
}

void MainWindow::change_slide_z_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, r.y, value});
	}
}

// SCALE
void MainWindow::change_slide_x_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{value, r.y, r.z});
	}
}

void MainWindow::change_slide_y_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, value, r.z});
	}
}

void MainWindow::change_slide_z_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->shape_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, r.y, value});
	}
}

void MainWindow::update_slide_position(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->shape_ptr()->id() == id) {
		_slide_x_position->setValue(pos.x);
		_slide_x_position_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_position->setValue(pos.y);
		_slide_y_position_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_position->setValue(pos.z);
		_slide_z_position_label->setText(std::to_string(int(pos.z)).c_str());
	}
}

void MainWindow::update_slide_rotation(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->shape_ptr()->id() == id) {
		_slide_x_rotation->setValue(pos.x);
		_slide_x_rotation_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_rotation->setValue(pos.y);
		_slide_y_rotation_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_rotation->setValue(pos.z);
		_slide_z_rotation_label->setText(std::to_string(int(pos.z)).c_str());
	}
}

void MainWindow::update_slide_scale(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->shape_ptr()->id() == id) {
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
	auto s = e->shape_ptr();
	update_slide_position(s->position(), s->id());
	update_slide_rotation(s->rotation(), s->id());
	update_slide_scale(s->scale(), s->id());
	if (s->name() == "Camera 0") {
		_delete_button->setEnabled(false);
	} else {
		_delete_button->setEnabled(true);
	}
}

Q_DECLARE_METATYPE(Entity_Item*)

void MainWindow::add_item_to_QListW(std::shared_ptr<Entity> shape_ptr) {
	auto s = new Entity_Item(shape_ptr);
	change_selected_entity(s);
	// Moche
	if (s->name() == "Camera 0") {
		_camera = s->shape_ptr();
	}
	QListWidgetItem *item = new QListWidgetItem(s->name().c_str());
	item->setData(100, QVariant::fromValue(s));
	_list->addItem(item);
	_list->setCurrentRow(_list->count()-1);
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
