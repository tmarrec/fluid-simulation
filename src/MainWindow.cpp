#include "config.h"
#include "MainWindow.h"

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMenu>
#include <QMenuBar>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QFileDialog>
#include <QtGlobal>
#include <QPixmap>

#include <iostream>
#include <sstream>

MainWindow::MainWindow()
	: _selected_entity{nullptr}
	, _last_mouse_x{0}
	, _last_mouse_y{0}
{
	setWindowTitle(TITLE);

	
	// Parametres OpenGL 
	QSurfaceFormat format;
	format.setVersion(4, 1);
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setSwapInterval(0);
	format.setSamples(8);
	QSurfaceFormat::setDefaultFormat(format);
	
	// Widget OpenGL
	_glw = new GLWidget(this);
	_glw->setFocusPolicy(Qt::ClickFocus);

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
	QAction *help = new QAction("Help", this);
	menuBar()->addAction(help);
	connect(help, &QAction::triggered, this, &MainWindow::show_help_box);

	// QMessageBox d'aide
	_help_box = new QMessageBox();
	_help_box->setWindowTitle(TITLE_HELP);
	_help_box->setFixedSize(500, 500);
	_help_box->setIconPixmap(QPixmap("images/works.png"));
	_help_box->setText("--- Tristan Marrec @ Paul Sabatier University @ 2020 ---");
	_help_box->setInformativeText("ZQSD : Camera movements\nShift/Ctrl : Camera UP/DOWN Y-axis\nMouse Click + Movements : Camera Yaw/Pitch\n\nSelect entity on the left panel for parameters");
	
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

	// Camera option
	_camera_box_group = camera_box();
	side_panel_l->addWidget(_camera_box_group);
	// Shaders
	_shaders_box_group = shaders_box();
	_position_box_group = position_box();
	_rotation_box_group = rotation_box();
	_scale_box_group = scale_box();

	side_panel_l->addWidget(_shaders_box_group);
	side_panel_l->addWidget(_position_box_group);
	side_panel_l->addWidget(_rotation_box_group);
	side_panel_l->addWidget(_scale_box_group);

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
	connect(add_model, &QAction::triggered, this, &MainWindow::search_model_file);

	connect(_check_box_camera_face, &QCheckBox::toggled, _openGL, &OpenGL::set_draw_fill);


}

MainWindow::~MainWindow() {

}

void MainWindow::show_help_box() const {
	_help_box->exec();
}

// Change le titre avec des infos (fps actuellement)
void MainWindow::update_title_infos(std::string infos) {
	std::string name = TITLE;
	std::string new_title = name + " | " + infos;
	setWindowTitle(new_title.c_str());
}

// Ouvre une fenetre pour selectionner un .obj et l'ajoute a la scene
void MainWindow::search_model_file() {
	QFileDialog dialog(this);
	dialog.setOption(QFileDialog::DontUseNativeDialog);
	dialog.setFileMode(QFileDialog::ExistingFiles);
	dialog.setViewMode(QFileDialog::Detail);
	dialog.setDirectory("models");
	dialog.setNameFilter("*.obj");
	QStringList file_names;
	if (dialog.exec()) {
		file_names = dialog.selectedFiles();
	}
	if (file_names.size() > 0) {
		_openGL->add_model(file_names.at(0).toUtf8().constData());
	}
}

// Sauvegarde les anciennes positions de la souris
void MainWindow::mousePressEvent(QMouseEvent *event) {
	_last_mouse_x = event->x();
	_last_mouse_y = event->y();
}

// Calcule le yaw et pitch de la camera en fonction
// des mouvements de la souris
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

// Deplace la caméra
void MainWindow::keyPressEvent(QKeyEvent *event) {
	auto key = event->key();
	auto p = _camera->position();
	auto camera = std::static_pointer_cast<Camera>(_camera);
	auto front = camera->front();
	auto up = camera->up();
	float speed = camera->speed();
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
			_camera->set_position(glm::vec3{p.x, p.y+speed, p.z});
			break;
		case Qt::Key_Control:
			_camera->set_position(glm::vec3{p.x, p.y-speed, p.z});
			break;
	}
}

// Enleve une entité du tree view a gauche
void MainWindow::delete_item_entities_tree_view() {
	_openGL->remove_entity(_selected_entity->entity_ptr());
	_list->takeItem(_list->currentRow());

	_list->setCurrentRow(_list->count()-1);
	auto new_item = _list->item(_list->count()-1);
	auto s = new_item->data(100).value<Entity_Item*>();
	change_selected_entity(s);
}	

///////////////// SLIDE TRANSLATION /////////////////
void MainWindow::change_slide_x_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{value, p.y, p.z});
	}
}
void MainWindow::change_slide_y_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, value, p.z});
	}
}
void MainWindow::change_slide_z_position(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto p = s->position();
		s->set_position(glm::vec3{p.x, p.y, value});
	}
}
void MainWindow::update_slide_position(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->entity_ptr()->id() == id) {
		_slide_x_position->setValue(pos.x);
		_slide_x_position_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_position->setValue(pos.y);
		_slide_y_position_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_position->setValue(pos.z);
		_slide_z_position_label->setText(std::to_string(int(pos.z)).c_str());
	}
}
/////////////// FIN SLIDE TRANSLATION ///////////////

////////////////// SLIDE ROTATIONS //////////////////
void MainWindow::change_slide_x_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{value, r.y, r.z});
	}
}
void MainWindow::change_slide_y_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, value, r.z});
	}
}
void MainWindow::change_slide_z_rotation(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->rotation();
		s->set_rotation(glm::vec3{r.x, r.y, value});
	}
}
void MainWindow::update_slide_rotation(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->entity_ptr()->id() == id) {
		_slide_x_rotation->setValue(pos.x);
		_slide_x_rotation_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_rotation->setValue(pos.y);
		_slide_y_rotation_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_rotation->setValue(pos.z);
		_slide_z_rotation_label->setText(std::to_string(int(pos.z)).c_str());
	}
}
//////////////// FIN SLIDE ROTATIONS ////////////////

////////////////// SLIDE SCALE ///////////////////
void MainWindow::change_slide_x_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{value, r.y, r.z});
	}
}
void MainWindow::change_slide_y_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, value, r.z});
	}
}
void MainWindow::change_slide_z_scale(int value) {
	if (_selected_entity != nullptr) {
		auto s = _selected_entity->entity_ptr();
		auto r = s->scale();
		s->set_scale(glm::vec3{r.x, r.y, value});
	}
}
void MainWindow::update_slide_scale(glm::vec3 pos, const unsigned long id) {
	if (_selected_entity != nullptr && _selected_entity->entity_ptr()->id() == id) {
		_slide_x_scale->setValue(pos.x);
		_slide_x_scale_label->setText(std::to_string(int(pos.x)).c_str());
		_slide_y_scale->setValue(pos.y);
		_slide_y_scale_label->setText(std::to_string(int(pos.y)).c_str());
		_slide_z_scale->setValue(pos.z);
		_slide_z_scale_label->setText(std::to_string(int(pos.z)).c_str());
	}
}
//////////////// FIN SLIDE SCALE /////////////////

// Change l'entité séléctionnée au clique sur un tree view
void MainWindow::on_item_clicked(QListWidgetItem *item) {
	auto s = item->data(100).value<Entity_Item*>();
	change_selected_entity(s);
}

// Récupère le dernier string d'un split de string avec un delim c
std::string MainWindow::last_split(std::string s, char c) const {
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);
	while (std::getline(tokenStream, token, c))
	{
		tokens.push_back(token);
	}
	return tokens.back();
}

// Au changement d'entité selectionné :
// Mets a jours les sliders
// Active/Desactives les buttons de propriétées
void MainWindow::change_selected_entity(Entity_Item* e) {
	_selected_entity = e;
	auto s = e->entity_ptr();
	update_slide_position(s->position(), s->id());
	update_slide_rotation(s->rotation(), s->id());
	update_slide_scale(s->scale(), s->id());

	if (s->type() == CAMERA) {
		_delete_button->setEnabled(false);
		_camera_box_group->setVisible(true);
	} else {
		_delete_button->setEnabled(true);
		_camera_box_group->setVisible(false);
	}

	if (s->type() == SHAPE) {
		auto vert = last_split(s->shader().vert_path(), '/');
		auto i = _combo_box_shaders_vert->findText(QString::fromStdString(vert));
		if (i != -1) {
			_combo_box_shaders_vert->setCurrentIndex(i);
		}
		auto frag = last_split(s->shader().frag_path(), '/');
		auto j = _combo_box_shaders_frag->findText(QString::fromStdString(frag));
		if (i != -1) {
			_combo_box_shaders_frag->setCurrentIndex(j);
		}
	}

	if (s->type() == CAMERA || s->type() == LIGHT) {
		_combo_box_shaders_frag->setEnabled(false);
		_combo_box_shaders_vert->setEnabled(false);
		_rotation_box_group->setEnabled(false);
		_scale_box_group->setEnabled(false);
	} else {
		_combo_box_shaders_frag->setEnabled(true);
		_combo_box_shaders_vert->setEnabled(true);
		_rotation_box_group->setEnabled(true);
		_scale_box_group->setEnabled(true);
	}
}

// Ajoute un item au tree view avec une Entity_Item qui stocke
// le pointeur de l'entitée
Q_DECLARE_METATYPE(Entity_Item*)
void MainWindow::add_item_to_QListW(std::shared_ptr<Entity> entity_ptr) {
	auto s = new Entity_Item(entity_ptr);
	change_selected_entity(s);
	if (s->entity_ptr()->type() == CAMERA) {
		_camera = s->entity_ptr();
	}
	QListWidgetItem *item = new QListWidgetItem(s->name().c_str());
	item->setData(100, QVariant::fromValue(s));
	_list->addItem(item);
	_list->setCurrentRow(_list->count()-1);
}

// Change le vert shader de l'entitée séléctionnée
void MainWindow::change_vert_shader(int i) {
	auto s = _selected_entity->entity_ptr();
	std::string file = _combo_box_shaders_vert->itemText(i).toUtf8().constData();
	std::string dir = "shaders/";
	std::string file_path = dir+file;

	auto frag_path = s->shader().frag_path();
	Shader shader {file_path.c_str(), frag_path.c_str()};
	s->set_shader(shader);
}

// Change le frag shader de l'entitée séléctionnée
void MainWindow::change_frag_shader(int i) {
	auto s = _selected_entity->entity_ptr();
	std::string file = _combo_box_shaders_frag->itemText(i).toUtf8().constData();
	std::string dir = "shaders/";
	std::string file_path = dir+file;

	auto vert_path = s->shader().vert_path();
	Shader shader {vert_path.c_str(), file_path.c_str()};
	s->set_shader(shader);
}

// Change la vitesse de la camera
void MainWindow::change_camera_speed(const QString & speed) {
	auto camera = std::static_pointer_cast<Camera>(_camera);
	camera->set_speed(speed.toFloat());
}

// Creer la box de propriétées de la camera
QGroupBox* MainWindow::camera_box() {
	QGroupBox *box = new QGroupBox("Camera options");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *camera_g_box = new QGroupBox();
	QVBoxLayout *camera_g_layout = new QVBoxLayout;

	_check_box_camera_face = new QCheckBox("Render faces");
	_check_box_camera_face->setChecked(true);
	camera_g_layout->addWidget(_check_box_camera_face);

	QGroupBox *camera_speed_box = new QGroupBox();
	QHBoxLayout *camera_speed_layout = new QHBoxLayout;
	auto title_label = new QLabel("Speed");
	camera_speed_layout->addWidget(title_label);
	_input_dialog_camera_speed = new QLineEdit("5");
	_input_dialog_camera_speed->setValidator(new QDoubleValidator(0, 100, 2, this));
	camera_speed_layout->addWidget(_input_dialog_camera_speed);
	camera_speed_box->setLayout(camera_speed_layout);


	camera_g_layout->addWidget(camera_speed_box);

	camera_g_box->setLayout(camera_g_layout);
	box_layout->addWidget(camera_g_box);
	box->setLayout(box_layout);

	connect(_input_dialog_camera_speed, QOverload<const QString &>::of(&QLineEdit::textChanged), this, &MainWindow::change_camera_speed);
	return box;
}


// Creer la box de propriétées des shaders
QGroupBox* MainWindow::shaders_box() {
	QGroupBox *box = new QGroupBox("Shaders");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *shaders_box = new QGroupBox();
	QHBoxLayout *shaders_layout = new QHBoxLayout;

	box->setMaximumHeight(80);
	QDir dir{"shaders"};

	_combo_box_shaders_vert = new QComboBox();
	QStringList verts = dir.entryList(QStringList() << "*.vert", QDir::Files);
	foreach(QString filename, verts) {
		_combo_box_shaders_vert->addItem(filename.toUtf8().constData());
	}
	shaders_layout->addWidget(_combo_box_shaders_vert);

	_combo_box_shaders_frag = new QComboBox();
	QStringList frags = dir.entryList(QStringList() << "*.frag", QDir::Files);
	foreach(QString filename, frags) {
		_combo_box_shaders_frag->addItem(filename.toUtf8().constData());
	}
	shaders_layout->addWidget(_combo_box_shaders_frag);
	shaders_box->setLayout(shaders_layout);

	box_layout->addWidget(shaders_box);
	
	connect(_combo_box_shaders_vert, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::change_vert_shader);
	connect(_combo_box_shaders_frag, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &MainWindow::change_frag_shader);
	
	box->setLayout(box_layout);
	return box;
}


// Creer la box de propriétées de la position
QGroupBox* MainWindow::position_box() {
	QGroupBox *box = new QGroupBox("Position");
	QVBoxLayout *box_layout = new QVBoxLayout;	

	QGroupBox *slide_x_box = new QGroupBox();
	QHBoxLayout *slide_x_layout = new QHBoxLayout;
	_slide_x_position_label = new QLabel("x");
	_slide_x_position = new QSlider(Qt::Orientation::Horizontal);
	_slide_x_position->setMinimumWidth(180);
	_slide_x_position->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
	_slide_x_position->setMinimum(-2000);
	_slide_x_position->setMaximum(2000);
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
	_slide_y_position->setMinimum(-2000);
	_slide_y_position->setMaximum(2000);
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
	_slide_z_position->setMinimum(-2000);
	_slide_z_position->setMaximum(2000);
	slide_z_layout->addWidget(_slide_z_position);
	slide_z_layout->addWidget(_slide_z_position_label);
	slide_z_box->setLayout(slide_z_layout);
	box_layout->addWidget(slide_z_box);
	connect(_slide_z_position, &QSlider::sliderMoved, this, &MainWindow::change_slide_z_position);
	
	box->setLayout(box_layout);
	return box;
}

// Creer la box de propriétées de la rotation
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

// Creer la box de propriétées de la scale
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
	_slide_x_scale->setMaximum(100);
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
	_slide_y_scale->setMaximum(100);
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
	_slide_z_scale->setMaximum(100);
	slide_z_layout->addWidget(_slide_z_scale);
	slide_z_layout->addWidget(_slide_z_scale_label);
	slide_z_box->setLayout(slide_z_layout);
	box_layout->addWidget(slide_z_box);
	connect(_slide_z_scale, &QSlider::sliderMoved, this, &MainWindow::change_slide_z_scale);
	
	box->setLayout(box_layout);
	return box;
}
