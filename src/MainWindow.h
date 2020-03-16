#pragma once

#define RND() static_cast<float>(rand())/static_cast<float>(RAND_MAX)

class GLWidget;
class Entity;

#include <QMainWindow>
#include <QListWidget>
#include <QLabel>
#include <QGroupBox>
#include <QSlider>
#include <QPushButton>
#include <QComboBox>
#include <QCheckBox>
#include <QLineEdit>
#include <QMessageBox>

#include "GLWidget.h"
#include "Entity_Item.h"
#include "Shader.h"

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();
	void add_item_to_QListW(std::shared_ptr<Entity> shape_ptr);
	void change_selected_entity(Entity_Item* e);
	void update_slide_position(glm::vec3 pos, const unsigned long id);
	void update_slide_rotation(glm::vec3 pos, const unsigned long id);
	void update_slide_scale(glm::vec3 pos, const unsigned long id);
	void update_title_infos(std::string infos);

private slots:
	void on_item_clicked(QListWidgetItem *item);
	void delete_item_entities_tree_view();

	void change_slide_x_position(int value);
	void change_slide_y_position(int value);
	void change_slide_z_position(int value);

	void change_slide_x_rotation(int value);
	void change_slide_y_rotation(int value);
	void change_slide_z_rotation(int value);

	void change_slide_x_scale(int value);
	void change_slide_y_scale(int value);
	void change_slide_z_scale(int value);

	void search_model_file();

	void keyPressEvent(QKeyEvent *event) override;
	void mouseMoveEvent(QMouseEvent *event) override;
	void mousePressEvent(QMouseEvent *event) override;

	void change_vert_shader(int i);
	void change_frag_shader(int i);
	void change_camera_speed(const QString &);

	void show_help_box() const;
	void three_points_scene();


private:
	QListWidget* _list;
	GLWidget* _glw;
	OpenGL* _openGL;

	Entity_Item* _selected_entity;
	int _last_mouse_x;
	int _last_mouse_y;

	QSlider *_slide_x_position;
	QSlider *_slide_y_position;
	QSlider *_slide_z_position;
	QSlider *_slide_x_rotation;
	QSlider *_slide_y_rotation;
	QSlider *_slide_z_rotation;
	QSlider *_slide_x_scale;
	QSlider *_slide_y_scale;
	QSlider *_slide_z_scale;

	QLabel *_slide_x_position_label;
	QLabel *_slide_y_position_label;
	QLabel *_slide_z_position_label;

	QLabel *_slide_x_rotation_label;
	QLabel *_slide_y_rotation_label;
	QLabel *_slide_z_rotation_label;
	
	QLabel *_slide_x_scale_label;
	QLabel *_slide_y_scale_label;
	QLabel *_slide_z_scale_label;

	std::shared_ptr<Entity> _camera;

	QGroupBox* position_box();
	QGroupBox* rotation_box();
	QGroupBox* scale_box();
	QGroupBox* shaders_box();
	QGroupBox* camera_box();
	QGroupBox* _position_box_group;
	QGroupBox* _rotation_box_group;
	QGroupBox* _scale_box_group;
	QGroupBox* _shaders_box_group;
	QGroupBox* _camera_box_group;

	QPushButton *_delete_button;

	QComboBox *_combo_box_shaders_vert;
	QComboBox *_combo_box_shaders_frag;

	QCheckBox *_check_box_camera_face;
	QLineEdit *_input_dialog_camera_speed;

	std::string last_split(std::string s, char c) const;
	QMessageBox *_help_box;

};
