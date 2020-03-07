#pragma once

class GLWidget;

#include <QMainWindow>
#include <QListWidget>
#include <QLabel>
#include <QGroupBox>

#include "GLWidget.h"
#include "Entity_Item.h"

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();
	void add_item_to_QListW(std::shared_ptr<Entity> shape_ptr);
	void change_selected_entity(Entity_Item* e);

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


private:
	QListWidget* _list;
	GLWidget* _glw;
	OpenGL* _openGL;

	QLabel *_slide_x_position_label;
	QLabel *_slide_y_position_label;
	QLabel *_slide_z_position_label;

	QLabel *_slide_x_rotation_label;
	QLabel *_slide_y_rotation_label;
	QLabel *_slide_z_rotation_label;
	
	QLabel *_slide_x_scale_label;
	QLabel *_slide_y_scale_label;
	QLabel *_slide_z_scale_label;

	Entity_Item* _selected_entity;

	QGroupBox* position_box();
	QGroupBox* rotation_box();
	QGroupBox* scale_box();
};
