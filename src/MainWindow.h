#pragma once

class GLWidget;

#include <QMainWindow>
#include <QListWidget>

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
	void change_slide_x(int value);
	void change_slide_y(int value);
	void change_slide_z(int value);

private:
	QListWidget* _list;
	GLWidget* _glw;
	OpenGL* _openGL;

	Entity_Item* _selected_entity;
};
