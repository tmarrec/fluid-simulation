#pragma once

class GLWidget;

#include <QMainWindow>
#include <QListWidget>

#include "GLWidget.h"

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();
	void add_item_to_QListW(uint id, std::string name);

private slots:
	void on_item_clicked(QListWidgetItem *item);
	void on_triangleButton_clicked();
	void on_rectangleButton_clicked();
	void add_triangle();
	void add_cube();
	void add_sphere();
	void change_slide_x(int value);
	void change_slide_y(int value);
	void change_slide_z(int value);

private:
	QListWidget *_list;
	GLWidget *_glw;
	OpenGL *_openGL;

	uint _selected_entity;
};
