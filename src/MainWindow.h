#pragma once

#include <QMainWindow>

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow();
	~MainWindow();

private slots:
	void on_triangleButton_clicked();
	void on_rectangleButton_clicked();

private:
	void add_item_to_QListW(std::string name);

};
