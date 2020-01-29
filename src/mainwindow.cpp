#include "mainwindow.h"
#include "./ui_mainwindow.h"

#include <iostream>

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
	, ui(new Ui::MainWindow)
{
	ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete ui;
}


void MainWindow::on_triangleButton_clicked()
{
	std::cout << "triangle" << std::endl;
}

void MainWindow::on_rectangleButton_clicked()
{
	std::cout << "rect" << std::endl;
}
