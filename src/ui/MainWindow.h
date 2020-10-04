#pragma once

#include <QtWidgets/QMainWindow>

#include "../System.h"

namespace Ui {
	class MainWindow;
}

class MainWindow : public QMainWindow, public System
{
Q_OBJECT

public:
	explicit MainWindow(MsgBus_ptr messageBus);
	void cout(std::string string) const override;
	void handleMessage(Message & msg) override;
	MsgBus_ptr messageBus() const;

	~MainWindow() final;

private:
	Ui::MainWindow * ui;

private slots:

};
