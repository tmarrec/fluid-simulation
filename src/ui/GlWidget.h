#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>

#include "../MessageBus.h"
#include "../System.h"

class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core, public System
{
Q_OBJECT
public:
	explicit GlWidget(QWidget *parent = nullptr);

	~GlWidget() final;

	void cout(std::string string) const;

public slots:
	void cleanup();

protected:
	void initializeGL() override;
	void paintGL() override;

	void keyPressEvent(QKeyEvent *event) override;

private:

};

