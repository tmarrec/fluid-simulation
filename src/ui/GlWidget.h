#pragma once

#include <QOpenGLWidget>
#include <QOpenGLFunctions_4_5_Core>
#include <QKeyEvent>
#include <chrono>
#include <cstdint>
#include <thread>

#include "../utils.h"

#include "../ECS.h"
#include "../TransformComponent.h"
#include "../CameraComponent.h"
#include "../DrawableComponent.h"

class Renderer;
class ECS_Manager;
class InputManager;
class Entity;
class Test;

using Renderer__ = std::shared_ptr<Renderer>; 
using ECS_Manager__ = std::shared_ptr<ECS_Manager>; 
using InputManager__ = std::shared_ptr<InputManager>; 

class GlWidget : public QOpenGLWidget, protected QOpenGLFunctions_4_5_Core
{
Q_OBJECT
public:
	explicit GlWidget(QWidget *parent = nullptr);
	Entity* getActiveCamera() const;

	~GlWidget();

public slots:
	void cleanup();

protected:
	void initializeGL() override;
	void paintGL() override;
	void resizeGL(int w, int h) override;

	void keyPressEvent(QKeyEvent *event) override;
	void keyReleaseEvent(QKeyEvent *event) override;

private:
	void _init();
	const Renderer__ _renderer;
	const ECS_Manager__ _ECS_manager;
	const InputManager__ _InputManager;
	Entity* _activeCamera; //TODO Change that

	std::uint64_t _frame_count = 0;
	std::uint64_t _start_timer_fps = 0;
	std::uint64_t _start_timer_frame = 0;
	double _deltaTime;

};

class Test
{
public:
	Test(ECS_Manager__ __ECS_manager)
	{
		_currentTime = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		_ECS_manager = __ECS_manager;
	}

	void operator()()
	{
		for (;;)
		{
			/*
			double newTime = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
			double frameTime = newTime - _currentTime;
			if (frameTime > (double)_tick_rate/100)
			{
				frameTime = (double)_tick_rate/100;
			}
			_currentTime = newTime;
			_accumulator += frameTime;

			while (_accumulator >= _deltaTime)
			{
				_time += _deltaTime;
				_accumulator -= _deltaTime;
				//ECS_manager->update(); //TODO a faire dans la gameloop
				std::cout << _time << " " << _deltaTime << std::endl;
				//std::this_thread::sleep_for(std::chrono::milliseconds(400));
				ECS_manager->update(_deltaTime);
			}
			*/

		}
	}

	double _time = 0.0;
	double _deltaTime = 0.01f;
	double _currentTime = 0.0f;
	double _accumulator = 0.0f;
	const std::uint8_t _tick_rate = 64;
	ECS_Manager__ _ECS_manager;
	InputManager__ _InputManager;

};
