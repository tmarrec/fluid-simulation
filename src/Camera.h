#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "Entity.h"

class Camera : public Entity {

public:
	Camera(glm::vec3 position, glm::vec3 rotation, glm::vec3 scale,
		float FOV, MainWindow * main_window);
	~Camera(void);
	
	void draw(glm::mat4 view, glm::mat4 projection, float delta_time,
				std::vector<std::shared_ptr<Entity>> lights) override;
	
	float FOV() const;
	float yaw() const;
	float pitch() const;
	glm::vec3 front() const;
	glm::vec3 up() const;
	glm::mat4 view() const;
	void set_front(glm::vec3 front);
	void set_yaw(float yaw);
	void set_pitch(float pitch);
	void set_shader(Shader shader) override;
	Shader & shader() override;
	Entity_Type type() override;
	float speed() const;
	void set_speed(float speed);

	void set_move_front(bool state);
	void set_move_back(bool state);
	void set_move_left(bool state);
	void set_move_right(bool state);
	void set_move_up(bool state);
	void set_move_down(bool state);

	bool move_front() const;
	bool move_back() const;
	bool move_left() const;
	bool move_right() const;
	bool move_up() const;
	bool move_down() const;

private:
	float _FOV;
	glm::vec3 _front;
	glm::vec3 _up;
	float _yaw;
	float _pitch;
	float _speed;

	bool _move_front;
	bool _move_back;
	bool _move_left;
	bool _move_right;
	bool _move_up;
	bool _move_down;
};
