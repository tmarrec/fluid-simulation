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

private:
	float _FOV;
	glm::vec3 _front;
	glm::vec3 _up;
	float _yaw;
	float _pitch;
	float _speed;
};
