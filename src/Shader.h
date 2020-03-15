#pragma once

#include <glm/glm.hpp>

#include <string>

class Shader {

public:
	Shader(const char* vert_path, const char* frag_path);
	~Shader();
	void use() const;
	void set_1i(const std::string &name, int value) const;	
	void set_1f(const std::string &name, float value) const;	
	void set_3f(const std::string &name, glm::vec3 values) const;
	void set_4f(const std::string &name, glm::vec4 values) const;
	void set_mat4(const std::string &name, glm::mat4 values) const;

	std::string vert_path() const;
	std::string frag_path() const;

private:
	unsigned int _id;
	std::string _vert_path;
	std::string _frag_path;

	void check_compilation(unsigned int shader, std::string type);
	int get_location(const std::string &name) const;
};
