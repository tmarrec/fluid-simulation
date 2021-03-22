#pragma once

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"
#include <glad/glad.h>

class Shader
{

public:
	void use() const;

    void setVert(const std::string &filename);
    void setFrag(const std::string &filename);
    void setGeo(const std::string &filename);

	void set1i(const std::string &name, const int value) const;
	void set1f(const std::string &name, const float value) const;
	void set4f(const std::string &name, const glm::vec4 values) const;
	void set3f(const std::string &name, const glm::vec3 values) const;
	void setMat4(const std::string &name, const glm::mat4 values) const;


private:
	GLuint _id = 0;
	std::string _vertPath;
	std::string _fragPath;
	std::string _geoPath;

	void init();
	static void checkCompilation(std::uint64_t shader, const std::string& type);
	int getLocation(const std::string &name) const;
};

