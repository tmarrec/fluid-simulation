#pragma once

#include <GL/gl.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"

class Shader {

public:

	Shader(const char* vertPath, const char* fragPath)
		: _vertPath{vertPath}
		, _fragPath{fragPath}
	{
		_init();	
	}
	
	void set_1i(const std::string &name, int value) const { glUniform1i(_getLocation(name), value); }
	void set_1f(const std::string &name, float value) const { glUniform1f(_getLocation(name), value); }
	void set_4f(const std::string &name, glm::vec4 values) const { glUniform4f(_getLocation(name), values.x, values.y, values.z, values.w); }
	void set_3f(const std::string &name, glm::vec3 values) const { glUniform3f(_getLocation(name), values.x, values.y, values.z); }
	void set_mat4(const std::string &name, glm::mat4 values) const { glUniformMatrix4fv(_getLocation(name), 1, GL_FALSE, glm::value_ptr(values)); }
	std::string vertPath() const { return _vertPath; }
	std::string fragPath() const { return _fragPath; }

	void use() const
	{
		glUseProgram(_id);
	}

	~Shader() {}

private:
	std::uint64_t _id;
	std::string _vertPath;
	std::string _fragPath;

	void _checkCompilation(std::uint64_t shader, std::string type)
	{
		GLint success;
		GLsizei infoLogLength = 0;
		if (type != "PROGRAM") {
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			char* infoLog = new char[infoLogLength];
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			if (!success) {
				glGetShaderInfoLog(shader, infoLogLength, NULL, infoLog);
				ERROR("Shader compilation error of type : " << type << "\n" << infoLog << "\n");
			}
		} else {
			glGetProgramiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			char* infoLog = new char[infoLogLength];
			glGetProgramiv(shader, GL_LINK_STATUS, &success);
			if (!success) {
				glGetProgramInfoLog(shader, infoLogLength, NULL, infoLog);
				ERROR("Program linking error of type : " << type << "\n" << infoLog << "\n");
			}
		}
	}

	int _getLocation(const std::string &name) const
	{
		int location = glGetUniformLocation(_id, name.c_str());
		if (location == -1) {
			WARNING("Cannot find uniform location : " << name);
		}	
		return location;
	}
	
	void _init()
	{
		// 1. retrieve the vertex/fragment source code from filePath
		std::string vertexCode;
		std::string fragmentCode;
		std::ifstream vShaderFile;
		std::ifstream fShaderFile;
		// ensure ifstream objects can throw exceptions:
		vShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
		fShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
		try
		{
			// open files
			vShaderFile.open(_vertPath);
			fShaderFile.open(_fragPath);
			std::stringstream vShaderStream, fShaderStream;
			// read file's buffer contents into streams
			vShaderStream << vShaderFile.rdbuf();
			fShaderStream << fShaderFile.rdbuf();
			// close file handlers
			vShaderFile.close();
			fShaderFile.close();
			// convert stream into string
			vertexCode   = vShaderStream.str();
			fragmentCode = fShaderStream.str();
		}
		catch (const std::ifstream::failure & e)
		{
			std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ" << std::endl;
		}

		const char* vShaderCode = vertexCode.c_str();
		const char* fShaderCode = fragmentCode.c_str();
		// 2. compile shaders
		std::uint64_t vertex, fragment;
		// vertex shader
		vertex = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertex, 1, &vShaderCode, NULL);
		glCompileShader(vertex);
		_checkCompilation(vertex, "VERTEX");
		// fragment Shader
		fragment = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragment, 1, &fShaderCode, NULL);
		glCompileShader(fragment);
		_checkCompilation(fragment, "FRAGMENT");
		// shader Program
		_id = glCreateProgram();
		glAttachShader(_id, vertex);
		glAttachShader(_id, fragment);
		glLinkProgram(_id);
		_checkCompilation(_id, "PROGRAM");
		// delete the shaders as they're linked into our program now and no longer necessary
		glDeleteShader(vertex);
		glDeleteShader(fragment);
	}
};
