#pragma once

#include <GL/gl.h>
#include <GL/glext.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"

class Shader {

public:

	Shader(const char* __vertPath, const char* __fragPath)
		: _vertPath { __vertPath }
		, _fragPath { __fragPath }
	{
		_init();	
	}

	Shader(const char* __vertPath, const char* __fragPath, const char* __geoPath)
		: _vertPath { __vertPath }
		, _fragPath { __fragPath }
		, _geoPath { __geoPath }
	{
		_init();	
	}
	
	void set_1i(const std::string &name, int value) const { glUniform1i(_getLocation(name), value); }
	void set_1f(const std::string &name, float value) const { glUniform1f(_getLocation(name), value); }
	void set_4f(const std::string &name, glm::vec4 values) const { glUniform4f(_getLocation(name), values.x, values.y, values.z, values.w); }
	void set_3f(const std::string &name, glm::vec3 values) const { glUniform3f(_getLocation(name), values.x, values.y, values.z); }
	void set_mat4(const std::string &name, glm::mat4 values) const { glUniformMatrix4fv(_getLocation(name), 1, GL_FALSE, glm::value_ptr(values)); }

	void use() const
	{
		glUseProgram(_id);
	}

	~Shader() {}

private:
	std::uint64_t _id;
	const std::string _vertPath;
	const std::string _fragPath;
	const std::string _geoPath;

	static void _checkCompilation(std::uint64_t shader, const std::string& type)
	{
		GLint success;
		GLsizei infoLogLength = 0;
		if (type != "PROGRAM")
		{
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			char* infoLog = new char[infoLogLength];
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			if (!success)
			{
				glGetShaderInfoLog(shader, infoLogLength, NULL, infoLog);
				ERROR("Shader compilation error of type : " << type << "\n" << infoLog << "\n");
			}
			delete[] infoLog;
		}
		else
		{
			glGetProgramiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
			char* infoLog = new char[infoLogLength];
			glGetProgramiv(shader, GL_LINK_STATUS, &success);
			if (!success)
			{
				glGetProgramInfoLog(shader, infoLogLength, NULL, infoLog);
				ERROR("Program linking error of type : " << type << "\n" << infoLog << "\n");
			}
			delete[] infoLog;
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
		std::string geoCode;
		std::ifstream vShaderFile;
		std::ifstream fShaderFile;
		std::ifstream gShaderFile;
		// ensure ifstream objects can throw exceptions:
		vShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
		fShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
		gShaderFile.exceptions (std::ifstream::failbit | std::ifstream::badbit);
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
			vertexCode = vShaderStream.str();
			fragmentCode = fShaderStream.str();
			if (!_geoPath.empty())
			{
				gShaderFile.open(_geoPath);
				std::stringstream gShaderStream;
				gShaderStream << gShaderFile.rdbuf();
				gShaderFile.close();
				geoCode = gShaderStream.str();
			}
		}
		catch (const std::ifstream::failure & e)
		{
			ERROR("Shader file not found");
		}

		const char* vShaderCode = vertexCode.c_str();
		const char* fShaderCode = fragmentCode.c_str();
		// 2. compile shaders
		std::uint64_t vertex, fragment;
		// Vertex shader
		vertex = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vertex, 1, &vShaderCode, NULL);
		glCompileShader(vertex);
		_checkCompilation(vertex, "VERTEX");
		// Fragment Shader
		fragment = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fragment, 1, &fShaderCode, NULL);
		glCompileShader(fragment);
		_checkCompilation(fragment, "FRAGMENT");
		// Geometry Shader
		std::uint64_t geometry = -1;
		if (!_geoPath.empty())
		{
			const char* gShaderCode = geoCode.c_str();
			geometry = glCreateShader(GL_GEOMETRY_SHADER);
			glShaderSource(geometry, 1, &gShaderCode, NULL);
			glCompileShader(geometry);
			_checkCompilation(geometry, "GEOMETRY");
		}
		// Shader Program linking
		_id = glCreateProgram();
		glAttachShader(_id, vertex);
		glAttachShader(_id, fragment);
		if (!_geoPath.empty())
		{
			glAttachShader(_id, geometry);
		}
		glLinkProgram(_id);
		_checkCompilation(_id, "PROGRAM");
		// Delete shaders as they're linked into our program
		glDeleteShader(vertex);
		glDeleteShader(fragment);
		if (!_geoPath.empty())
		{
			glDeleteShader(geometry);
		}
	}
};
