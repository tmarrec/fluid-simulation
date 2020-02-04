#pragma once

#include <string>

class Shader {

public:
	Shader(const char* vert_path, const char* frag_path);
	~Shader();
	void use();
	void set_bool(const std::string &name, bool value) const;	
	void set_int(const std::string &name, int value) const;	
	void set_float(const std::string &name, float value) const;	

private:
	unsigned int _id;

	void check_compilation(unsigned int shader, std::string type);
};
