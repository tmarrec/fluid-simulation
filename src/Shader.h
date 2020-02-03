#pragma once

class Shader {

public:
	Shader(const char* vert_path, const char* frag_path);
	~Shader();
	void use();

private:
	unsigned int _id;

};
