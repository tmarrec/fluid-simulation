#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "./glad/glad.h"
#include "./glm/gtc/type_ptr.hpp"
#include "./utils.h"

class Shader
{
 public:
    void use() const;

    void setVert(const std::string &filename);
    void setFrag(const std::string &filename);

    void set1f(const std::string &name, const float value) const;
    void set3f(const std::string &name, const glm::vec3 values) const;
    void set2f(const std::string &name, const glm::vec2 values) const;
    void setMat4(const std::string &name, const glm::mat4 values) const;


 private:
    GLuint _id = 0;
    std::string _vertPath;
    std::string _fragPath;

    void init();
    static void checkCompilation(std::uint64_t shader, const std::string& type);
    int getLocation(const std::string &name) const;
};

