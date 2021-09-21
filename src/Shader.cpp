#include "Shader.h"

void Shader::use() const
{
    glUseProgram(_id);
}

void Shader::setVert(const std::string &filename)
{
    _vertPath = filename;
    init();
}

void Shader::setFrag(const std::string &filename)
{
    _fragPath = filename;
    init();
}

void Shader::set1f(const std::string &name, const float value) const
{
    glUniform1f(getLocation(name), value);
}

void Shader::set3f(const std::string &name, const glm::vec3 values) const
{
    glUniform3f(getLocation(name), values.x, values.y, values.z);
}

void Shader::set2f(const std::string &name, const glm::vec2 values) const
{
    glUniform2f(getLocation(name), values.x, values.y);
}

void Shader::setMat4(const std::string &name, const glm::mat4 values) const
{
    glUniformMatrix4fv(getLocation(name), 1, GL_FALSE, glm::value_ptr(values));
}

void Shader::checkCompilation(std::uint64_t shader, const std::string& type)
{
    GLint success;
    GLsizei infoLogLength = 0;
    if (type != "PROGRAM")
    {
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success)
        {
            char* infoLog = new char[infoLogLength];
            glGetShaderInfoLog(shader, infoLogLength, NULL, infoLog);
            ERROR("Shader compilation error of type : " <<
                    type << "\n" << infoLog << "\n");
        }
    }
    else
    {
        glGetProgramiv(shader, GL_INFO_LOG_LENGTH, &infoLogLength);
        glGetProgramiv(shader, GL_LINK_STATUS, &success);
        if (!success)
        {
            char* infoLog = new char[infoLogLength];
            glGetProgramInfoLog(shader, infoLogLength, NULL, infoLog);
            ERROR("Program linking error of type : " <<
                    type << "\n" << infoLog << "\n");
        }
    }
}

GLint Shader::getLocation(const std::string &name) const
{
    const GLint location = glGetUniformLocation(_id, name.c_str());
    if (location == -1) {
        WARNING("Cannot find uniform location : " << name);
    }
    return location;
}

void Shader::init()
{
    // Could create trouble..
    glDeleteProgram(_id);

    // Retrieve the vertex/fragment source code
    std::string vertexCode;
    std::string fragmentCode;
    std::ifstream vShaderFile;
    std::ifstream fShaderFile;
    std::ifstream gShaderFile;

    // Check if files exists
    vShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    fShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    gShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    try
    {
        if (!_vertPath.empty())
        {
            vShaderFile.open(_vertPath);
            std::stringstream vShaderStream;
            vShaderStream << vShaderFile.rdbuf();
            vShaderFile.close();
            vertexCode = vShaderStream.str();
        }
        if (!_fragPath.empty())
        {
            fShaderFile.open(_fragPath);
            std::stringstream fShaderStream;
            fShaderStream << fShaderFile.rdbuf();
            fShaderFile.close();
            fragmentCode = fShaderStream.str();
        }
    }
    catch (const std::ifstream::failure & e)
    {
        ERROR("Shader file not found");
    }

    // Compile shaders
    // Vertex shader
    GLuint vertex = 0;
    if (!_vertPath.empty())
    {
        const char* vShaderCode = vertexCode.c_str();
        vertex = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex, 1, &vShaderCode, NULL);
        glCompileShader(vertex);
        checkCompilation(vertex, "VERTEX");
    }
    // Fragment Shader
    GLuint fragment = 0;
    if (!_fragPath.empty())
    {
        const char* fShaderCode = fragmentCode.c_str();
        fragment = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment, 1, &fShaderCode, NULL);
        glCompileShader(fragment);
        checkCompilation(fragment, "FRAGMENT");
    }

    // Shader Program linking
    _id = glCreateProgram();
    if (!_vertPath.empty())
    {
        glAttachShader(_id, vertex);
    }
    if (!_fragPath.empty())
    {
        glAttachShader(_id, fragment);
    }

    glLinkProgram(_id);
    checkCompilation(_id, "PROGRAM");

    // Delete shaders as they are linked into our program
    if (!_vertPath.empty())
    {
        glDeleteShader(vertex);
    }
    if (!_fragPath.empty())
    {
        glDeleteShader(fragment);
    }
}
