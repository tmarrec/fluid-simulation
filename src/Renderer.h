#pragma once

#include "Shader.h"
#include "Components.h"
#include "types.h"
#include "utils.h"
#include "Window.h"

#include <cstdint>
#include <memory>
#include <chrono>
#include <deque>

struct FrameBuffer
{
    GLuint FBO;
    GLuint RBO;
    GLuint VAO;
    GLuint VBO;
    GLuint texture;
    std::shared_ptr<Shader> shader;
};

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
    void prePass();
    void endPass() const;
    void initMesh(Mesh& mesh) const;
    void drawMesh(Mesh& mesh) const;
    void freeMesh(Mesh& mesh) const;
    void applyMaterial(Material& material, Camera& camera, Transform& transform) const;
    void initMaterial(Material& material) const;
    void initTexture3D(const std::vector<std::uint8_t>& texture, const std::uint32_t textureGL) const;
    void writeImg(const std::uint32_t iteration) const;

private:
	std::shared_ptr<Window> _window = nullptr;
    FrameBuffer _screenbuffer {};
    FrameBuffer _raymarchingbuffer {};
    WindowInfos _windowInfos {};
    std::uint32_t _iterations = 0;

    void _initFrameBuffer(FrameBuffer& framebuffer, std::string vert, std::string frag);
};

