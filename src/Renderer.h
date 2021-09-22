#pragma once

#include <vector>
#include <memory>
#include <string>
#include <filesystem>

#include "./glm/glm.hpp"
#include "./glm/gtc/matrix_transform.hpp"
#include "./glm/gtc/type_ptr.hpp"
#include "./glm/ext/matrix_clip_space.hpp"
#include "./glm/gtx/string_cast.hpp"

#include "./types.h"
#include "./utils.h"
#include "./Shader.h"
#include "./Window.h"
#include "./config.h"

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
    void init();
    void prePass();
    void endPass() const;
    void raymarchPass() const;
    void initMesh(Mesh& mesh) const;
    void drawMesh(const Mesh& mesh) const;
    void freeMesh(Mesh& mesh) const;
    void applyMaterial(
            const Material& material,
            const Camera& camera,
            const Transform& transform
        ) const;
    void initMaterial(Material& material) const;
    void initTexture3D(
            const std::vector<std::uint8_t>& texture,
            const std::uint32_t textureGL
        ) const;
    void initTexture2D(
            const std::vector<std::uint8_t>& texture,
            const std::uint32_t textureGL
        ) const;
    void writeImg(const std::uint32_t iteration) const;
    void setLineWidth(const float width) const;

 private:
    std::shared_ptr<Window> _window = nullptr;
    FrameBuffer _screenbuffer {};
    FrameBuffer _raymarchingbuffer {};

    void initFrameBuffer(
            FrameBuffer& framebuffer,
            std::string vert,
            std::string frag
        );
};

