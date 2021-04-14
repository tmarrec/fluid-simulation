#pragma once

#include "Shader.h"
#include "Window.h"
#include "Components.h"
#include "utils.h"
#include <cstdint>

class Renderer
{
public:
	void init(std::shared_ptr<Window> window);
    void prePass() const;
    void initMesh(Mesh& mesh) const;
    void drawMesh(Mesh& mesh) const;
    void freeMesh(Mesh& mesh) const;
    void applyMaterial(Material& material, Camera& camera, Transform& transform) const;
    void initMaterial(Material& material) const;
    void updateTexture(const std::vector<std::uint8_t>& texture, const std::uint32_t textureGL) const;

private:
	std::shared_ptr<Window> _window = nullptr;
};

