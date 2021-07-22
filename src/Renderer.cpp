#include "Renderer.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

void Renderer::init(const WindowInfos& windowInfos)
{
    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)))
    {
        ERROR("Failed to initialize glad");
    }

    _windowInfos = windowInfos;
    glViewport(0, 0, _windowInfos.x, _windowInfos.y);
    glEnable(GL_DEPTH_TEST);

    _initFrameBuffer(_screenbuffer, "shaders/screen.vert", "shaders/screen.frag");
    _initFrameBuffer(_raymarchingbuffer, "shaders/screen.vert", "shaders/raymarch.frag");
}


void Renderer::prePass()
{
    glBindFramebuffer(GL_FRAMEBUFFER, _screenbuffer.FBO);
    glEnable(GL_DEPTH_TEST);
    glClearColor(0.15f, 0.15f, 0.15f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void Renderer::endPass() const
{
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
    //glClear(GL_COLOR_BUFFER_BIT);
	glDisable(GL_DEPTH_TEST);

    // Screenbuffer step
	_screenbuffer.shader->use();
	glBindVertexArray(_screenbuffer.VAO);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, _screenbuffer.texture);
    glDrawArrays(GL_TRIANGLES, 0, 6);


    // Raymarching step
	_raymarchingbuffer.shader->use();
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glBindVertexArray(_raymarchingbuffer.VAO);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, _raymarchingbuffer.texture);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glDisable(GL_BLEND);
}

void Renderer::writeImg(const std::uint32_t iteration) const
{
    GLsizei nbChannels = 3;
    GLsizei stride = nbChannels * _windowInfos.x;
    stride += (stride % 4) ? (4 - stride % 4) : 0;
    GLsizei bufferSize = stride * _windowInfos.y;
    std::vector<std::uint8_t> buffer(bufferSize);
    glPixelStorei(GL_PACK_ALIGNMENT, 4);
    glReadBuffer(GL_FRONT);
    glReadPixels(0, 0, _windowInfos.x, _windowInfos.y, GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
    stbi_flip_vertically_on_write(true);
    std::string path = "result/";
    path += std::to_string(iteration);
    path += ".png";
    stbi_write_png(path.c_str(), _windowInfos.x, _windowInfos.y, nbChannels, buffer.data(), stride);
    stbi_write_png("last.png", _windowInfos.x, _windowInfos.y, nbChannels, buffer.data(), stride);

}

void Renderer::drawMesh(Mesh& mesh) const
{
    GLenum renderMode = GL_TRIANGLES;
    switch (mesh.renderMode)
    {
        case LINES:
            renderMode = GL_LINES;
            break;
        case TRIANGLES:
            renderMode = GL_TRIANGLES;
            break;
    }
    glBindVertexArray(mesh.VAO);
    glDrawElements(renderMode, mesh.indices.size(), GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
}

void Renderer::_initFrameBuffer(FrameBuffer& framebuffer, std::string vert, std::string frag)
{
    float quadVertices[] =
	{
		-1.0f,  1.0f,  0.0f, 1.0f,
		-1.0f, -1.0f,  0.0f, 0.0f,
		 1.0f, -1.0f,  1.0f, 0.0f,
		-1.0f,  1.0f,  0.0f, 1.0f,
		 1.0f, -1.0f,  1.0f, 0.0f,
		 1.0f,  1.0f,  1.0f, 1.0f
	};
	glGenFramebuffers(1, &framebuffer.FBO);
	glGenTextures(1, &framebuffer.texture);
	glGenRenderbuffers(1, &framebuffer.RBO);
	glGenVertexArrays(1, &framebuffer.VAO);
	glGenBuffers(1, &framebuffer.VBO);
	glBindVertexArray(framebuffer.VAO);
	glBindBuffer(GL_ARRAY_BUFFER, framebuffer.VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)0);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(2*sizeof(float)));
	glBindVertexArray(0);
    Shader shaderProgram {};
    shaderProgram.setVert(vert);
    shaderProgram.setFrag(frag);
	framebuffer.shader = std::make_shared<Shader>(shaderProgram);

    // Do this at each resize
    framebuffer.shader->use();
	//framebuffer.shader->set1i("screenTexture", 0);
    
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffer.FBO);
	// Color attachment texture
	glBindTexture(GL_TEXTURE_2D, framebuffer.texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, _windowInfos.x, _windowInfos.y, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffer.texture, 0);

	// Renderbuffer for depth and stencil
	glBindRenderbuffer(GL_RENDERBUFFER, framebuffer.RBO);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, _windowInfos.x, _windowInfos.y);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, framebuffer.RBO);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		ERROR("Framebuffer is not complete");
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, _windowInfos.x, _windowInfos.y);
}

void Renderer::applyMaterial(Material& material, Camera& camera, Transform& transform) const
{
    auto& shader = material.shader;
    if (!material.noShader)
    {
        shader.use();
    }
    if (material.hasTexture)
    {
        if (!material.is2D)
        {
            glBindTexture(GL_TEXTURE_3D, material.texture);
            _raymarchingbuffer.shader->use();
            _raymarchingbuffer.shader->set2f("u_resolution", {_windowInfos.x, _windowInfos.y});
            _raymarchingbuffer.shader->set3f("u_eyePos", camera.transform.position);
            _raymarchingbuffer.shader->set3f("u_eyeFront", camera.front);
            _raymarchingbuffer.shader->set1f("u_eyeFOV", camera.FOV);
            _raymarchingbuffer.shader->set1f("u_absorption", material.absorption);
            _raymarchingbuffer.shader->set3f("u_lightIntensity", material.lightIntensity);
        }
        else
        {
            glBindTexture(GL_TEXTURE_2D, material.texture);
        }
    }
    if (material.noShader)
    {
        return;
    }

    glm::mat4 model {1.0f};
    model = glm::translate(model, transform.position);
    model = glm::rotate(model, glm::radians(transform.rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
    model = glm::rotate(model, glm::radians(transform.rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
    model = glm::rotate(model, glm::radians(transform.rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
    model = glm::scale(model, glm::vec3{transform.scale});
    const glm::mat4 projection = glm::infinitePerspective(glm::radians(camera.FOV), static_cast<float>(_windowInfos.x)/_windowInfos.y, 0.1f);
    const glm::mat4 view = glm::lookAt(camera.transform.position, camera.transform.position+camera.front, camera.up);


    shader.setMat4("model", model);
    shader.setMat4("view", view);
    shader.setMat4("projection", projection);
}

void Renderer::initMesh(Mesh& mesh) const
{
    glGenVertexArrays(1, &mesh.VAO);

    glGenBuffers(1, &mesh.VBO);
    glGenBuffers(1, &mesh.NBO);
    glGenBuffers(1, &mesh.EBO);
    glGenBuffers(1, &mesh.TBO);

    glBindVertexArray(mesh.VAO);
        glBindBuffer(GL_ARRAY_BUFFER, mesh.VBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			mesh.vertices.size()*sizeof(GLfloat),
			mesh.vertices.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(0);

		glBindBuffer(GL_ARRAY_BUFFER, mesh.NBO);
		glBufferData(
			GL_ARRAY_BUFFER,
			mesh.normals.size()*sizeof(GLfloat),
			mesh.normals.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
		glEnableVertexAttribArray(1);

		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.EBO);
		glBufferData(
			GL_ELEMENT_ARRAY_BUFFER,
			mesh.indices.size()*sizeof(GLuint),
			mesh.indices.data(),
			GL_STATIC_DRAW
		);
		glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLuint), (GLvoid*)nullptr);
		glEnableVertexAttribArray(2);

        if (!mesh.is2D)
        {
            std::vector<float> texCoords =
            {
                0.0f, 0.0f, 0.0f,
                0.0f, 0.0f, 1.0f,
                1.0f, 0.0f, 1.0f,
                1.00f, 0.0f, 0.0f,
                0.0f, 1.0f, 0.0f,
                0.0f, 1.0f, 1.0f,
                1.0f, 1.0f, 1.0f,
                1.00f, 1.0f, 0.0f,
            };
            glBindBuffer(GL_ARRAY_BUFFER, mesh.TBO);
            glBufferData(
                GL_ARRAY_BUFFER,
                texCoords.size()*sizeof(GLfloat),
                texCoords.data(),
                GL_STATIC_DRAW
            );
            glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), (GLvoid*)nullptr);
            glEnableVertexAttribArray(3);
        }
        else
        {
            std::vector<float> texCoords =
            {
                0.0f, 0.0f, 
                0.0f, 1.0f,
                1.0f, 1.0f,
                1.0f, 0.0f,
            };
            glBindBuffer(GL_ARRAY_BUFFER, mesh.TBO);
            glBufferData(
                GL_ARRAY_BUFFER,
                texCoords.size()*sizeof(GLfloat),
                texCoords.data(),
                GL_STATIC_DRAW
            );
            glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), (GLvoid*)nullptr);
            glEnableVertexAttribArray(3);
        }

    glBindVertexArray(0);
    mesh.initialized = true;
}

void Renderer::initMaterial(Material& material) const
{
    material.hasTexture = true;
    glGenTextures(1, &material.texture);
}

void Renderer::initTexture3D(const std::vector<std::uint8_t>& texture, const std::uint32_t textureGL) const
{
    std::uint32_t size = (std::uint32_t)std::cbrt(static_cast<std::uint32_t>(texture.size()));
    glBindTexture(GL_TEXTURE_3D, textureGL);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_RED, size, size, size, 0, GL_RED, GL_UNSIGNED_BYTE, texture.data());
    glGenerateMipmap(GL_TEXTURE_3D);
}

void Renderer::initTexture2D(const std::vector<std::uint8_t>& texture, const std::uint32_t textureGL) const
{
    std::uint32_t size = (std::uint32_t)sqrt(texture.size()/3);
    glBindTexture(GL_TEXTURE_2D, textureGL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, size, size, 0, GL_RGB, GL_UNSIGNED_BYTE, texture.data());
    glGenerateMipmap(GL_TEXTURE_2D);
}

void Renderer::freeMesh(Mesh& mesh) const
{
    glDeleteVertexArrays(1, &mesh.VAO);
    glDeleteBuffers(1, &mesh.VBO);
    glDeleteBuffers(1, &mesh.NBO);
    glDeleteBuffers(1, &mesh.EBO);
    glDeleteBuffers(1, &mesh.TBO);
}
