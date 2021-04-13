#include "Renderer.h"

void Renderer::init(std::shared_ptr<Window> window)
{
	_window = window;

    if (!gladLoadGLLoader(reinterpret_cast<GLADloadproc>(glfwGetProcAddress)))
    {
        ERROR("Failed to initialize glad");
    }

    const WindowInfos windowInfos = _window->windowInfos();
    glViewport(0, 0, windowInfos.x, windowInfos.y);
    glEnable(GL_DEPTH_TEST);
}

void Renderer::prePass() const
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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

void Renderer::applyMaterial(Material& material, Camera& camera, Transform& transform) const
{
    if (material.hasTexture)
    {
        glBindTexture(GL_TEXTURE_2D, material.texture);
    }

    glm::mat4 model {1.0f};
    model = glm::translate(model, transform.position);
    model = glm::rotate(model, glm::radians(transform.rotation.x), glm::vec3(1.0f, 0.0f, 0.0f));
    model = glm::rotate(model, glm::radians(transform.rotation.y), glm::vec3(0.0f, 1.0f, 0.0f));
    model = glm::rotate(model, glm::radians(transform.rotation.z), glm::vec3(0.0f, 0.0f, 1.0f));
    model = glm::scale(model, glm::vec3{transform.scale});
    const WindowInfos windowInfos = _window->windowInfos();
    const glm::mat4 projection = glm::infinitePerspective(glm::radians(camera.FOV), static_cast<float>(windowInfos.x)/windowInfos.y, 0.1f);
    const glm::mat4 view = glm::lookAt(camera.transform.position, camera.transform.position+camera.front, camera.up);

    auto& shader = material.shader;
    shader.use();
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

    glBindVertexArray(0);
    mesh.initialized = true;
}

void Renderer::initMaterial(Material& material) const
{
    material.hasTexture = true;
    glGenTextures(1, &material.texture);
    glBindTexture(GL_TEXTURE_2D, material.texture);
}

void Renderer::updateTexture(const std::vector<std::uint8_t>& texture) const
{
    std::uint32_t size = (std::uint32_t)sqrt(texture.size()/3);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glad_glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
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
