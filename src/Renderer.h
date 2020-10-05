#pragma once

#include "System.h"
#include "Shader.h"

class Renderer : public System
{
public:
	Renderer(MsgBus_ptr messageBus);
	void cout(std::string string) const override;
	void handleMessage(Message & msg) override;

	void initGl(int width, int height) const;
	void resizeGl(int width, int height) const;

	void initDrawable(Message & msg);
	void freeDrawable(Message & msg);
	void draw(Message & msg);

private:
	void _useShader(Message & msg);
	glm::mat4 _getModel(Message & msg);
};

