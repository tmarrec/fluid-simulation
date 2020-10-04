#pragma once

#include "System.h"

#include <unordered_map>

class Renderer : public System
{
public:
	Renderer(MsgBus_ptr messageBus);
	void cout(std::string string) const override;
	void handleMessage(Message & msg) override;

	void initGl(int width, int height) const;
	void resizeGl(int width, int height) const;
	void initDrawable(std::uint64_t componentID, Shape shape);
	void freeDrawable(std::uint64_t componentID, Shape shape);
	void draw(std::uint64_t componentID, Shape shape);

private:
	std::unordered_map<std::uint64_t, GLObjects> _GLObjects; 

};
