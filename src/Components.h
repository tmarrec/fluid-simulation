#pragma once

#include "ECS.h"

class PositionComponent : public Component
{
private:
	float _xPos;
	float _yPos;

public:
	float x() const
	{
		return _xPos;
	}

	float y() const
	{
		return _yPos;
	}

	void init() override
	{
		_xPos = 0;
		_yPos = 0;
	}

	void update() override
	{
		_xPos++;
		_yPos++;
	}

	void setPos(float x, float y)
	{
		_xPos = x;
		_yPos = y;
	}
	
};
