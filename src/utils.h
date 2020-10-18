#pragma once

#include <iostream>

#include <glm/glm/glm.hpp>
#include <glm/glm/gtc/matrix_transform.hpp>
#include <glm/glm/gtc/type_ptr.hpp>
#include <glm/glm/gtx/string_cast.hpp>
#include <glm/glm/gtx/rotate_vector.hpp>

#include <GL/gl.h>

#include <vector>
#include <algorithm>

struct Range
{
	double start;
	double end;
};

#define WARNING(m) \
	std::cerr << m << std::endl;

#define ERROR(m) \
	std::cerr << m << std::endl; \
	exit(2);

#ifndef NDEBUG
#define ASSERT(c, m) \
	if (!(c)) { \
		std::cerr << std::endl; \
		std::cerr << "\033[41m\033[1m[ASSERT ERROR]\033[49m\033[0m " << std::endl; \
		std::cerr << "\033[1mFILE\033[0m    : " << __FILE__ << std::endl; \
		std::cerr << "\033[1mFUNCTION\033[0m: " << __func__ << std::endl; \
		std::cerr << "\033[1mLINE\033[0m    : " << __LINE__ << std::endl; \
		std::cerr << "\033[1mMESSAGE\033[0m : " << m << std::endl; \
		std::cerr << std::endl; \
		exit(3); \
	}
#else
#define ASSERT(c, m) \
		do { \
		} while(0)
#endif

#define PRINT_TITLE() \
	std::cout << std::endl << "\033[96m\033[1m"; \
	std::cout << "                          __                                                                            \n                         /\\ \\                                                      __                   \n  ___    ___   __  __  __\\ \\ \\____    ___   __  __              __    ___      __ /\\_\\    ___      __   \n /'___\\ / __`\\/\\ \\/\\ \\/\\ \\\\ \\ '__`\\  / __`\\/\\ \\/\\ \\  _______  /'__`\\/' _ `\\  /'_ `\\/\\ \\ /' _ `\\  /'__`\\ \n/\\ \\__//\\ \\L\\ \\ \\ \\_/ \\_/ \\\\ \\ \\L\\ \\/\\ \\L\\ \\ \\ \\_\\ \\/\\______\\/\\  __//\\ \\/\\ \\/\\ \\L\\ \\ \\ \\/\\ \\/\\ \\/\\  __/ \n\\ \\____\\ \\____/\\ \\___x___/' \\ \\_,__/\\ \\____/\\/`____ \\/______/\\ \\____\\ \\_\\ \\_\\ \\____ \\ \\_\\ \\_\\ \\_\\ \\____\\\n \\/____/\\/___/  \\/__//__/    \\/___/  \\/___/  `/___/> \\        \\/____/\\/_/\\/_/\\/___L\\ \\/_/\\/_/\\/_/\\/____/\n                                                /\\___/                         /\\____/                  \n                                                \\/__/                          \\_/__/                   \n" << std::endl; \
	std::cout << "\033[39m\033[0m"; \
	for (unsigned char i = 0; i < 3; ++i) \
	{ \
		std::cout << std::endl; \
	}
