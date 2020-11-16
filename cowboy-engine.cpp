#include <iostream>

#include "Coordinator.h"
#include "utils.h"

Coordinator ECS;

int main()
{
	PRINT_TITLE();

	ECS.Init();

	return EXIT_SUCCESS;
}