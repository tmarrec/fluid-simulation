#pragma once
#include <queue>
#include <array>
#include <iostream>

#include "../types.h"
#include "../utils.h"

class EntityManager
{
public:
	EntityManager()
	{
		// Initialize the queue with all possible entity IDs
		for (Entity entity = 0; entity < MAX_ENTITIES; ++entity)
		{
			mAvailableEntities.push(entity);
		}
	}

	Entity CreateEntity()
	{
		ASSERT(mLivingEntityCount < MAX_ENTITIES, "Too many entities in existence.");

		// Take an ID from the front of the queue
		Entity id = mAvailableEntities.front();
		mAvailableEntities.pop();
		++mLivingEntityCount;

		return id;
	}

	void DestroyEntity(Entity entity)
	{
		ASSERT(entity < MAX_ENTITIES, "Entity out of range.");

		// Invalidate the destroyed entity's signature
		mSignatures[entity].reset();

		// Put the destroyed ID at the back of the queue
		mAvailableEntities.push(entity);
		--mLivingEntityCount;
	}

	void SetSignature(Entity entity, Signature signature)
	{
		ASSERT(entity < MAX_ENTITIES, "Entity out of range.");

		// Put this entity's signature into the array
		mSignatures[entity] = signature;
	}

	Signature GetSignature(Entity entity)
	{
		ASSERT(entity < MAX_ENTITIES, "Entity out of range.");

		// Get this entity's signature from the array
		return mSignatures[entity];
	}

private:
	// Queue of unused entity IDs
	std::queue<Entity> mAvailableEntities {};

	// Array of signatures where the index corresponds to the entity ID
	std::array<Signature, MAX_ENTITIES> mSignatures {};

	// Total living entities - used to keep limits on how many exist
	uint32_t mLivingEntityCount {};
};
