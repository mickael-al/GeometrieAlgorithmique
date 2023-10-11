#ifndef __BLOCK__
#define __BLOCK__

#include "MetaData.hpp"

class Block
{
public:
	Block(ptrClass pc,LoadPrefab*lp,glm::vec3 position = glm::vec3(0.0f),glm::vec3 euler = glm::vec3(0.0f), glm::vec3 scale = glm::vec3(1.0f),bool collider = true);
	Model* getMainModel() const;
	std::vector<CollisionBody*> getCollision() const;
	void UpdateTransform();
	void DebugCollider();
	~Block();
private:
	ptrClass m_pc;
	LoadPrefab* m_lp;
	std::vector<Model*> m_model;
	std::vector<CollisionBody*> m_collision;
};

#endif //!__BLOCK__