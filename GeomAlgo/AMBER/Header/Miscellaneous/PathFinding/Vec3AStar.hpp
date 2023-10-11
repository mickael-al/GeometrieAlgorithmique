#ifndef __VEC3_ASTAR__
#define __VEC3_ASTAR__

#include "AStar.hpp"
#include "glm/glm.hpp"

class Vec3AStar : public AStarNode
{
public:
	glm::vec3 position;
public:	
	Vec3AStar(){}
	~Vec3AStar(){}
	float distanceTo(AStarNode* node) const
	{
		Vec3AStar * nodeC = (Vec3AStar*)node;

		return glm::distance(position, nodeC->position);
	}	
};

#endif // !__VEC3_ASTAR__
