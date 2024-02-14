#ifndef __BOUDING_BOX__
#define __BOUDING_BOX__

#include <string>

namespace Ge
{
	class Model;
	class Materials;
}
using namespace Ge;

struct BoundingBox
{
	Model* box;
	Materials* material;
	std::string name;
};

#include <vector>
#include "glm/glm.hpp"
#include <glm/gtc/quaternion.hpp>
struct BoneNode
{
	Model* segments;
	std::string* name;
	BoundingBox* bb;
	glm::vec3 A;
	glm::vec3 B;
	glm::vec3 position;
	glm::vec3 eulerAngle;
	std::vector<BoneNode*> child;
	Model* root_vertex;
	bool parent = true;
};

#endif //!__BOUDING_BOX__