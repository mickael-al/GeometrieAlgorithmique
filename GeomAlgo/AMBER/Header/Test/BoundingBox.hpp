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

#endif //!__BOUDING_BOX__