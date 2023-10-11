#ifndef __G_NODE__
#define __G_NODE__

#include <vector>
#include <unordered_map>
#include <string>
class GAction;

struct GNode
{
	std::vector<const GAction*> openNode;
	const GAction* action;
	std::unordered_map<std::string, float> ressource;
	unsigned short cost = 0;
	GNode* parent;
};

#endif //!__G_NODE__