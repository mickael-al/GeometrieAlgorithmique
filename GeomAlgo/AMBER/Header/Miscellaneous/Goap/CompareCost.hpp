#ifndef __COMPARE_COST__
#define __COMPARE_COST__

#include "GNode.hpp"

struct CompareCost
{
	bool operator()(const GNode* const a, const GNode* const b)
	{
		if (a->openNode.size() == b->openNode.size())
		{
			return a->cost > b->cost;
		}
		return a->openNode.size() > b->openNode.size();
	}
};

#endif //!__COMPARE_COST__