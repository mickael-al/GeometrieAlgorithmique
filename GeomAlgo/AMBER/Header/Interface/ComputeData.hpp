#ifndef __COMPUTE_DATA__
#define __COMPUTE_DATA__

#include "vulkan/vulkan.h"

class ComputeData
{
public:
	virtual VkDescriptorSetLayout getDescriptorSetLayout() = 0;
	virtual VkDescriptorSet getDescriptorSet() = 0;
};

#endif //!__COMPUTE_DATA__