#ifndef __COMPUTE_IMAGE__
#define __COMPUTE_IMAGE__

#include "VulkanMisc.hpp"
#include "Debug.hpp"
#include "BufferManager.hpp"
#include "ComputeData.hpp"

namespace Ge
{
	class ComputeImage : public ComputeData
	{
	public:
		ComputeImage(VulkanMisc* vM, uint32_t width, uint32_t height, VkFormat format);
		ComputeImage(VulkanMisc* vM, VkImageView image);
		VkDescriptorSetLayout getDescriptorSetLayout();
		VkDescriptorSet getDescriptorSet();
		VmaBufferImage getVmaBufferImage();
		~ComputeImage();
	private:
		VulkanMisc* vulkanM;
		VkDescriptorSetLayout m_DescriptorSetLayout;
		VkDescriptorPool m_DescriptorPool;
		VkDescriptorSet m_DescriptorSets;
		VkDescriptorType m_DescriptorType;
		VmaBufferImage m_vmaImageBuffer;
		VkImageView m_imageView;
		bool outImage = false;
	};
}

#endif //!__COMPUTE_BUFFER__