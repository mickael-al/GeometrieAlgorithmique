#include "ComputeImage.hpp"

namespace Ge
{
	ComputeImage::ComputeImage(VulkanMisc* vM, uint32_t width, uint32_t height, VkFormat format)
	{
		vulkanM = vM;

		VkDescriptorSetLayoutBinding uboLayoutBinding{};
		uboLayoutBinding.binding = 0;
		uboLayoutBinding.descriptorCount = 1;
		uboLayoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		uboLayoutBinding.pImmutableSamplers = nullptr;
		uboLayoutBinding.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

		VkDescriptorSetLayoutCreateInfo layoutInfo{};
		layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
		layoutInfo.bindingCount = 1;
		layoutInfo.pBindings = &uboLayoutBinding;

		if (vkCreateDescriptorSetLayout(vM->str_VulkanDeviceMisc->str_device, &layoutInfo, nullptr, &m_DescriptorSetLayout) != VK_SUCCESS)
		{
			Debug::Error("Echec de la creation du descriptor set layout");
		}

		VkDescriptorPoolSize poolSize{};
		poolSize.type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		poolSize.descriptorCount = 1;

		VkDescriptorPoolCreateInfo poolInfo{};
		poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
		poolInfo.poolSizeCount = 1;
		poolInfo.pPoolSizes = &poolSize;
		poolInfo.maxSets = static_cast<uint32_t>(vM->str_VulkanSwapChainMisc->str_swapChainImages.size());

		if (vkCreateDescriptorPool(vM->str_VulkanDeviceMisc->str_device, &poolInfo, nullptr, &m_DescriptorPool) != VK_SUCCESS)
		{
			Debug::Error("Echec de la creation d'un descriptor pool");
		}

		VkDescriptorSetAllocateInfo allocInfo{};
		allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
		allocInfo.descriptorPool = m_DescriptorPool;
		allocInfo.descriptorSetCount = 1;
		allocInfo.pSetLayouts = &m_DescriptorSetLayout;

		if (vkAllocateDescriptorSets(vM->str_VulkanDeviceMisc->str_device, &allocInfo, &m_DescriptorSets) != VK_SUCCESS)
		{
			Debug::Error("Echec de l'allocation du descriptor sets");
		}

		if (!BufferManager::createImageBuffer(width, height, VK_IMAGE_TYPE_2D, 1, 1, VK_SAMPLE_COUNT_1_BIT, format, VK_IMAGE_TILING_OPTIMAL, VK_IMAGE_USAGE_TRANSFER_SRC_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT | VK_IMAGE_USAGE_STORAGE_BIT | VK_IMAGE_USAGE_SAMPLED_BIT, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, m_vmaImageBuffer, 0, vM))
		{
			Debug::Error("Echec de la creation d'un image buffer");
		}

		VkImageViewCreateInfo viewInfo{};
		viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
		viewInfo.image = m_vmaImageBuffer.image;
		viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
		viewInfo.format = format;
		//viewInfo.components = { VK_COMPONENT_SWIZZLE_R };
		viewInfo.subresourceRange = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 };
		viewInfo.subresourceRange.layerCount = 1;

		if (vkCreateImageView(vM->str_VulkanDeviceMisc->str_device, &viewInfo, nullptr, &m_imageView) != VK_SUCCESS)
		{
			Debug::Error("Echec de la creation d'une image vue");
		}

		VkDescriptorImageInfo imageInfo{};
		imageInfo.imageLayout = VK_IMAGE_LAYOUT_GENERAL; // VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL
		imageInfo.imageView = m_imageView;

		VkWriteDescriptorSet descriptorWrites{};
		descriptorWrites.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
		descriptorWrites.dstSet = m_DescriptorSets;
		descriptorWrites.dstBinding = 0;
		descriptorWrites.dstArrayElement = 0;
		descriptorWrites.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		descriptorWrites.descriptorCount = 1;
		descriptorWrites.pImageInfo = &imageInfo;

		vkUpdateDescriptorSets(vM->str_VulkanDeviceMisc->str_device, 1, &descriptorWrites, 0, nullptr);
		outImage = false;
	}

	ComputeImage::ComputeImage(VulkanMisc* vM, VkImageView image)
	{
		vulkanM = vM;

		VkDescriptorSetLayoutBinding uboLayoutBinding{};
		uboLayoutBinding.binding = 0;
		uboLayoutBinding.descriptorCount = 1;
		uboLayoutBinding.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		uboLayoutBinding.pImmutableSamplers = nullptr;
		uboLayoutBinding.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;

		VkDescriptorSetLayoutCreateInfo layoutInfo{};
		layoutInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
		layoutInfo.bindingCount = 1;
		layoutInfo.pBindings = &uboLayoutBinding;

		if (vkCreateDescriptorSetLayout(vM->str_VulkanDeviceMisc->str_device, &layoutInfo, nullptr, &m_DescriptorSetLayout) != VK_SUCCESS)
		{
			Debug::Error("Echec de la creation du descriptor set layout");
		}

		VkDescriptorPoolSize poolSize{};
		poolSize.type = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		poolSize.descriptorCount = 1;

		VkDescriptorPoolCreateInfo poolInfo{};
		poolInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
		poolInfo.poolSizeCount = 1;
		poolInfo.pPoolSizes = &poolSize;
		poolInfo.maxSets = static_cast<uint32_t>(vM->str_VulkanSwapChainMisc->str_swapChainImages.size());

		if (vkCreateDescriptorPool(vM->str_VulkanDeviceMisc->str_device, &poolInfo, nullptr, &m_DescriptorPool) != VK_SUCCESS)
		{
			Debug::Error("Echec de la creation d'un descriptor pool");
		}

		VkDescriptorSetAllocateInfo allocInfo{};
		allocInfo.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
		allocInfo.descriptorPool = m_DescriptorPool;
		allocInfo.descriptorSetCount = 1;
		allocInfo.pSetLayouts = &m_DescriptorSetLayout;

		if (vkAllocateDescriptorSets(vM->str_VulkanDeviceMisc->str_device, &allocInfo, &m_DescriptorSets) != VK_SUCCESS)
		{
			Debug::Error("Echec de l'allocation du descriptor sets");
		}

		VkDescriptorImageInfo imageInfo{};
		imageInfo.imageLayout = VK_IMAGE_LAYOUT_GENERAL; // VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL
		imageInfo.imageView = image;

		VkWriteDescriptorSet descriptorWrites{};
		descriptorWrites.sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
		descriptorWrites.dstSet = m_DescriptorSets;
		descriptorWrites.dstBinding = 0;
		descriptorWrites.dstArrayElement = 0;
		descriptorWrites.descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_IMAGE;
		descriptorWrites.descriptorCount = 1;
		descriptorWrites.pImageInfo = &imageInfo;

		vkUpdateDescriptorSets(vM->str_VulkanDeviceMisc->str_device, 1, &descriptorWrites, 0, nullptr);
		outImage = true;
	}

	VkDescriptorSetLayout ComputeImage::getDescriptorSetLayout()
	{
		return m_DescriptorSetLayout;
	}

	VkDescriptorSet ComputeImage::getDescriptorSet()
	{
		return m_DescriptorSets;
	}

	VmaBufferImage ComputeImage::getVmaBufferImage()
	{
		return m_vmaImageBuffer;
	}

	ComputeImage::~ComputeImage()
	{
		if (!outImage)
		{
			vkDestroyImageView(vulkanM->str_VulkanDeviceMisc->str_device, m_imageView, nullptr);
			BufferManager::destroyImageBuffer(m_vmaImageBuffer);
		}
		vkDestroyDescriptorPool(vulkanM->str_VulkanDeviceMisc->str_device, m_DescriptorPool, nullptr);
		vkDestroyDescriptorSetLayout(vulkanM->str_VulkanDeviceMisc->str_device, m_DescriptorSetLayout, nullptr);
	}
}