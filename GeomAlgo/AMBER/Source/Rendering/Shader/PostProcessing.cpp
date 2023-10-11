#include "PostProcessing.hpp"
#include "ShaderUtil.hpp"
#include "Debug.hpp"

namespace Ge
{
    void transitionImageLayout(VkImage image, VkFormat format, VkImageLayout oldLayout, VkImageLayout newLayout, uint32_t mipLevels, uint32_t layerCount, VulkanMisc* vM)
    {
        VkCommandBuffer commandBuffer = BufferManager::beginSingleTimeCommands(vM);

        VkImageMemoryBarrier barrier{};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = oldLayout;
        barrier.newLayout = newLayout;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = mipLevels;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = layerCount;

        VkPipelineStageFlags sourceStage;
        VkPipelineStageFlags destinationStage;

        if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL)
        {
            barrier.srcAccessMask = 0;
            barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
            destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        }
        else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL)
        {
            barrier.srcAccessMask = 0;
            barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;

            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
            destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        }
        else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_GENERAL)
        {
            barrier.srcAccessMask = 0;
            barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
            destinationStage = VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT;
        }
        else if (oldLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL)
        {
            barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
            barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

            sourceStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
            destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
        }
        else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED && newLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL)
        {
            barrier.srcAccessMask = 0;
            barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

            sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
            destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        }
        else if (oldLayout == VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL && newLayout == VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL)
        {
            barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
            barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;

            sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
            destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        }
        else
        {
            Debug::Error("Transition d'organisation non supportée !");
            return;
        }

        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );

        BufferManager::endSingleTimeCommands(commandBuffer, vM);
    }

    bool PostProcessing::initialize(VulkanMisc* vM)
    {
        vulkanM = vM;

        m_sizeBufferData = new ComputeBuffer(vM, 2, sizeof(uint32_t));
        m_bloomBlur = new ComputeImage(vM, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, VkFormat::VK_FORMAT_B8G8R8A8_UNORM);
        transitionImageLayout(m_bloomBlur->getVmaBufferImage().image, VK_FORMAT_B8G8R8A8_UNORM, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL, 1, 1, vM);
        m_blur = new ComputeImage(vM, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, VkFormat::VK_FORMAT_B8G8R8A8_UNORM);
        transitionImageLayout(m_blur->getVmaBufferImage().image, VK_FORMAT_B8G8R8A8_UNORM, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL, 1, 1, vM);

        m_imageBase = new ComputeImage(vulkanM, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, VkFormat::VK_FORMAT_B8G8R8A8_UNORM);
        transitionImageLayout(m_imageBase->getVmaBufferImage().image, VK_FORMAT_B8G8R8A8_UNORM, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL, 1, 1, vM);

        m_chromaticAberration = new ComputeImage(vulkanM, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, VkFormat::VK_FORMAT_B8G8R8A8_UNORM);
        transitionImageLayout(m_chromaticAberration->getVmaBufferImage().image, VK_FORMAT_B8G8R8A8_UNORM, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_GENERAL, 1, 1, vM);

        //m_depthImage = new ComputeImage(vulkanM, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width, vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, VkFormat::VK_FORMAT_R32_SFLOAT);
        //transitionImageLayout(m_depthImage->getVmaBufferImage().image, VK_FORMAT_R32_SFLOAT, VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, 1, vM);

        std::vector<ComputeData*> cd;
        cd.push_back(m_imageBase);
        cd.push_back(m_bloomBlur);
        cd.push_back(m_sizeBufferData);
        m_cshaderTone = new ComputeShader(vulkanM, "../Shader/tonemapping.comp.spv", cd);

        std::vector<ComputeData*> bloomcd;
        bloomcd.push_back(m_bloomBlur);
        bloomcd.push_back(m_blur);
        bloomcd.push_back(m_sizeBufferData);
        m_cshaderBlur = new ComputeShader(vulkanM, "../Shader/bloom.comp.spv", bloomcd);


        std::vector<ComputeData*> mixcd;
        mixcd.push_back(m_blur);
        mixcd.push_back(m_imageBase);
        mixcd.push_back(m_sizeBufferData);
        m_cshaderMix = new ComputeShader(vulkanM, "../Shader/bloom_mix.comp.spv", mixcd);

        std::vector<ComputeData*> cacd;
        cacd.push_back(m_imageBase);
        cacd.push_back(m_chromaticAberration);
        cacd.push_back(m_sizeBufferData);
        m_cshaderChromaticAberration = new ComputeShader(vulkanM, "../Shader/chromatic_aberation.comp.spv", cacd);

        std::vector<ComputeData*> dofcd;
        dofcd.push_back(m_chromaticAberration);
        dofcd.push_back(m_imageBase);
        dofcd.push_back(m_sizeBufferData);

        m_cshaderDepthOfField = new ComputeShader(vulkanM, "../Shader/depth_of_field.comp.spv", dofcd);
        return true;
    }

    void PostProcessing::execute(VkCommandBuffer commandBuffer, uint32_t index)
    {
        m_sizeBufferData->SetData((void*)(&vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent), 0, 0, 2);
        VkImageMemoryBarrier barrier{};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = vulkanM->str_VulkanCommandeBufferMisc->str_colorImage;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;
        VkPipelineStageFlags sourceStage;
        VkPipelineStageFlags destinationStage;
        barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;

        sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );

        VkImageResolve imageResolve = {};
        imageResolve.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageResolve.srcSubresource.mipLevel = 0;
        imageResolve.srcSubresource.baseArrayLayer = 0;
        imageResolve.srcSubresource.layerCount = 1;
        imageResolve.srcOffset = { 0, 0, 0 };
        imageResolve.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageResolve.dstSubresource.mipLevel = 0;
        imageResolve.dstSubresource.baseArrayLayer = 0;
        imageResolve.dstSubresource.layerCount = 1;
        imageResolve.dstOffset = { 0, 0, 0 };
        imageResolve.extent.width = vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width;
        imageResolve.extent.height = vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height;
        imageResolve.extent.depth = 1;

        vkCmdResolveImage(
            commandBuffer,
            vulkanM->str_VulkanCommandeBufferMisc->str_colorImage, // The multisampled image (mss8)
            VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
            m_imageBase->getVmaBufferImage().image, // The destination image (mss1)
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            1,
            &imageResolve
        );

        int x, y, z;
        ShaderUtil::CalcWorkSize(vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width * vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height, &x, &y, &z);
        m_cshaderTone->dispatch(commandBuffer, x, y, z);
        m_cshaderBlur->dispatch(commandBuffer, x, y, z);
        m_cshaderBlur->swapBuffer(0, 1);
        m_cshaderBlur->dispatch(commandBuffer, x, y, z);
        m_cshaderBlur->swapBuffer(0, 1);

        m_cshaderMix->dispatch(commandBuffer, x, y, z);
        m_cshaderChromaticAberration->dispatch(commandBuffer, x, y, z);
        m_cshaderDepthOfField->dispatch(commandBuffer, x, y, z);

        barrier = {};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = vulkanM->str_VulkanSwapChainMisc->str_swapChainImages[index];
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;

        barrier.srcAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

        sourceStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );



        barrier = {};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_GENERAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = m_imageBase->getVmaBufferImage().image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;

        barrier.srcAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_READ_BIT;

        sourceStage = VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT;
        destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );


        VkImageCopy copyRegion = {};
        copyRegion.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        copyRegion.srcSubresource.mipLevel = 0;
        copyRegion.srcSubresource.baseArrayLayer = 0;
        copyRegion.srcSubresource.layerCount = 1;
        copyRegion.srcOffset = { 0, 0, 0 };

        copyRegion.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        copyRegion.dstSubresource.mipLevel = 0;
        copyRegion.dstSubresource.baseArrayLayer = 0;
        copyRegion.dstSubresource.layerCount = 1;
        copyRegion.dstOffset = { 0, 0, 0 };

        copyRegion.extent.width = vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.width;
        copyRegion.extent.height = vulkanM->str_VulkanSwapChainMisc->str_swapChainExtent.height;
        copyRegion.extent.depth = 1;

        vkCmdCopyImage(
            commandBuffer,
            m_imageBase->getVmaBufferImage().image, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
            vulkanM->str_VulkanSwapChainMisc->str_swapChainImages[index], VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            1,
            &copyRegion
        );



        barrier = {};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = vulkanM->str_VulkanSwapChainMisc->str_swapChainImages[index];
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;

        barrier.srcAccessMask = VK_PIPELINE_STAGE_TRANSFER_BIT;
        barrier.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

        sourceStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        destinationStage = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );



        barrier = {};
        barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        barrier.oldLayout = VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL;
        barrier.newLayout = VK_IMAGE_LAYOUT_GENERAL;
        barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        barrier.image = m_imageBase->getVmaBufferImage().image;
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        barrier.subresourceRange.baseMipLevel = 0;
        barrier.subresourceRange.levelCount = 1;
        barrier.subresourceRange.baseArrayLayer = 0;
        barrier.subresourceRange.layerCount = 1;

        barrier.srcAccessMask = VK_ACCESS_TRANSFER_READ_BIT;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;

        sourceStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        destinationStage = VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT;
        vkCmdPipelineBarrier(
            commandBuffer,
            sourceStage, destinationStage,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );
    }

    void PostProcessing::release()
    {
        delete m_imageBase;
        delete m_chromaticAberration;
        delete m_cshaderTone;
        delete m_cshaderBlur;
        delete m_cshaderMix;
        delete m_cshaderChromaticAberration;
        delete m_cshaderDepthOfField;
        //delete m_depthImage;
        delete m_sizeBufferData;
        delete m_bloomBlur;
        delete m_blur;
    }
}