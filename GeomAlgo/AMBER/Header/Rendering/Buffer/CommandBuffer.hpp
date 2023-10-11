#ifndef __ENGINE_COMMAND_BUFFER__
#define __ENGINE_COMMAND_BUFFER__

#include "BufferManager.hpp"
#include "PointeurClass.hpp"
#include "Model.hpp"
#include "PostProcessing.hpp"
#include "ShadowManager.hpp"
#include "VulkanMisc.hpp"
#include <vector>

namespace Ge
{
    class CommandBuffer
    {
    private:
        friend class RenderingEngine;
        bool initialize(ShadowManager* shadowManager, PostProcessing postProcessing,VulkanMisc * vM, ptrClass * ptrC);
        VmaBuffer createInstanceBuffer(std::vector<Model*> models);
        void release();
    private:
        VulkanMisc * vulkanM;
		std::vector<VkCommandBuffer> m_commandBuffers;
        std::vector<VmaBuffer> m_instancedBuffer;
        VmaBuffer m_instancedBufferSkybox;
    };
}

#endif //__ENGINE_COMMAND_BUFFER__