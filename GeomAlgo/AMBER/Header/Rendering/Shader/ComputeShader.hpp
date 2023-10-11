#ifndef __COMPUTE_SHADER__
#define __COMPUTE_SHADER__

#include "VulkanMisc.hpp"
#include <vector>
#include "ComputeBuffer.hpp"

namespace Ge
{
	class ComputeShader
	{
	public:
		ComputeShader(VulkanMisc* vM, const std::string& shaderPath, const std::vector<ComputeData*>& buffers);
		void dispatch(uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);
		void dispatch(VkCommandBuffer commandBuffer,uint32_t groupCountX, uint32_t groupCountY, uint32_t groupCountZ);
		void swapBuffer(size_t index1, size_t index2);
		void setComputeData(size_t index, ComputeData * cd);
		ComputeData* getComputeData(size_t index) const;
		~ComputeShader();
	private:
		static VkShaderModule LoadShader(const std::string& filename, VkDevice device, VulkanMisc* vM);
	private:
		VulkanMisc* vulkanM;
		VkPipeline m_Pipeline;
		VkPipelineLayout m_PipelineLayout;
		std::vector<ComputeData*> m_Buffers;
	};
}

#endif //!__COMPUTE_SHADER__