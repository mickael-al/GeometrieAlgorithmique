#ifndef __COMPUTE_SHADER_MANAGER__
#define __COMPUTE_SHADER_MANAGER__

#include "VulkanMisc.hpp"
#include "ComputeBuffer.hpp"
#include "ComputeImage.hpp"
#include "ComputeShader.hpp"

namespace Ge
{
	class ComputeShaderManager final
	{
	public:
		bool initialize(VulkanMisc* vM);
		void release();
		ComputeBuffer* createComputeBuffer(int baseCount, int sizeofR);
		ComputeImage* createComputeImage(uint32_t width, uint32_t height, VkFormat format);
		ComputeShader* createComputeShader(const std::string& shaderPath, const std::vector<ComputeData*>& buffers);
	private:
		VulkanMisc* vulkanM;
	};
}

#endif //!__COMPUTE_SHADER_MANAGER__