#include "ComputeShaderManager.hpp"

namespace Ge
{
	bool ComputeShaderManager::initialize(VulkanMisc* vM)
	{
		vulkanM = vM;
		return true;
	}

	void ComputeShaderManager::release()
	{

	}

	ComputeBuffer* ComputeShaderManager::createComputeBuffer(int baseCount, int sizeofR)
	{
		ComputeBuffer* cb = new ComputeBuffer(vulkanM, baseCount, sizeofR);
		return cb;
	}

	ComputeImage* ComputeShaderManager::createComputeImage(uint32_t width, uint32_t height, VkFormat format)
	{
		ComputeImage * ci = new ComputeImage(vulkanM, width, height, format);
		return ci;
	}

	ComputeShader* ComputeShaderManager::createComputeShader(const std::string& shaderPath, const std::vector<ComputeData*>& buffers)
	{
		ComputeShader* cs = new ComputeShader(vulkanM, shaderPath, buffers);
		return cs;
	}
}