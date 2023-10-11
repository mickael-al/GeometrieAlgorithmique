#ifndef __POST_PROCESSING__
#define __POST_PROCESSING__

#include "VulkanMisc.hpp"
#include "ComputeShader.hpp"
#include "ComputeImage.hpp"

namespace Ge
{
	class PostProcessing final
	{
	public:
		bool initialize(VulkanMisc* vM);
		void execute(VkCommandBuffer commandBuffer, uint32_t index);		
		void release();
	private:
		VulkanMisc* vulkanM;
		ComputeBuffer * m_sizeBufferData;
		ComputeImage* m_bloomBlur;
		ComputeImage* m_blur;
		ComputeImage* m_imageBase;
		ComputeImage* m_chromaticAberration;
		ComputeImage* m_depthImage;
		ComputeShader* m_cshaderTone;
		ComputeShader* m_cshaderBlur;
		ComputeShader* m_cshaderMix;
		ComputeShader* m_cshaderChromaticAberration;
		ComputeShader* m_cshaderDepthOfField;
	};
}

#endif //!__POST_PROCESSING__