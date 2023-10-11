#ifndef __ENGINE_MATERIALS__
#define __ENGINE_MATERIALS__

#include "Textures.hpp"
#include "Debug.hpp"
#include "VulkanMisc.hpp"
#include "BufferManager.hpp"
#include "UniformBufferMaterial.hpp"
#include "imgui-cmake/Header/imgui.h"
#include "GraphiquePipeline.hpp"
#include "Component.hpp"

namespace Ge
{
	class Materials : Component
	{
	public:
		Materials(int index, VulkanMisc * vM);
		void setColor(glm::vec3 color);
		void setMetallic(float metal);
		void setRoughness(float roughness);
		void setNormal(float normal);
		void setOclusion(float ao);
		void setAlbedoTexture(Textures * albedo);
		void setNormalTexture(Textures * normal);
		void setMetallicTexture(Textures * metallic);
		void setRoughnessTexture(Textures * roughness);
		void setOclusionTexture(Textures * oclu);
		glm::vec3 getColor() const;
		float getMetallic() const;
		float getRoughness() const;
		float getNormal() const;
		float getOclusion() const;
		Textures * getAlbedoTexture() const;
		Textures * getNormalTexture() const;
		Textures * getMetallicTexture() const;
		Textures * getRoughnessTexture() const;
		Textures * getOclusionTexture() const;
		VkBuffer getUniformBuffers() const;
		void updateUniformBufferMaterial();
		int getIndex() const;
		void setIndex(int i);
		void majTextureIndex();
		bool getShadowCast() const;
		void setShadowCast(bool state);
		glm::vec2 getOffset() const;
		glm::vec2 getTilling() const;
		void setOffset(glm::vec2 off);
		void setTilling(glm::vec2 tilling);
		void setPipeline(GraphiquePipeline * p);
		void setPipelineIndex(int p);
		int getPipelineIndex() const;
		void onGUI();
		~Materials();
	private:
		UniformBufferMaterial m_ubm{};
		VulkanMisc * vulkanM;
		VmaBuffer m_vmaUniformBuffer;
		Textures * m_albedoMap;
		Textures * m_normalMap;
		Textures * m_metallicMap;
		Textures * m_RoughnessMap;
		Textures * m_aoMap;
		int m_pipelineIndex;
		float m_color[3];
		float m_offset[2];
		float m_tilling[2];
		int m_index = 0;
	};
}

#endif//__ENGINE_MATERIALS__