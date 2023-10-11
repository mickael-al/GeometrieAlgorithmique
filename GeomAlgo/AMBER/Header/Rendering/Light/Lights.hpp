#ifndef __ENGINE_LIGHT__
#define __ENGINE_LIGHT__

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#include "glm/glm.hpp"
#include "glm/gtc/quaternion.hpp"
#include "glm/gtx/euler_angles.hpp"
#include "glm/common.hpp"
#include "UniformBufferLight.hpp"
#include "Debug.hpp"
#include "VulkanMisc.hpp"
#include "BufferManager.hpp"
#include "GObject.hpp"
#include "LightData.hpp"
#include "Shadow.hpp"

namespace Ge
{
	class Lights : public GObject
	{
	public:
		Lights(int index, VulkanMisc *vM);
		void setColors(glm::vec3 color);
		glm::vec3 getColors() const;		
		float getRange() const;
		float getSpotAngle() const;
		int getStatus() const; //Statut directional spotlight pointlight
		int getIndex() const;
		VkBuffer getUniformBuffers() const;
		bool getShadow() const;
		void setRange(float r);
		void setSpotAngle(float r);
		void setIndex(int i);		
		void updateUniformBufferLight();
		void mapMemory();
		void onGUI() override;
		void setShadow(bool state);		
		virtual ~Lights();
	protected:
		UniformBufferLight m_ubl{};
		VulkanMisc *vMisc;
		VmaBuffer m_vmaUniformBuffer;
		VmaBuffer m_vmaOffScreenShadowBuffer;
		Shadow* m_shadowData;
		LightData m_lightData;
		float m_nearPlane = 1.0f;
		float m_farPlane = 7.5f;
		int m_index = 0;
		bool m_shadow = false;
	};
}

#endif //__ENGINE_LIGHT__
