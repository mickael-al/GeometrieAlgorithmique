#ifndef __ENGINE_MATERIAL_MANAGER__
#define __ENGINE_MATERIAL_MANAGER__

#include "VulkanMisc.hpp"
#include "Manager.hpp"
#include "Materials.hpp"
#include "TextureManager.hpp"

namespace Ge
{
    class MaterialManager final: public Manager
    {
    public:
        bool initialize(VulkanMisc *vM);
        void release();
        Materials *createMaterial();
        std::vector<Materials*> & loadMltMaterial(const char* path, bool filter, TextureManager * tm);
        void destroyMaterial(Materials *material);
		void initDescriptor(VulkanMisc * vM);
		void updateDescriptor();    
		static Materials * getDefaultMaterial();
    private:
        friend class RenderingEngine;
        void destroyElement();
    private:
		static Materials * defaultMaterial;
        std::vector<Materials *> m_materials;
        std::vector<Materials *> m_destroymaterials;
        VulkanMisc *vulkanM;        
    };
}

#endif //__ENGINE_MATERIAL_MANAGER__