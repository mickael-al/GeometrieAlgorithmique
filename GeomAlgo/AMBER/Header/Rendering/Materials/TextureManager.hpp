#ifndef __ENGINE_TEXTURE_MANAGER__
#define __ENGINE_TEXTURE_MANAGER__

#include "VulkanMisc.hpp"
#include "Textures.hpp"
#include "TextureCubeMap.hpp"
#include "Manager.hpp"

namespace Ge
{
    class TextureManager final : public Manager
    {
    private:
        friend class RenderingEngine;
        bool initialize(VulkanMisc *vM);
        void release();                
    public:
        Textures * createTexture(const char * path,bool filter = true);
        Textures* createTextureImGUI(const char* path, bool filter = true);
        Textures * createTexture(const std::vector<unsigned char> & imageData, bool filter = true);
        Textures * createTexture(unsigned char* pixel, int tw, int th, bool filter = true);
		TextureCubeMap * createTextureCubeMap(const char * path, bool filter = true);
		void destroyTexture(Textures * texture);
        void destroyImGUITexture(Textures* texture);
		void destroyTextureCubeMap(TextureCubeMap * texture);
		void initDescriptor(VulkanMisc * vM);
		void updateDescriptor();
		static std::vector<unsigned char *> convertCubMap(unsigned char * pixel, int tw, int th);
		static TextureCubeMap * GetNullCubeMap();
    private:
        VulkanMisc *vulkanM;
        std::vector<Textures *> m_textures;
        std::vector<Textures*> m_ImGUItextures;
		std::vector<TextureCubeMap *> m_texturesCube;
		Textures * nullTexture;
        Textures* normalTexture;
		static TextureCubeMap * s_nullTextureCubeMap;
    };
}

#endif //__ENGINE_TEXTURE_MANAGER__