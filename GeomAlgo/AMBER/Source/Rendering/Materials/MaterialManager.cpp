#include "MaterialManager.hpp"
#include "ModelManager.hpp"
#include "tinyobjloader/tiny_obj_loader.h"
#include "Debug.hpp"

namespace Ge
{
	Materials * MaterialManager::defaultMaterial = nullptr;
    bool MaterialManager::initialize(VulkanMisc *vM)
    {
        vulkanM = vM;
		defaultMaterial = createMaterial();
        Debug::INITSUCCESS("MaterialManager");
        return true;
    }

	Materials * MaterialManager::getDefaultMaterial()
	{
		return defaultMaterial;
	}

    Materials * MaterialManager::createMaterial()
    {     
		Materials * material = new Materials(m_materials.size(), vulkanM);		
		m_materials.push_back(material);
		vulkanM->str_VulkanDescriptor->materialCount = m_materials.size();
        updateDescriptor();
		return material;
    }

    std::vector<Materials*> & MaterialManager::loadMltMaterial(const char* path,bool filter, TextureManager* tm)
    {
        std::vector<Materials*> mats;
        tinyobj::ObjReader reader;

        if (!reader.ParseFromFile(path))
        {
            if (!reader.Error().empty())
            {
                Debug::Error("TinyObjLoader error: %s", reader.Error().c_str());
            }
            return mats;
        }
        const std::vector<tinyobj::material_t>& materials = reader.GetMaterials();

        // Get the shapes from the reader
        const std::vector<tinyobj::shape_t>& shapes = reader.GetShapes();

        const char* lastSlash = strrchr(path, '/');
        char* result = nullptr;
        if (lastSlash != nullptr)
        {
            ptrdiff_t position = lastSlash - path;
            result = new char[position + 2];
            strncpy(result, path, position + 1);
            result[position + 1] = '\0';
        }
        if (result == nullptr)
        {
            return mats;
        }

        std::map<std::string, Textures*> m_textures;

        for (size_t shapeIdx = 0; shapeIdx < shapes.size(); shapeIdx++)
        {
            const tinyobj::shape_t& shape = shapes[shapeIdx];

            // Assuming each shape uses a single material
            if (shape.mesh.material_ids.size() > 0)
            {
                int materialIndex = shape.mesh.material_ids[0];
                if (materialIndex >= 0 && materialIndex < materials.size())
                {
                    const tinyobj::material_t& material = materials[materialIndex];

                    if (m_textures.find(material.diffuse_texname) == m_textures.end())
                    {
                        std::string npath = result + material.diffuse_texname;
                        m_textures[material.diffuse_texname] = tm->createTexture(npath.c_str(), filter);
                    }

                    mats.push_back(createMaterial());
                    mats.back()->setAlbedoTexture(m_textures[material.diffuse_texname]);
                }
            }
        }

        delete[] result;

        vulkanM->str_VulkanDescriptor->materialCount = m_materials.size();
        updateDescriptor();
        return mats;
    }

    void MaterialManager::destroyMaterial(Materials *material)
    {
		m_materials.erase(std::remove(m_materials.begin(), m_materials.end(), material), m_materials.end());
		int mat_pi = material->getPipelineIndex();
        //delete(material);
        m_destroyElement = true;
        m_destroymaterials.push_back(material);
		for (int i = 0; i < m_materials.size(); i++)
		{
			m_materials[i]->setIndex(i);
		}
		std::vector<Model *> all_models = ModelManager::GetModels();
		for (int i = 0; i < all_models.size();i++)
		{
			all_models[i]->majMaterialIndex(mat_pi);
		}
        vulkanM->str_VulkanDescriptor->materialCount = m_materials.size();
        updateDescriptor();
		vulkanM->str_VulkanDescriptor->recreateCommandBuffer = true;
    }

    void MaterialManager::destroyElement()
    {
        if (m_destroyElement)
        {
            for (int i = 0; i < m_destroymaterials.size(); i++)
            {
                delete m_destroymaterials[i];
            }
            m_destroymaterials.clear();
            m_destroyElement = false;
        }
    }

    void MaterialManager::updateDescriptor()
    {
        std::vector<VkDescriptorBufferInfo> bufferInfoMaterial{};
		VkDescriptorBufferInfo bufferIM{};
		for (int i = 0; i < m_materials.size(); i++)
		{
            bufferIM.buffer = m_materials[i]->getUniformBuffers();
			bufferIM.offset = 0;
			bufferIM.range = sizeof(UniformBufferMaterial);
			bufferInfoMaterial.push_back(bufferIM);
        }
        m_descriptor[0]->updateCount(vulkanM,m_materials.size(),bufferInfoMaterial);
    }

    void MaterialManager::release()
    {
        for (int i = 0; i < m_materials.size();i++)
		{
			delete (m_materials[i]);
		}
		m_materials.clear();
        for (int i = 0; i < m_destroymaterials.size(); i++)
        {
            delete m_destroymaterials[i];
        }
        m_destroymaterials.clear();
		for (int i = 0; i < m_descriptor.size(); i++)
		{
			delete m_descriptor[i];
		}
		m_descriptor.clear();
        Debug::RELEASESUCCESS("MaterialManager");
    }

	void MaterialManager::initDescriptor(VulkanMisc * vM)
	{
		if (m_descriptor.size() == 0)
		{
			m_descriptor.push_back(new Descriptor(vM, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1));
		}
	}

}