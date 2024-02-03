#include "ModelManager.hpp"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tinyobjloader/tiny_obj_loader.h"
#include <glm/gtx/normal.hpp>
#include "MaterialManager.hpp"
#include "OpenFBX/src/ofbx.h"
#include <algorithm>

namespace Ge
{
	std::vector<Model*> ModelManager::m_models;
	std::map<ShapeBuffer*, std::map<int, std::vector<Model*>>> ModelManager::m_instancing;
	bool ModelManager::initiliaze(VulkanMisc *vM)
	{
		vulkanM = vM;		
		if (!BufferManager::createBuffer(sizeof(UniformBufferObject), VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, m_vmaUniformBuffers, vM->str_VulkanDeviceMisc))
		{
			Debug::Error("Echec de la creation d'un uniform buffer");
			return false;
		}
		updateDescriptor();
		Debug::INITSUCCESS("ModelManager");
		return true;
	}

	void ModelManager::updateDescriptor()
	{
		std::vector<VkDescriptorBufferInfo> bufferInfoModel{};
		VkDescriptorBufferInfo bufferIM{};
        for (int i = 0 ; i < m_models.size();i++)
        {
            bufferIM.buffer = m_models[i]->getUniformBuffers();
			bufferIM.offset = 0;
			bufferIM.range = sizeof(UniformBufferObject);
			bufferInfoModel.push_back(bufferIM);
        }
		if(m_models.size() == 0)
		{
			bufferIM.buffer = m_vmaUniformBuffers.buffer;
			bufferIM.offset = 0;
			bufferIM.range = sizeof(UniformBufferObject);
			bufferInfoModel.push_back(bufferIM);
			m_descriptor[0]->updateCount(vulkanM, 1, bufferInfoModel);			
		}
		else
		{			
			m_descriptor[0]->updateCount(vulkanM, m_models.size(), bufferInfoModel);
		}
		vulkanM->str_VulkanSwapChainMisc->str_descriptorSetLayoutModel = m_descriptor[0]->getDescriptorSetLayout();
		vulkanM->str_VulkanSwapChainMisc->str_descriptorSetModel = m_descriptor[0]->getDescriptorSets();
	}

	void ModelManager::release()
	{
		for (int i = 0; i < m_shapeBuffers.size();i++)
		{
			delete (m_shapeBuffers[i]);
		}
		m_shapeBuffers.clear();
		for (int i = 0; i < m_models.size(); i++)
		{
			delete (m_models[i]);
		}
		m_models.clear();
		m_instancing.clear();
		BufferManager::destroyBuffer(m_vmaUniformBuffers);
		for (int i = 0; i < m_descriptor.size(); i++)
		{
			delete m_descriptor[i];
		}
		m_descriptor.clear();
		Debug::RELEASESUCCESS("ModelManager");
	}

	Model * ModelManager::createModel(ShapeBuffer *buffer, std::string nom)
	{
		if (buffer == nullptr)
		{
			Debug::Warn("Le buffer n'existe pas");
			return nullptr;
		}
		Model * Mesh = new Model(buffer, m_models.size(), vulkanM);
		m_instancing[buffer][0].push_back(Mesh);
		Mesh->setName(nom);
		Mesh->setMaterial(MaterialManager::getDefaultMaterial());
		m_models.push_back(Mesh);
		vulkanM->str_VulkanDescriptor->modelCount = m_models.size();		
		updateDescriptor();
		vulkanM->str_VulkanDescriptor->recreateCommandBuffer = true;
		vulkanM->str_VulkanDescriptor->recreateShadowPipeline = true;
		return Mesh;
	}

	void ModelManager::destroyModel(Model *model)
	{
		m_destroyElement = true;
		m_destroymodels.erase(std::remove(m_destroymodels.begin(), m_destroymodels.end(), model), m_destroymodels.end());
		m_destroymodels.push_back(model);
		vulkanM->str_VulkanDescriptor->recreateCommandBuffer = true;
		vulkanM->str_VulkanDescriptor->recreateShadowPipeline = true;
	}

	void ModelManager::destroyBuffer(ShapeBuffer *buffer)
	{
		m_destroyElement = true;
		m_destroyshapeBuffers.erase(std::remove(m_destroyshapeBuffers.begin(), m_destroyshapeBuffers.end(), buffer), m_destroyshapeBuffers.end());
		m_destroyshapeBuffers.push_back(buffer);
		vulkanM->str_VulkanDescriptor->recreateCommandBuffer = true;
		vulkanM->str_VulkanDescriptor->recreateShadowPipeline = true;
	}

	void ModelManager::destroyElement()
	{		
		if (m_destroyElement)
		{
			ShapeBuffer* sb;
			int pi;
			for (int j = 0; j < m_destroymodels.size(); j++)
			{
				m_models.erase(std::remove(m_models.begin(), m_models.end(), m_destroymodels[j]), m_models.end());
				sb = m_destroymodels[j]->getShapeBuffer();
				pi = m_destroymodels[j]->getMaterial()->getPipelineIndex();
				m_instancing[sb][pi].erase(std::remove(m_instancing[sb][pi].begin(), m_instancing[sb][pi].end(), m_destroymodels[j]), m_instancing[sb][pi].end());
				auto& modelVector = m_instancing[sb][pi];
				if (modelVector.empty())
				{
					m_instancing[sb].erase(pi);
					if (m_instancing[sb].empty())
					{
						m_instancing.erase(sb);
					}
				}
				delete (m_destroymodels[j]);
			}

			for (int j = 0; j < m_destroyshapeBuffers.size(); j++)
			{
				for (int i = 0; i < m_models.size(); i++)
				{
					if (m_models[i]->getShapeBuffer() == m_destroyshapeBuffers[j])
					{
						Model * m = m_models[i];
						m_models.erase(std::remove(m_models.begin(), m_models.end(), m), m_models.end());
						delete (m);
						i--;
					}
				}
				auto it = m_instancing.find(m_destroyshapeBuffers[j]);
				if (it != m_instancing.end()) 
				{
					m_instancing.erase(it);
				}
				m_shapeBuffers.erase(std::remove(m_shapeBuffers.begin(), m_shapeBuffers.end(), m_destroyshapeBuffers[j]), m_shapeBuffers.end());
				delete (m_destroyshapeBuffers[j]);
			}

			m_destroyshapeBuffers.clear();
			m_destroymodels.clear();
			for (int i = 0; i < m_models.size(); i++)
			{
				m_models[i]->setIndexUbo(i);
			}
			vulkanM->str_VulkanDescriptor->modelCount = m_models.size();
			updateDescriptor();
			m_destroyElement = false;
		}
	}

	std::map<ShapeBuffer*, std::map<int, std::vector<Model*>>> & ModelManager::GetModelInstancing()
	{
		return m_instancing;
	}

	std::vector<Model*> & ModelManager::GetModels()
	{
		return m_models;
	}

	void ModelManager::updateInstanced(ShapeBuffer* sb, Model* m,int pi_Start,int pi_End,VulkanMisc * vm)
	{
		if (pi_Start == pi_End)
		{
			return;
		}
		auto& modelVector = m_instancing[sb][pi_Start];
		modelVector.erase(std::remove(modelVector.begin(), modelVector.end(), m), modelVector.end());

		m_instancing[sb][pi_End].push_back(m);

		if (modelVector.empty())
		{
			m_instancing[sb].erase(pi_Start);
			if (m_instancing[sb].empty())
			{
				m_instancing.erase(sb);
			}
		}
		vm->str_VulkanDescriptor->recreateCommandBuffer = true;
	}

	void ModelManager::initDescriptor(VulkanMisc * vM)
	{
		if (m_descriptor.size() == 0)
		{
			m_descriptor.push_back(new Descriptor(vM, VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1));
			vM->str_VulkanSwapChainMisc->str_descriptorSetLayoutModel = m_descriptor[0]->getDescriptorSetLayout();
			vM->str_VulkanSwapChainMisc->str_descriptorSetModel = m_descriptor[0]->getDescriptorSets();
		}
	}

	void ModelManager::printModelInfo(const char *path)
	{
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string warn, err;
		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		glm::vec3 normalResult;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path))
		{
			Debug::Warn("%s  %s", nullptr, warn.c_str(), err.c_str());
		}

		std::unordered_map<Vertex, uint32_t> uniqueVertices{};

		for (const auto &shape : shapes)
		{
			Debug::Log("%s", path);
			for (const auto &index : shape.mesh.indices)
			{
				Vertex vertex{};

				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2] };

				vertex.texCoord = {
					attrib.texcoords[2 * index.texcoord_index + 0],
					1.0f - attrib.texcoords[2 * index.texcoord_index + 1] };

				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2] };

				if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.push_back(vertex);
				}

				indices.push_back(uniqueVertices[vertex]);
			}
		}
		Debug::Info("POS: %d", vertices.size());
		for (int i = 0; i < vertices.size(); i++)
		{
			std::cout << vertices[i].pos.x << "," << vertices[i].pos.y << "," << vertices[i].pos.z << ",";
		}
		std::cout << std::endl;
		Debug::Info("TEXCORD: %d", vertices.size());
		for (int i = 0; i < vertices.size(); i++)
		{
			std::cout << vertices[i].texCoord.x << "," << vertices[i].texCoord.y << ",";
		}
		std::cout << std::endl;
		Debug::Info("NORMAL: %d", vertices.size());
		for (int i = 0; i < vertices.size(); i++)
		{
			std::cout << vertices[i].normal.x << "," << vertices[i].normal.y << "," << vertices[i].normal.z << ",";
		}
		std::cout << std::endl;
		Debug::Info("INDICE: %d", indices.size());
		for (int i = 0; i < indices.size(); i++)
		{
			std::cout << indices[i] << ",";
		}
		std::cout << std::endl;
	}

	ShapeBuffer * ModelManager::allocateBuffer(float * pos, float * texCord, float * normal, unsigned int * indice, unsigned vertexSize, unsigned indiceSize)
	{
		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices(indice, indice + indiceSize);

		vertices.reserve(vertexSize/3);
		indices.reserve(vertexSize);
		for (int i = 0; i < vertexSize; i++)
		{
			Vertex vertex{};

			vertex.pos = {
				pos[3 * i + 0],
				pos[3 * i + 1],
				pos[3 * i + 2] };

			vertex.texCoord = {
				texCord[2 * i + 0],
				texCord[2 * i + 1] };

			vertex.normal = {
				normal[3 * i + 0],
				normal[3 * i + 1],
				normal[3 * i + 2] };
			vertices.push_back(vertex);
		}

		uint32_t index0, index1, index2;
		glm::vec3 edge1, edge2, tangent;
		glm::vec2 deltaUV1, deltaUV2;
		float r;
		for (size_t i = 0; i < indices.size(); i += 3)
		{
			index0 = indices[i];
			index1 = indices[i + 1];
			index2 = indices[i + 2];

			Vertex& vertex0 = vertices[index0];
			Vertex& vertex1 = vertices[index1];
			Vertex& vertex2 = vertices[index2];

			edge1 = vertex1.pos - vertex0.pos;
			edge2 = vertex2.pos - vertex0.pos;

			deltaUV1 = vertex1.texCoord - vertex0.texCoord;
			deltaUV2 = vertex2.texCoord - vertex0.texCoord;

			r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

			vertex0.tangents = tangent;
			vertex1.tangents = tangent;
			vertex2.tangents = tangent;
		}

		ShapeBuffer *buffer = new ShapeBuffer(vertices, indices, vulkanM);
		m_shapeBuffers.push_back(buffer);
		return buffer;
	}

	void ModelManager::ComputationTangent(std::vector<Vertex> &vertices)
	{
		glm::vec3 edge1;
		glm::vec3 edge2;
		glm::vec2 deltaUV1;
		glm::vec2 deltaUV2;
		glm::vec3 tangents;
		
		float r;
		for (int i = 0; i+2 < vertices.size(); i+=3)
		{
			edge1 = vertices[i + 1].pos - vertices[i].pos;
			edge2 = vertices[i + 2].pos - vertices[i].pos;

			deltaUV1 = vertices[i + 1].texCoord - vertices[i].texCoord;
			deltaUV2 = vertices[i + 2].texCoord - vertices[i].texCoord;

			r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangents.x = r * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
			tangents.y = r * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
			tangents.z = r * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);

			tangents = glm::normalize(tangents);

			vertices[i + 0].tangents = tangents;
			vertices[i + 1].tangents = tangents;
			vertices[i + 2].tangents = tangents;
		}
	}

	std::vector<ShapeBuffer*> ModelManager::allocateBuffers(const char* path,bool normal_recalculate)
	{
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::vector<ShapeBuffer*> sbvec;
		std::string warn, err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path))
		{
			Debug::Warn("%s  %s", warn.c_str(), err.c_str());
			return sbvec;
		}

		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		std::unordered_map<Vertex, uint32_t> uniqueVertices;

		vertices.reserve(attrib.vertices.size() / 3);
		indices.reserve(attrib.vertices.size());

		for (const auto& shape : shapes)
		{
			vertices.clear();
			indices.clear();
			uniqueVertices.clear();
			for (const auto& index : shape.mesh.indices)
			{
				Vertex vertex{};

				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.texCoord = {
					attrib.texcoords[2 * index.texcoord_index + 0],
					1.0f - attrib.texcoords[2 * index.texcoord_index + 1]
				};

				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2]
				};

				vertex.color = { 1, 1, 1 };
				vertex.tangents = { 0, 0, 0 };

				if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.emplace_back(vertex);
				}

				indices.emplace_back(uniqueVertices[vertex]);
			}

			uint32_t index0, index1, index2;
			glm::vec3 edge1, edge2, tangent;
			glm::vec2 deltaUV1, deltaUV2;
			float r;
			for (size_t i = 0; i < indices.size(); i += 3)
			{
				index0 = indices[i];
				index1 = indices[i + 1];
				index2 = indices[i + 2];

				Vertex& vertex0 = vertices[index0];
				Vertex& vertex1 = vertices[index1];
				Vertex& vertex2 = vertices[index2];

				edge1 = vertex1.pos - vertex0.pos;
				edge2 = vertex2.pos - vertex0.pos;

				deltaUV1 = vertex1.texCoord - vertex0.texCoord;
				deltaUV2 = vertex2.texCoord - vertex0.texCoord;

				r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

				tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

				vertex0.tangents = tangent;
				vertex1.tangents = tangent;
				vertex2.tangents = tangent;
				if (normal_recalculate)
				{
					glm::vec3 normal = glm::normalize(glm::cross(vertex1.pos- vertex0.pos, vertex2.pos - vertex0.pos));
					vertex0.normal = normal;
					vertex1.normal = normal;
					vertex2.normal = normal;
				}
			}
			ShapeBuffer* buffer = new ShapeBuffer(vertices, indices, vulkanM);
			m_shapeBuffers.push_back(buffer);
			sbvec.push_back(buffer);
		}

		return sbvec;
	}

	std::vector<ShapeBuffer*> ModelManager::allocateFBXBuffer(const char* path, bool normal_recalculate, std::vector<int> m_loadIdMesh)
	{
		std::vector<ShapeBuffer*> sbvec;
		FILE* fp = fopen(path, "rb");

		if (!fp)
		{
			Debug::Warn("cannot load : %s ", path);
			return sbvec;
		}
		fseek(fp, 0, SEEK_END);
		long file_size = ftell(fp);
		fseek(fp, 0, SEEK_SET);
		auto* content = new ofbx::u8[file_size];
		fread(content, 1, file_size, fp);

		ofbx::LoadFlags flags =
			ofbx::LoadFlags::TRIANGULATE |
			//		ofbx::LoadFlags::IGNORE_MODELS |
			ofbx::LoadFlags::IGNORE_BLEND_SHAPES |
			ofbx::LoadFlags::IGNORE_CAMERAS |
			ofbx::LoadFlags::IGNORE_LIGHTS |
			//		ofbx::LoadFlags::IGNORE_TEXTURES |
			ofbx::LoadFlags::IGNORE_SKIN |
			ofbx::LoadFlags::IGNORE_BONES |
			ofbx::LoadFlags::IGNORE_PIVOTS |
			//		ofbx::LoadFlags::IGNORE_MATERIALS |
			ofbx::LoadFlags::IGNORE_POSES |
			ofbx::LoadFlags::IGNORE_VIDEOS |
			ofbx::LoadFlags::IGNORE_LIMBS |
			//		ofbx::LoadFlags::IGNORE_MESHES |
			ofbx::LoadFlags::IGNORE_ANIMATIONS;

		ofbx::IScene* scene = ofbx::load((ofbx::u8*)content, file_size, (ofbx::u16)flags);

		int size = (m_loadIdMesh.size() > 0) ? m_loadIdMesh.size() : scene->getMeshCount();

		for (int i = 0; i < size; ++i)
		{
			std::vector<Vertex> vertices;
			std::vector<uint32_t> indices;
			std::unordered_map<Vertex, uint32_t> uniqueVertices;

			const ofbx::Geometry* geom = nullptr;
			if (m_loadIdMesh.size() > 0)
			{
				geom = scene->getMesh(m_loadIdMesh[i]% scene->getMeshCount())->getGeometry();
			}
			else
			{
				geom = scene->getMesh(i)->getGeometry();
			}			

			vertices.reserve(geom->getVertexCount());
			indices.reserve(geom->getIndexCount());
			const ofbx::Vec3* position = geom->getVertices();
			const ofbx::Vec2* uv = geom->getUVs();
			const ofbx::Vec3* normals = geom->getNormals();
			const ofbx::Vec3* tangents = geom->getTangents();
			const ofbx::Vec4* colors = geom->getColors();

			for (int j = 0; j < geom->getVertexCount(); ++j)
			{								
				Vertex vertex{};
				vertex.pos.x = position[j].x;
				vertex.pos.y = position[j].y;
				vertex.pos.z = position[j].z;

				vertex.texCoord.x = uv[j].x;
				vertex.texCoord.y = 1.0f - uv[j].y;
				
				if (normals != nullptr && !normal_recalculate)
				{
					vertex.normal.x = normals[j].x;
					vertex.normal.y = normals[j].y;
					vertex.normal.z = normals[j].z;
				}

				if (tangents != nullptr)
				{
					vertex.tangents.x = tangents[j].x;
					vertex.tangents.y = tangents[j].y;
					vertex.tangents.z = tangents[j].z;
				}
				else
				{					
					vertex.tangents = { 0, 0, 0 };
				}

				if (colors != nullptr)
				{
					vertex.color.x = colors[j].x;
					vertex.color.y = colors[j].y;
					vertex.color.z = colors[j].z;
				}
				else
				{
					vertex.color = { 1, 1, 1 };
				}
				
				if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.emplace_back(vertex);
				}

				indices.emplace_back(uniqueVertices[vertex]);
			}

			if (tangents == nullptr || normal_recalculate || normals == nullptr)
			{
				uint32_t index0, index1, index2;
				glm::vec3 edge1, edge2, tangent;
				glm::vec2 deltaUV1, deltaUV2;
				float r;
				for (size_t i = 0; i < indices.size(); i += 3)
				{
					index0 = indices[i];
					index1 = indices[i + 1];
					index2 = indices[i + 2];

					Vertex& vertex0 = vertices[index0];
					Vertex& vertex1 = vertices[index1];
					Vertex& vertex2 = vertices[index2];

					if (tangents == nullptr)
					{
						edge1 = vertex1.pos - vertex0.pos;
						edge2 = vertex2.pos - vertex0.pos;

						deltaUV1 = vertex1.texCoord - vertex0.texCoord;
						deltaUV2 = vertex2.texCoord - vertex0.texCoord;

						r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

						tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

						vertex0.tangents = tangent;
						vertex1.tangents = tangent;
						vertex2.tangents = tangent;
					}
					if (normal_recalculate || normals == nullptr)
					{
						glm::vec3 normal = glm::normalize(glm::cross(vertex1.pos - vertex0.pos, vertex2.pos - vertex0.pos));
						vertex0.normal = normal;
						vertex1.normal = normal;
						vertex2.normal = normal;
					}
				}
			}
			ShapeBuffer* buffer = new ShapeBuffer(vertices, indices, vulkanM);
			m_shapeBuffers.push_back(buffer);
			sbvec.push_back(buffer);
		}

		delete[] content;
		fclose(fp);

		return sbvec;
	}
	/*glm::vec3 barycentric(glm::vec2 v0, glm::vec2 v1, glm::vec2 v2, glm::vec2 p)
	{
		// On calcule les vecteurs du triangle
		glm::vec2 e0 = v1 - v0;
		glm::vec2 e1 = v2 - v0;
		glm::vec2 e2 = p - v0;

		// On calcule les déterminants
		float d00 = glm::dot(e0, e0);
		float d01 = glm::dot(e0, e1);
		float d11 = glm::dot(e1, e1);
		float d20 = glm::dot(e2, e0);
		float d21 = glm::dot(e2, e1);
		float invDenom = 1.0 / (d00 * d11 - d01 * d01);

		// On calcule les coordonnées barycentriques
		float v = (d11 * d20 - d01 * d21) * invDenom;
		float w = (d00 * d21 - d01 * d20) * invDenom;
		float u = 1.0 - v - w;

		// On renvoie le résultat sous forme de vec3
		return glm::vec3(u, v, w);
	}*/

	glm::vec3 barycentric(glm::vec2 A, glm::vec2 B, glm::vec2 C, glm::vec2 P)
	{
		// Calcul des côtés du triangle
		glm::vec2 v0 = B - A;
		glm::vec2 v1 = C - A;
		glm::vec2 v2 = P - A;

		// Calcul du déterminant principal
		float det = v0.x * v1.y - v1.x * v0.y;

		// Calcul des coordonnées barycentriques
		float u = (v1.y * v2.x - v1.x * v2.y) / det;
		float v = (v0.x * v2.y - v0.y * v2.x) / det;
		float w = 1.0 - u - v;

		// Retourne les coordonnées barycentriques
		return glm::vec3(u, v, w);
	}

	bool isPointInsideBoundingBox(const glm::vec3& point, const glm::vec3& min, const glm::vec3& max) 
	{
		return (point.x >= min.x && point.x <= max.x &&
			point.y >= min.y && point.y <= max.y &&
			point.z >= min.z && point.z <= max.z);
	}

	ShapeBuffer * ModelManager::allocateBufferBone(const char* path, std::unordered_map<BoundingBox*, int> m_map, bool normal_recalculate)
	{
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string warn, err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path))
		{
			Debug::Warn("%s  %s", warn.c_str(), err.c_str());
			return nullptr;
		}

		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		std::unordered_map<Vertex, uint32_t> uniqueVertices;

		vertices.reserve(attrib.vertices.size() / 3);
		indices.reserve(attrib.vertices.size());

		for (const auto& shape : shapes)
		{
			for (const auto& index : shape.mesh.indices)
			{
				Vertex vertex{};

				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.texCoord = {
					attrib.texcoords[2 * index.texcoord_index + 0],
					1.0f - attrib.texcoords[2 * index.texcoord_index + 1]
				};

				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2]
				};

				vertex.color = { 1, 1, 1 };
				vertex.tangents = { 0, 0, 0 };

				if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.emplace_back(vertex);
					for (const auto& bb : m_map)
					{
						glm::vec3 min_pos = bb.first->box->getPosition() - bb.first->box->getScale() / 2.0f;
						glm::vec3 max_pos = bb.first->box->getPosition() + bb.first->box->getScale() / 2.0f;
						if (isPointInsideBoundingBox(vertex.pos, min_pos, max_pos))
						{
							vertex.color.x = (int)bb.second;
						}
					}
				}

				indices.emplace_back(uniqueVertices[vertex]);
			}
		}

		uint32_t index0, index1, index2;
		glm::vec3 edge1, edge2, tangent;
		glm::vec2 deltaUV1, deltaUV2;
		float r;
		for (size_t i = 0; i < indices.size(); i += 3)
		{
			index0 = indices[i];
			index1 = indices[i + 1];
			index2 = indices[i + 2];

			Vertex& vertex0 = vertices[index0];
			Vertex& vertex1 = vertices[index1];
			Vertex& vertex2 = vertices[index2];			

			edge1 = vertex1.pos - vertex0.pos;
			edge2 = vertex2.pos - vertex0.pos;

			deltaUV1 = vertex1.texCoord - vertex0.texCoord;
			deltaUV2 = vertex2.texCoord - vertex0.texCoord;

			r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

			vertex0.tangents = tangent;
			vertex1.tangents = tangent;
			vertex2.tangents = tangent;
			if (normal_recalculate)
			{
				glm::vec3 normal = glm::normalize(glm::cross(vertex1.pos - vertex0.pos, vertex2.pos - vertex0.pos));
				vertex0.normal = normal;
				vertex1.normal = normal;
				vertex2.normal = normal;
			}
		}

		ShapeBuffer* buffer = new ShapeBuffer(vertices, indices, vulkanM);
		m_shapeBuffers.push_back(buffer);
		return buffer;
	}



	ShapeBuffer * ModelManager::allocateBufferWire(const char* path, bool normal_recalculate)
	{
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string warn, err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path))
		{
			Debug::Warn("%s  %s", warn.c_str(), err.c_str());
			return nullptr;
		}

		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		std::unordered_map<Vertex, uint32_t> uniqueVertices;

		vertices.reserve(attrib.vertices.size() / 3);
		indices.reserve(attrib.vertices.size());

		for (const auto& shape : shapes)
		{
			for (const auto& index : shape.mesh.indices)
			{
				Vertex vertex{};

				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.texCoord = {
					attrib.texcoords[2 * index.texcoord_index + 0],
					1.0f - attrib.texcoords[2 * index.texcoord_index + 1]
				};

				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2]
				};

				vertex.color = { 1, 1, 1 };
				vertex.tangents = { 0, 0, 0 };

				//if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.emplace_back(vertex);
				}

				indices.emplace_back(uniqueVertices[vertex]);
			}
		}

		uint32_t index0, index1, index2;
		glm::vec3 edge1, edge2, tangent;
		glm::vec2 deltaUV1, deltaUV2;
		float r;
		for (size_t i = 0; i < indices.size(); i += 3)
		{
			index0 = indices[i];
			index1 = indices[i + 1];
			index2 = indices[i + 2];

			Vertex& vertex0 = vertices[index0];
			Vertex& vertex1 = vertices[index1];
			Vertex& vertex2 = vertices[index2];

			vertex0.color = glm::vec3(1,0,0);
			vertex1.color = glm::vec3(0,1,0);
			vertex2.color = glm::vec3(0,0,1);

			edge1 = vertex1.pos - vertex0.pos;
			edge2 = vertex2.pos - vertex0.pos;

			deltaUV1 = vertex1.texCoord - vertex0.texCoord;
			deltaUV2 = vertex2.texCoord - vertex0.texCoord;

			r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

			vertex0.tangents = tangent;
			vertex1.tangents = tangent;
			vertex2.tangents = tangent;
			if (normal_recalculate)
			{
				glm::vec3 normal = glm::normalize(glm::cross(vertex1.pos - vertex0.pos, vertex2.pos - vertex0.pos));
				vertex0.normal = normal;
				vertex1.normal = normal;
				vertex2.normal = normal;
			}
		}

		ShapeBuffer* buffer = new ShapeBuffer(vertices, indices, vulkanM);
		m_shapeBuffers.push_back(buffer);
		return buffer;
	}

	ShapeBuffer* ModelManager::allocateBuffer(const char* path, bool normal_recalculate)
	{
		tinyobj::attrib_t attrib;
		std::vector<tinyobj::shape_t> shapes;
		std::vector<tinyobj::material_t> materials;
		std::string warn, err;

		if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, path))
		{
			Debug::Warn("%s  %s", warn.c_str(), err.c_str());
			return nullptr;
		}

		std::vector<Vertex> vertices;
		std::vector<uint32_t> indices;
		std::unordered_map<Vertex, uint32_t> uniqueVertices;

		vertices.reserve(attrib.vertices.size() / 3);
		indices.reserve(attrib.vertices.size());

		for (const auto& shape : shapes)
		{
			for (const auto& index : shape.mesh.indices)
			{
				Vertex vertex{};

				vertex.pos = {
					attrib.vertices[3 * index.vertex_index + 0],
					attrib.vertices[3 * index.vertex_index + 1],
					attrib.vertices[3 * index.vertex_index + 2]
				};

				vertex.texCoord = {
					attrib.texcoords[2 * index.texcoord_index + 0],
					1.0f - attrib.texcoords[2 * index.texcoord_index + 1]
				};

				vertex.normal = {
					attrib.normals[3 * index.normal_index + 0],
					attrib.normals[3 * index.normal_index + 1],
					attrib.normals[3 * index.normal_index + 2]
				};

				vertex.color = { 1, 1, 1 };
				vertex.tangents = { 0, 0, 0 };

				if (uniqueVertices.count(vertex) == 0)
				{
					uniqueVertices[vertex] = static_cast<uint32_t>(vertices.size());
					vertices.emplace_back(vertex);
				}

				indices.emplace_back(uniqueVertices[vertex]);
			}
		}

		uint32_t index0, index1, index2;
		glm::vec3 edge1, edge2, tangent;
		glm::vec2 deltaUV1, deltaUV2;
		float r;
		for (size_t i = 0; i < indices.size(); i += 3)
		{
			index0 = indices[i];
			index1 = indices[i + 1];
			index2 = indices[i + 2];

			Vertex& vertex0 = vertices[index0];
			Vertex& vertex1 = vertices[index1];
			Vertex& vertex2 = vertices[index2];

			edge1 = vertex1.pos - vertex0.pos;
			edge2 = vertex2.pos - vertex0.pos;

			deltaUV1 = vertex1.texCoord - vertex0.texCoord;
			deltaUV2 = vertex2.texCoord - vertex0.texCoord;

			 r = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

			tangent = glm::normalize(r * (deltaUV2.y * edge1 - deltaUV1.y * edge2));

			vertex0.tangents = tangent;
			vertex1.tangents = tangent;
			vertex2.tangents = tangent;
			if (normal_recalculate)
			{
				glm::vec3 normal = glm::normalize(glm::cross(vertex1.pos - vertex0.pos, vertex2.pos - vertex0.pos));
				vertex0.normal = normal;
				vertex1.normal = normal;
				vertex2.normal = normal;
			}
		}

		ShapeBuffer* buffer = new ShapeBuffer(vertices, indices, vulkanM);
		m_shapeBuffers.push_back(buffer);
		return buffer;
	}
}