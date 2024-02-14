#ifndef __ENGINE_MODEL_MANAGER___
#define __ENGINE_MODEL_MANAGER___

#include "VulkanMisc.hpp"
#include "Debug.hpp"
#include <map>
#include "ShapeBuffer.hpp"
#include "Model.hpp"
#include <unordered_map>
#include "Manager.hpp"
#include "UniformBufferObject.hpp"
#include "BoundingBox.hpp"

namespace Ge
{
    class ModelManager final : public Manager
    {
    public:
        bool initiliaze(VulkanMisc *vM);
        void release();
        ShapeBuffer * allocateBuffer(const char *path, bool normal_recalculate = false);
        ShapeBuffer * allocateBufferWire(const char* path, bool normal_recalculate = false);
        ShapeBuffer * allocateBufferBone(const char* path, std::unordered_map<BoundingBox*, int> map, bool normal_recalculate = false);
        ShapeBuffer * allocateBufferRig(Model* currentmesh, ShapeBuffer* emptyShape, std::vector<BoneNode*> bones, Model** m);
        void createRoot(BoneNode* bn, ShapeBuffer * emptyShape);
        std::vector<ShapeBuffer*> allocateFBXBuffer(const char* path, bool normal_recalculate = false,std::vector<int> m_loadIdMesh = std::vector<int>());
        std::vector<ShapeBuffer*> allocateBuffers(const char* path, bool normal_recalculate = false);
		ShapeBuffer *allocateBuffer(float * pos, float * texCord, float * normal, unsigned int * indice, unsigned vertexSize, unsigned indiceSize);
		void printModelInfo(const char *path);
        Model *createModel(ShapeBuffer *buffer, std::string nom = "Model");
        void destroyModel(Model *model);
        void destroyBuffer(ShapeBuffer *buffer);
        void updateDescriptor();		        
		void initDescriptor(VulkanMisc * vM) override;
        static std::map<ShapeBuffer*, std::map<int, std::vector<Model*>>> & GetModelInstancing();
        static std::vector<Model*> & GetModels();
        static void updateInstanced(ShapeBuffer* sb, Model* m, int pi_Start, int pi_End, VulkanMisc* vm);
        void ComputationTangent(std::vector<Vertex>& vertices);
	private:
		friend class RenderingEngine;
		void destroyElement();
    private:		
        static std::map<ShapeBuffer*, std::map<int, std::vector<Model*>>> m_instancing;
        static std::vector<Model*> m_models;
        std::vector<ShapeBuffer *> m_shapeBuffers;        
		std::vector<Model *> m_destroymodels;
		std::vector<ShapeBuffer *> m_destroyshapeBuffers;
        VmaBuffer m_vmaUniformBuffers;
        VulkanMisc *vulkanM;
    };
}

#endif //__ENGINE_MODEL_MANAGER___