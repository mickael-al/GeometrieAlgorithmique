#ifndef __METADATA__
#define __METADATA__

#include "GameEngine.hpp"

struct ColliderTransform
{
	glm::vec3 collider_size;
	glm::vec3 collider_pos;
	glm::vec3 collider_eul;
};

struct MetaPrefab
{
	int load_mesh_type = 0;//obj = 0 && fbx == 1
	float mesh_resize_scale;
	std::vector<std::string> mesh;
	std::vector<std::vector<std::string>> texture;
	std::vector<ColliderTransform> collider_transform;
	float metallic;
	float roughness;
};

struct LoadPrefab
{
	std::vector<ShapeBuffer*> sb;
	std::vector<Materials*> mat;
	MetaPrefab meta;
};

struct PipelinePrefab
{
	std::vector<int> id;
	int typePipeline;
};

#endif // !__METADATA__