#ifndef __PREFAB_LOADER__
#define __PREFAB_LOADER__

#include "MetaData.hpp"
#include "Block.hpp"

class PrefabLoader : public ImguiBlock
{
public:
	PrefabLoader(MetaPrefab * mp,PipelinePrefab * pp,int sizem,int sizep);
	std::vector<LoadPrefab*> getPrefabs() const;

	void preRender(VulkanMisc* vM);
	void render(VulkanMisc* vM);
	void GenerateCode();

	~PrefabLoader();
private:
	ptrClass m_pc;
	std::vector<LoadPrefab*> m_prefab;
	std::map<std::string, Textures*> m_mapTexture;

	int m_id_Prefab = 0;
	bool m_collider = true;
	bool m_cam_position = true;
	bool m_generate_code = false;
	int m_listboxCurrentItem = 0;
	std::vector<const char*> cnames;
	std::vector<std::string> cnames_string;
	std::vector<Block*> m_temp_block;
	std::vector<int> m_temp_id;
	std::string m_start_code = "m_blocks.push_back(";
	std::string m_code = "";
	std::string m_end_code = ");\n";

	float m_scale_arrow = 0.2f;
	Model* m_arrow;
};

#endif //!__PREFAB_LOADER__