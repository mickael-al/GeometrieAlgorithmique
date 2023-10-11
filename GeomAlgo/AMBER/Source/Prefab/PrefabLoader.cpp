#include "PrefabLoader.hpp"

PrefabLoader::PrefabLoader(MetaPrefab * mp, PipelinePrefab* pp, int sizem, int sizep)
{
	m_pc = GameEngine::getPtrClass();
	int numPrefabs = sizem;
	if (numPrefabs == 0)
	{
		return;
	}
	std::map<std::string, std::vector<ShapeBuffer*>> mapShape;
	for (int i = 0; i < numPrefabs; i++)
	{
		LoadPrefab* prefab = new LoadPrefab();
		prefab->sb.reserve(mp[i].mesh.size());
		prefab->mat.reserve(mp[i].texture.size());
		for (int j = 0; j < mp[i].mesh.size(); j++)
		{
			if (mapShape.find(mp[i].mesh[j]) != mapShape.end())
			{
				prefab->sb.insert(prefab->sb.end(), mapShape[mp[i].mesh[j]].begin(), mapShape[mp[i].mesh[j]].end());
			}
			else
			{
				if (mp[i].load_mesh_type == 0)
				{
					prefab->sb.push_back(m_pc.modelManager->allocateBuffer(mp[i].mesh[j].c_str()));
					mapShape[mp[i].mesh[j]] = prefab->sb;
				}
				else if (mp[i].load_mesh_type == 1)
				{
					std::vector<ShapeBuffer*> sbl = m_pc.modelManager->allocateFBXBuffer(mp[i].mesh[j].c_str());
					prefab->sb.insert(prefab->sb.end(), sbl.begin(), sbl.end());
					mapShape[mp[i].mesh[j]] = prefab->sb;
				}
			}			
		}
		for (int s = 0; s < mp[i].texture.size(); s++)
		{
			prefab->mat.push_back(m_pc.materialManager->createMaterial());
			prefab->mat[s]->setMetallic(mp[i].metallic);
			prefab->mat[s]->setRoughness(mp[i].roughness);
			bool findPipeline = false;
			for (int z = 0; z < sizep && !findPipeline; z++)
			{				
				for (int h = 0; h < pp[z].id.size() && !findPipeline; h++)
				{
					if (pp[z].id[h] == i)
					{
						prefab->mat[s]->setPipelineIndex(pp[z].typePipeline);
						findPipeline = true;
					}
				}
			}
			for (int t = 0; t < mp[i].texture[s].size(); t++)
			{
				if (!mp[i].texture[s][t].empty())
				{
					Textures* texture = nullptr;
					if (m_mapTexture.find(mp[i].texture[s][t]) != m_mapTexture.end())
					{
						texture = m_mapTexture[mp[i].texture[s][t]];
					}
					else
					{
						texture = m_pc.textureManager->createTexture(mp[i].texture[s][t].c_str(), false);
					}
					switch (t)
					{
					case 0:
						prefab->mat[s]->setAlbedoTexture(texture);
						break;
					case 1:
						prefab->mat[s]->setRoughnessTexture(texture);
						break;
					case 2:
						prefab->mat[s]->setMetallicTexture(texture);
						break;
					case 3:
						prefab->mat[s]->setNormalTexture(texture);
						break;
					default:
						break;
					}
					m_mapTexture[mp[i].texture[s][t]] = texture;
				}
			}
		}
		prefab->meta = mp[i];
		m_prefab.push_back(prefab);
	}
}

std::vector<LoadPrefab*> PrefabLoader::getPrefabs() const
{
	return m_prefab;
}

void PrefabLoader::preRender(VulkanMisc* vM)
{
	ShapeBuffer * sb = m_pc.modelManager->allocateBuffer("../Model/arrow.obj");
	m_arrow = m_pc.modelManager->createModel(sb);
	Materials* mat = m_pc.materialManager->createMaterial();
	mat->setAlbedoTexture(m_pc.textureManager->createTexture("../Asset/Debug/arrow.png",false));
	m_arrow->setMaterial(mat);	
	m_arrow->setScale(glm::vec3(m_scale_arrow));
}

void PrefabLoader::GenerateCode()
{
	m_code = "";
	for (int i = 0; i < m_temp_block.size(); i++)
	{
		glm::vec3 pos = m_temp_block[i]->getMainModel()->getPosition();
		glm::vec3 euler = m_temp_block[i]->getMainModel()->getEulerAngles();
		glm::vec3 scale = m_temp_block[i]->getMainModel()->getScale() / m_prefab[m_temp_id[i]]->meta.mesh_resize_scale;
		m_code += m_start_code + "new Block(m_pc, lp[" + std::to_string(m_temp_id[i]) + "], glm::vec3(" + std::to_string(pos.x) + "," + std::to_string(pos.y) + "," + std::to_string(pos.z) + ")" +
			", glm::vec3("
			+
			((euler.x == euler.y && euler.x == euler.z)
				?
				std::to_string(euler.x) + ")"

				:
				std::to_string(euler.x) + "," + std::to_string(euler.y) + "," + std::to_string(euler.z) + ")"
				)
			+

			", glm::vec3("


			+
			((scale.x == scale.y && scale.x == scale.z)
				?
				std::to_string(scale.x) + ")"
				:
				std::to_string(scale.x) + "," + std::to_string(scale.y) + "," + std::to_string(scale.z) + ")"
				)
			+
				(m_temp_block[i]->getCollision().size() > 0 ? "" : ",false")
			+
			")" + m_end_code;

	}
	Debug::TempFile(m_code);
}

void PrefabLoader::render(VulkanMisc* vM)
{
	ImGuiWindowFlags window_flags;
	window_flags = 0;
	ImGui::SetNextWindowBgAlpha(0.55f);
	if (ImGui::Begin("Prefab Loader", nullptr, window_flags))
	{
		//ImGui::SetWindowSize("Prefab Loader", ImVec2(vM->str_VulkanSwapChainMisc->str_swapChainExtent.width / 4.2f, vM->str_VulkanSwapChainMisc->str_swapChainExtent.height / 1.4f));
		ImGui::TextColored(ImVec4(0.2f, 1, 0.2f, 1), "Prefabs\n\n");
		ImGui::Text("Prefabs count : %d", m_prefab.size());

		if (ImGui::Button("-"))
		{
			m_id_Prefab = std::max(m_id_Prefab - 1, 0);
		}

		ImGui::SameLine();

		if (ImGui::Button("+"))
		{
			m_id_Prefab = std::min(m_id_Prefab + 1, static_cast<int>(m_prefab.size() - 1));
		}
		if (ImGui::SliderFloat("Arrow Scale", &m_scale_arrow, 0, 2))
		{
			m_arrow->setScale(glm::vec3(m_scale_arrow));
		}

		if (ImGui::DragInt("ID", &m_id_Prefab, 1, 0, m_prefab.size()))
		{
			m_id_Prefab = std::min(m_id_Prefab + 1, static_cast<int>(m_prefab.size() - 1));
		}
		ImGui::Text(m_prefab[m_id_Prefab]->meta.mesh[0].c_str());
		ImGui::Checkbox("Collider", &m_collider);
		ImGui::Checkbox("Camera spawn position", &m_cam_position);
		ImGui::Checkbox("Generate Code Menu", &m_generate_code);

		if (m_generate_code)
		{
			ImGui::InputText("Start Code ", const_cast<char*>(m_start_code.c_str()), m_start_code.size());
			ImGui::InputText("End Code ", const_cast<char*>(m_end_code.c_str()), m_end_code.size());
			if (ImGui::Button("Generate"))
			{
				GenerateCode();
			}

			ImGui::InputTextMultiline("Generated Code", const_cast<char*>(m_code.c_str()), m_code.size(), ImVec2(-1, ImGui::GetTextLineHeight() * 50));
		}
		else
		{

			if (ImGui::Button("Create Block"))
			{
				GenerateCode();
				m_temp_block.push_back(new Block(m_pc, m_prefab[m_id_Prefab], m_cam_position ? m_pc.cameraManager->getCurrentCamera()->getPosition() : glm::vec3(0.0f), glm::vec3(0), glm::vec3(1), m_collider));				
				m_temp_id.push_back(m_id_Prefab);
				cnames_string.push_back(m_prefab[m_id_Prefab]->meta.mesh[0] + " " + std::to_string(cnames_string.size()));
				cnames.push_back(cnames_string[cnames_string.size() - 1].c_str());
				m_listboxCurrentItem = std::max(m_listboxCurrentItem, static_cast<int>(m_temp_block.size()) - 1);
			}
			if (ImGui::Button("Duplicate Block") && m_temp_block.size() > 0 && m_listboxCurrentItem < m_temp_block.size() && m_listboxCurrentItem >= 0)
			{
				GenerateCode();
				Model* m = m_temp_block[m_listboxCurrentItem]->getMainModel();
				int did = m_temp_id[m_listboxCurrentItem];
				m_temp_block.push_back(new Block(m_pc, m_prefab[did], m->getPosition(), m->getEulerAngles(), m->getScale()/ m_prefab[did]->meta.mesh_resize_scale, m_temp_block[m_listboxCurrentItem]->getCollision().size() > 0));
				m_temp_id.push_back(did);
				cnames_string.push_back(m_prefab[did]->meta.mesh[0] + " " + std::to_string(cnames_string.size()));
				cnames.push_back(cnames_string[cnames_string.size() - 1].c_str());
				m_listboxCurrentItem = std::max(m_listboxCurrentItem, static_cast<int>(m_temp_block.size()) - 1);
			}
			if (ImGui::Button("Destroy Block") && m_temp_block.size() > 0 && m_listboxCurrentItem < m_temp_block.size() && m_listboxCurrentItem >= 0)
			{
				GenerateCode();
				delete m_temp_block[m_listboxCurrentItem];
				m_temp_block.erase(m_temp_block.begin() + m_listboxCurrentItem);
				m_temp_id.erase(m_temp_id.begin() + m_listboxCurrentItem);
				cnames_string.erase(cnames_string.begin() + m_listboxCurrentItem);
				cnames.erase(cnames.begin() + m_listboxCurrentItem);
				m_listboxCurrentItem = std::min(m_listboxCurrentItem, static_cast<int>(m_temp_block.size()) - 1);
			}
			
			ImGui::PushItemWidth(500);
			ImGui::ListBox("###Block", &m_listboxCurrentItem, cnames.data(), cnames.size(), 8);
			ImGui::PopItemWidth();

			if (m_temp_block.size() > 0)
			{
				ImGui::BeginChild("");
				glm::vec3 pos = m_temp_block[m_listboxCurrentItem]->getMainModel()->getPosition();
				if(ImGui::DragFloat3("Slow Position", (float*)&pos, 0.001f))
				{
					m_temp_block[m_listboxCurrentItem]->getMainModel()->setPosition(pos);
				}
				m_temp_block[m_listboxCurrentItem]->getMainModel()->onGUI();
				m_temp_block[m_listboxCurrentItem]->UpdateTransform();
				m_arrow->setPosition(m_temp_block[m_listboxCurrentItem]->getMainModel()->getPosition());				
				m_arrow->setEulerAngles(m_temp_block[m_listboxCurrentItem]->getMainModel()->getEulerAngles());
				ImGui::EndChild();
			}
		}
	}
	ImGui::End();
}

PrefabLoader::~PrefabLoader()
{
	for (int i = 0; i < m_prefab.size(); i++)
	{
		for (int j = 0; j < m_prefab[i]->sb.size(); j++)
		{
			m_pc.modelManager->destroyBuffer(m_prefab[i]->sb[j]);
		}
		m_prefab[i]->sb.clear();
		for (int j = 0; j < m_prefab[i]->mat.size(); j++)
		{
			m_pc.materialManager->destroyMaterial(m_prefab[i]->mat[j]);
		}
		m_prefab[i]->sb.clear();
		delete m_prefab[i];
	}
	m_prefab.clear();
	for (auto it = m_mapTexture.begin(); it != m_mapTexture.end(); ++it) 
	{
		m_pc.textureManager->destroyTexture(it->second);
	}
	m_mapTexture.clear();
}