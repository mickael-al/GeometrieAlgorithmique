#include "GraphiquePipelineManager.hpp"

namespace Ge
{
	std::vector<GraphiquePipeline *> GraphiquePipelineManager::m_graphiquePipeline;
	bool GraphiquePipelineManager::initialize(VulkanMisc *vM)
	{
		vulkanM = vM;
		for (int i = 0; i < m_fileNameShaders.size(); i++)
		{
			createPipeline(m_fileNameShaders[i]->Frag, m_fileNameShaders[i]->Vert, m_fileNameShaders[i]->back, m_fileNameShaders[i]->multiSampling, m_fileNameShaders[i]->transparency, m_fileNameShaders[i]->cullMode);
		}
		if (m_fileNameShaders.size() == 0)
		{
			GPData* gpd = new GPData();
			createPipeline("../Shader/frag.spv", "../Shader/vert.spv");
		}
		return true;
	}

	void GraphiquePipelineManager::release()
	{
		m_fileNameShaders.clear();
		for (int i = 0; i < m_graphiquePipeline.size();i++)
		{			
			m_fileNameShaders.push_back(m_graphiquePipeline[i]->getShaderPair());
			delete (m_graphiquePipeline[i]);
		}
		m_graphiquePipeline.clear();
	}

	GraphiquePipeline * GraphiquePipelineManager::createPipeline(const std::string &frag, const std::string &vert,bool back,bool multiS,bool transparency,int cullmode)
	{
		ShaderPair * sp = new ShaderPair(frag, vert, back,multiS, transparency, cullmode);
		GraphiquePipeline * gp = new GraphiquePipeline(vulkanM, sp);
		m_graphiquePipeline.push_back(gp);				
		return gp;
	}

	/*GPData* GraphiquePipelineManager::createPipelineD(const std::string& frag, const std::string& vert, bool back = false, bool multiS = true, bool transparency = false, int cullmode = 1)
	{
		GPData * gpd = new GPData();
		gpd->gp = createPipeline(frag, vert, back,multiS, transparency,cullmode);
		m_gpdata.push_back(gpd);
		return gpd;
	}*/


	std::vector<GraphiquePipeline *> & GraphiquePipelineManager::GetPipelines()
	{
		return m_graphiquePipeline;
	}

	void GraphiquePipelineManager::destroyPipeline(GraphiquePipeline * pipeline)
	{
		m_graphiquePipeline.erase(std::remove(m_graphiquePipeline.begin(), m_graphiquePipeline.end(), pipeline), m_graphiquePipeline.end());
		delete (pipeline);
	}
}