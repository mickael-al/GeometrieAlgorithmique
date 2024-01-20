#include "BoneManager.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/intersect.hpp>
#include <iostream>
#include <unordered_set>
#include <functional>

BoneManager::BoneManager(Model* plane)
{
	m_plane = plane;
	m_planeShowHide = true;
}

glm::vec3 BoneManager::directionToRotation(glm::vec3 direction)
{
	direction = glm::normalize(direction);
	glm::quat quaternion = glm::rotation(glm::vec3(1.0f, 0.0f, 0.0f), direction);
	return glm::degrees(glm::eulerAngles(quaternion));
}

Model* BoneManager::createSegment3D(glm::vec3 p1, glm::vec3 p2, Materials* mat)
{
	glm::vec3 pos = (p1 + p2) / 2.0f;
	float scale = glm::distance(p1, p2);
	Model* m = m_pc.modelManager->createModel(m_sb);
	m->setMaterial(mat == nullptr ? m_segmentMat : mat);
	m->setPosition(pos);
	m->setScale(glm::vec3(scale, m_size, m_size));
	glm::vec3 direction = p2 - p1;
	m->setEulerAngles(directionToRotation(direction));
	return m;
}


void BoneManager::start()
{
	m_pc = GameEngine::getPtrClass();
	m_pc.lightManager->createPointLight(glm::vec3(10.0f, 5.0f, -10.0f), glm::vec3(1.0f));
	m_pc.lightManager->createPointLight(glm::vec3(-10.0f, 5.0f, 10.0f), glm::vec3(1.0f));
	m_pc.lightManager->createPointLight(glm::vec3(-10.0f, 5.0f, -10.0f), glm::vec3(1.0f));
	m_pc.lightManager->createPointLight(glm::vec3(10.0f, 5.0f, 10.0f), glm::vec3(1.0f));

	m_cam3D = m_pc.cameraManager->getCurrentCamera();
	m_cam2D = m_pc.cameraManager->createCamera();
	m_cam2D->setOrtho(true);
	m_cam2D->setPosition(glm::vec3(0, 10, 0));
	m_cam2D->setEulerAngles(glm::vec3(-90.0f, 0.0f, 0.0f));
	m_cam2D->setOrthoSize(15.0f);

	m_sb = m_pc.modelManager->allocateBuffer("../Model/cube.obj");
	GraphiquePipeline * gp_unlit = m_pc.graphiquePipelineManager->createPipeline("../Shader/frag_unlit.spv", "../Shader/vert_unlit.spv");
	GraphiquePipeline * gp_wire = m_pc.graphiquePipelineManager->createPipeline("../Shader/shader_wireframe_fs.spv", "../Shader/vert.spv");

	Model * cube = m_pc.modelManager->createModel(m_sb);
	cube->setScale(glm::vec3(0));
	m_pointMat = m_pc.materialManager->createMaterial();
	m_pointMat->setColor(glm::vec3(0.0f, 0.0f, 1.0f));
	m_pointMat->setMetallic(0.7f);
	m_pointMat->setRoughness(0.15f);
	m_pointMat->setPipeline(gp_unlit);
}

void BoneManager::fixedUpdate()
{

}

void BoneManager::update()
{
	if (m_cameraChange)
	{
		m_cameraChange = false;
		m_priority = !m_priority;
		m_cam3D->setPriority(m_priority ? 0 : 1);
		m_cam2D->setPriority(m_priority ? 1 : 0);
		m_pc.cameraManager->updatePriorityCamera();
	}
	if (m_planeActive)
	{
		m_planeActive = false;
		m_plane->setScale(m_planeShowHide ? glm::vec3(10.0f, 1.0f, 10.0f) : glm::vec3(0.0f));
	}
	if (m_modelCurrent)
	{
		m_modelCurrent->setScale(glm::vec3(m_sliderScale));
		if (m_delete)
		{
			m_pc.modelManager->destroyModel(m_modelCurrent);
			m_pc.modelManager->destroyBuffer(m_shapeCurrent);
			m_pc.materialManager->destroyMaterial(m_matCurrent);
			m_delete = false;
			m_modelCurrent = nullptr;
		}
	}
	if (m_create)
	{
		m_shapeCurrent = m_pc.modelManager->allocateBufferWire(selectedItem == 0 ? "../Assets/spaghet.obj" : selectedItem == 1 ? "../Assets/MIddlePoly.obj" : "../Assets/UltraLowPoly.obj");
		m_modelCurrent = m_pc.modelManager->createModel(m_shapeCurrent);		
		m_matCurrent = m_pc.materialManager->createMaterial();
		m_matCurrent->setPipelineIndex(3);
		m_modelCurrent->setMaterial(m_matCurrent);
		m_create = false;
		m_delete = false;
	}	
}
	}

	if (m_create_box)
	{
		BoundingBox* new_box = new BoundingBox();
		new_box->box = m_pc.modelManager->createModel(m_sb);
		m_bb.push_back(new_box);
		new_box->name = "Box" + std::to_string(m_bb.size());
		m_pointListbb.push_back(new_box->name.c_str());
		m_create_box = false;
	}

	if (m_delete_box)
	{
		BoundingBox* current_box = m_bb[m_selectedItemBB];
		Debug::Log("%d", m_selectedItemBB);
		m_pc.modelManager->destroyModel(m_bb[m_selectedItemBB]->box);
		m_bb.erase(std::remove(m_bb.begin(), m_bb.end(), m_bb[m_selectedItemBB]), m_bb.end());
		m_pointListbb.erase(m_pointListbb.begin() + m_selectedItemBB);
		m_selectedItemBB = 0;
		m_delete_box = false;
	}

}

void BoneManager::stop()
{

}

void BoneManager::preRender(VulkanMisc* vM)
{

}

void BoneManager::render(VulkanMisc* vM)
{
	ImGuiIO& io = ImGui::GetIO();
	ImGuiWindowFlags window_flags = ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoFocusOnAppearing;

	if (ImGui::Begin("Extrusion", &m_open, window_flags))
	{
		ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), "Extrusion");
		ImVec2 windowPos = ImGui::GetWindowPos();
		ImVec2 windowSize = ImGui::GetWindowSize();
		ImVec2 rectMin = windowPos;
		ImVec2 rectMax = ImVec2(windowPos.x + windowSize.x, windowPos.y + windowSize.y);
		m_isMouseOverUI = ImGui::IsMouseHoveringRect(rectMin, rectMax);
		if (ImGui::Checkbox("Change Camera", &m_cameraChange)) {}
		if (ImGui::Checkbox("Show/Hide Plane", &m_planeShowHide))
		{
			m_planeActive = true;
		}
		ImGui::Combo("Model", &selectedItem, "Spaghet\0MIddlePoly\0UltraLowPoly\0");
		ImGui::DragFloat("Scale", &m_sliderScale,0.05f);
		if (ImGui::Button("Create Mesh"))
		{
			m_create = true;		
			m_delete = true;
		}		
		if (ImGui::Button("Delete Mesh"))
		{
			m_delete = true;
		}
		if (m_bb.size() > 0)
		{
			ImGui::ListBox("Liste de BoundingBox", &m_selectedItemBB, m_pointListbb.data(), m_pointListbb.size());
			m_bb[m_selectedItemBB]->box->onGUI();
		}
		if (ImGui::Button("Create Box"))
		{
			m_create_box = true;
		}

		if (ImGui::Button("Delete Box"))
		{
			m_delete_box = true;
		}
		ImGui::Spacing();
	}
	ImGui::End();
}

void BoneManager::onGUI()
{

}