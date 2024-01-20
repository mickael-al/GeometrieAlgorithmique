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
	GraphiquePipeline* gp_unlit = m_pc.graphiquePipelineManager->createPipeline("../Shader/frag_unlit.spv", "../Shader/vert_unlit.spv");

	Model* cube = m_pc.modelManager->createModel(m_sb);
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
	if (m_create)
	{
		m_shapeCurrent = m_pc.modelManager->allocateBuffer(selectedItem == 0 ? "../Assets/spaghet.obj" : selectedItem == 1 ? "../Assets/MIddlePoly.obj" : "../Assets/UltraLowPoly.obj");
		m_modelCurrent = m_pc.modelManager->createModel(m_shapeCurrent);
		m_create = false;
	}

}

void BoneManager::createBones(std::vector<glm::vec3> cloud_point)
{
	glm::vec3 bary;

	for (int i = 0; i < cloud_point.size(); i++)
	{
		bary += cloud_point[i];
	}
	bary = bary / float(cloud_point.size());
	for (int i = 0; i < cloud_point.size(); i++)
	{
		cloud_point[i] -= bary;
	}

	glm::mat3 covarMat = covarianceMatrix(cloud_point);

	glm::vec3 properVec = cloud_point[0];
	float maxDelta = glm::vec3(covarMat * cloud_point[0]).x/ cloud_point[0].x;

	for (int i = 1; i < cloud_point.size(); i++)
	{
		float delta = glm::vec3(covarMat * cloud_point[i]).x / cloud_point[i].x;
		if (delta > maxDelta)
		{
			maxDelta = delta;
			properVec = cloud_point[i];
		}
	}

	std::vector<glm::vec3> proj_point;
	for (int i = 0; i < cloud_point.size(); i++)
	{
		proj_point.push_back(((cloud_point[i] * properVec) / (glm::abs(properVec) * glm::abs(properVec))) * properVec);
	}

	glm::vec3 BMin = glm::vec3(0);
	glm::vec3 CMax = glm::vec3(0);
	for (int i = 0; i < proj_point.size(); i++)
	{
		if (glm::dot(properVec, proj_point[i]) < 0)
		{
			if (glm::distance(proj_point[i], glm::vec3(0)) > glm::distance(BMin, glm::vec3(0)))
			{
				BMin = proj_point[i];
			}
		}
		else
		{
			if (glm::distance(proj_point[i], glm::vec3(0)) > glm::distance(CMax, glm::vec3(0)))
			{
				CMax = proj_point[i];
			}
		}
	}
}

glm::mat3 BoneManager::covarianceMatrix(std::vector<glm::vec3> cloud_point)
{
	float moyX = 0.0f, moyY = 0.0f, moyZ = 0.0f, varX = 0.0f, varY = 0.0f, varZ = 0.0f, coXY = 0.0f, coXZ = 0.0f, coYZ = 0.0f;
	for (int i = 0; i < cloud_point.size(); i++)
	{
		moyX += cloud_point[i].x;
		moyY += cloud_point[i].y;
		moyZ += cloud_point[i].z;
	}
	moyX /= cloud_point.size();
	moyY /= cloud_point.size();
	moyZ /= cloud_point.size();

	for (int i = 0; i < cloud_point.size(); i++)
	{
		varX += std::pow(cloud_point[i].x - moyX, 2);
		varY += std::pow(cloud_point[i].y - moyY, 2);
		varZ += std::pow(cloud_point[i].z - moyZ, 2);
		coXY += (cloud_point[i].x - moyX) * (cloud_point[i].y - moyY);
		coXZ += (cloud_point[i].y - moyY) * (cloud_point[i].z - moyZ);
		coYZ += (cloud_point[i].x - moyX) * (cloud_point[i].z - moyZ);
	}
	varX /= cloud_point.size();
	varY /= cloud_point.size();
	varZ /= cloud_point.size();

	coXY /= cloud_point.size();
	coXZ /= cloud_point.size();
	coYZ /= cloud_point.size();

	return glm::mat3(varX, coXY, coXZ,
					 coXY, varY, coYZ,
					 coXZ, coYZ, varZ);
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
		if (ImGui::Button("Create"))
		{
			m_create = true;			
		}
		if (m_bb.size() > 0)
		{
			ImGui::ListBox("Liste de BoundingBox", &m_selectedItemBB, m_pointListbb.data(), m_pointListbb.size());
		}
		if (ImGui::Button("Create Box"))
		{
			m_pointListbb.push_back(std::to_string(m_bb.size()).c_str());
		}

		ImGui::Spacing();
	}
	ImGui::End();
}

void BoneManager::onGUI()
{

}