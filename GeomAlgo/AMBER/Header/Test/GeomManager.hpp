#ifndef __GEOMETRIE_MANAGER__
#define __GEOMETRIE_MANAGER__

#include "ImguiBlock.hpp"
#include "Behaviour.hpp"
#include "GameEngine.hpp"
#include <string>
#include <vector>

class GeomManager : public ImguiBlock, public Behaviour
{
public:
	GeomManager(Model* plane);
	void start();
	void fixedUpdate();
	void update();
	void stop();
	void onGUI();
	void preRender(VulkanMisc* vM);
	void render(VulkanMisc* vM);
	Model* createSegment(glm::vec2 p1, glm::vec2 p2);
	glm::vec3 directionToRotation(glm::vec3 direction);
	float getAngle(glm::vec2 a, glm::vec2 b);
	void marcheJarvis();
	void grahamScan();
	void Triangulation2D();
private:
	ptrClass m_pc;
	Camera* m_cam2D;
	Camera* m_cam3D;

	bool m_open = false;
	bool m_cameraChange = false;
	bool m_planeShowHide = false;
	bool m_convexHullChange = false;
	bool m_convexHullCheck = true;
	bool m_planeActive = true;
	bool m_priority = false;
	bool m_cloudPoint = false;
	bool m_clearCloudPoint = false;
	bool m_rotate = false;
	bool m_convexHull = false;
	bool m_marcheJarvis = false;
	bool m_cloudButton = false;
	bool m_isMouseOverUI = false;


	float m_size = 0.05f;
	int m_control_point = -1;
	int m_courbe_point = 0;
	int m_listboxCurrentItem = 0;
	int m_listboxCurrentItemExtrusion = 0;
	std::vector<const char*> cnames;
	std::vector<const char*> cExtrusionNames;
	std::vector<std::string> valueS;
	std::vector<std::string> valueSExtrusion;

	std::vector<Model*> m_points_clouds;
	std::vector<PointLight*> m_points_light_clouds;

	std::vector<Model*> m_segments;

	ShapeBuffer* m_sb = nullptr;
	Materials* m_pointMat = nullptr;
	Materials* m_segmentMat = nullptr;	
	std::vector<Model*> m_extrusionList;
	Model* m_plane;
};

#endif//!__GEOMETRIE_MANAGER__