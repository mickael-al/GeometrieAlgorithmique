#ifndef __BONE_MANAGER__
#define __BONE_MANAGER__

#include "ImguiBlock.hpp"
#include "Behaviour.hpp"
#include "GameEngine.hpp"
#include <string>
#include <vector>

struct BoundingBox
{
	Model* box;
	Materials* material;
	std::string name;
};

class BoneManager : public ImguiBlock, public Behaviour
{
public:
	BoneManager(Model* plane);
	void start();
	void fixedUpdate();
	void update();
	void stop();
	void onGUI();
	void preRender(VulkanMisc* vM);
	void render(VulkanMisc* vM);
	glm::vec3 directionToRotation(glm::vec3 direction);
	Model* createSegment3D(glm::vec3 p1, glm::vec3 p2, Materials* mat = nullptr);
private:
	ptrClass m_pc;
	Camera* m_cam2D;
	Camera* m_cam3D;
	
	bool m_open = false;
	bool m_cameraChange = false;
	bool m_planeActive = true;
	bool m_priority = false;
	bool m_planeShowHide = false;
	bool m_isMouseOverUI = false;
	bool m_create = false;
	bool m_create_box = false;
	bool m_delete_box = false;
	bool m_delete = false;
	float m_size = 0.05f;
	int selectedItem = 0;

	std::vector<BoundingBox*> m_bb;
	std::vector<const char*> m_pointListbb;
	ShapeBuffer* m_shapeBuffer_box = nullptr;
	int m_selectedItemBB = 0;
	float m_sliderScale = 1.0f;
	Model* m_modelCurrent = nullptr;
	Materials * m_matCurrent = nullptr;
	ShapeBuffer* m_shapeCurrent = nullptr;
	ShapeBuffer* m_sb = nullptr;
	Materials* m_pointMat = nullptr;
	Materials* m_segmentMat = nullptr;
	Model* m_plane;
};

#endif//!__BONE_MANAGER__