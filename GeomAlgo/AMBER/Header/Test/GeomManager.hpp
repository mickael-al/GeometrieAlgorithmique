#ifndef __GEOMETRIE_MANAGER__
#define __GEOMETRIE_MANAGER__

#include "ImguiBlock.hpp"
#include "Behaviour.hpp"
#include "GameEngine.hpp"
#include <string>
#include <vector>

struct Sommet3D;
struct Triangle;
struct Edge3D;
struct Face;
struct Tetra;
struct VoronoiEdge;

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
	Model* createSegment(glm::vec2 p1, glm::vec2 p2, Materials* mat = nullptr);
	Model* createSegment3D(glm::vec3 p1, glm::vec3 p2, Materials* mat = nullptr);
	void updateSegment(glm::vec2 p1, glm::vec2 p2, Model* m);
	std::vector<glm::uvec2> grahamScan(std::vector<glm::vec2> point2D);
	std::vector<VoronoiEdge> calculateVoronoiDiagram(const std::vector<Triangle>& delaunayTriangles, std::vector<glm::vec2> points2d);
	glm::vec3 directionToRotation(glm::vec3 direction);
	float getAngle(glm::vec2 a, glm::vec2 b);
	void marcheJarvis();
	void grahamScan();
	void Triangulation2D();
	void DelaunayTriangulation2D();
	void DelaunayTriangulation2DTest();
	void DelaunayNoyauxTriangulation();
	void Voronoi();
	void ConvexHull3D();
private:
	ptrClass m_pc;
	Camera* m_cam2D;
	Camera* m_cam3D;

	bool m_open = false;
	bool m_cameraChange = false;
	bool m_planeShowHide = false;
	bool m_convexHullChange = false;
	bool m_convexHull3DChange = false;
	bool m_moveVoronoi = false;
	bool m_delaunayNoyaux = false;
	bool m_convexHullCheck = true;
	bool m_convexHull3DCheck = false;
	bool m_delaunayNoyauxCheck = false;
	bool m_planeActive = true;
	bool m_priority = false;
	bool m_cloudPoint = false;
	bool m_clearCloudPoint = false;
	bool m_rotate = false;
	bool m_convexHull = false;
	bool m_convexHull3D = false;
	bool m_marcheJarvis = false;
	bool m_cloudButton = false;
	bool m_isMouseOverUI = false;


	float m_size = 0.05f;
	float m_height = 0.0f;
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
	std::vector<Model*> m_segmentsVoronoi;
	std::vector<glm::uvec3> currentTriangle;
	std::vector<Sommet3D*> sommets;
	std::vector<Edge3D*> edges;
	std::vector<Face*> faces;
	std::vector<Tetra*> tetra;

	ShapeBuffer* m_sb = nullptr;
	Materials* m_pointMat = nullptr;
	Materials* m_segmentMat = nullptr;	
	Materials* m_segmentMat2 = nullptr;
	Materials* m_segmentMat3 = nullptr;
	std::vector<Model*> m_extrusionList;
	Model* m_plane;
};

#endif//!__GEOMETRIE_MANAGER__