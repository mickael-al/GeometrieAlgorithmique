#include "GeomManager.hpp"
#include "CircularList.hpp"
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/epsilon.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/intersect.hpp>
#include <iostream>

GeomManager::GeomManager(Model* plane)
{
	m_plane = plane;
	m_planeShowHide = true;
}

glm::vec3 GeomManager::directionToRotation(glm::vec3 direction)
{
	direction = glm::normalize(direction);
	glm::quat quaternion = glm::rotation(glm::vec3(1.0f, 0.0f, 0.0f), direction);
	return glm::degrees(glm::eulerAngles(quaternion));
}

Model* GeomManager::createSegment(glm::vec2 p1, glm::vec2 p2)
{
	glm::vec2 pos = (p1 + p2) / 2.0f;
	float scale = glm::distance(p1, p2);
	Model* m = m_pc.modelManager->createModel(m_sb);
	m->setMaterial(m_segmentMat);
	m->setPosition(glm::vec3(pos.x, 0.0f, pos.y));
	m->setScale(glm::vec3(scale, m_size, m_size));
	glm::vec2 direction = p2 - p1;
	m->setEulerAngles(directionToRotation(glm::vec3(direction.x, 0.0f, direction.y)));
	return m;
}


float GeomManager::getAngle(glm::vec2 a, glm::vec2 b)
{
	float dot = glm::dot(a, b);
	float magA = glm::length(a);
	float magB = glm::length(b);
	float angle = std::acos(dot / (magA * magB));
	glm::vec3 crossProduct = glm::cross(glm::vec3(a, 0.0f), glm::vec3(b, 0.0f));
	if (crossProduct.z < 0.0f) {
		angle = 2.0f * glm::pi<float>() - angle;
	}

	return angle;
}

bool onSegment(glm::vec2 p, glm::vec2 q, glm::vec2 r)
{
	if (q.x <= max(p.x, r.x) && q.x >= min(p.x, r.x) && q.y <= max(p.y, r.y) && q.y >= min(p.y, r.y))	
	{
		return true;
	}

	return false;
}


int orientation(glm::vec2 p, glm::vec2 q, glm::vec2 r)
{
	int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0;

	return (val > 0) ? 1 : 2;
}


bool SegmentsIntersect(glm::vec2 p1, glm::vec2 q1, glm::vec2 p2, glm::vec2 q2)
{
	int o1 = orientation(p1, q1, p2);
	int o2 = orientation(p1, q1, q2);
	int o3 = orientation(p2, q2, p1);
	int o4 = orientation(p2, q2, q1);

	if (o1 != o2 && o3 != o4)
		return true;

	if (o1 == 0 && onSegment(p1, p2, q1)) return true;
	if (o2 == 0 && onSegment(p1, q2, q1)) return true;
	if (o3 == 0 && onSegment(p2, p1, q2)) return true;
	if (o4 == 0 && onSegment(p2, q1, q2)) return true;

	return false;
}

const double EPSILON = 1e-12;

struct Edge {
	int a, b;

	Edge(int a, int b) : a(a), b(b) {}
};

struct Point {
	double x, y;

	Point(double x, double y) : x(x), y(y) {}
};

double crossProduct(const Point& A, const Point& B) {
	return A.x * B.y - A.y * B.x;
}

bool insideCircle(const Point& A, const Point& B, const Point& C, const Point& P) {
	double ax = A.x - P.x;
	double ay = A.y - P.y;
	double bx = B.x - P.x;
	double by = B.y - P.y;
	double cx = C.x - P.x;
	double cy = C.y - P.y;

	double a2 = ax * ax + ay * ay;
	double b2 = bx * bx + by * by;
	double c2 = cx * cx + cy * cy;

	return (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) >= EPSILON;
}

std::vector<Edge> triangulate(std::vector<Point>& points) {
	int n = points.size();
	std::vector<Edge> edges;

	// Sorting the points by x-coordinate
	std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
		return a.x <= b.x;
		});

	// Creating the lower hull
	for (int i = 0; i < n; i++) {
		while (edges.size() >= 2) {
			int j = edges.size() - 2;
			int k = edges.size() - 1;
			Point A = points[edges[j].a];
			Point B = points[edges[j].b];
			Point C = points[edges[k].b];

			if (crossProduct(Point(B.x - A.x, B.y - A.y), Point(C.x - B.x, C.y - B.y)) > 0) {
				break;
			}

			edges.pop_back();
		}
		edges.push_back(Edge(edges.size(), i));
	}
	int lower = edges.size();

	// Creating the upper hull
	for (int i = n - 2, t = lower + 1; i >= 0; i--) {
		while (edges.size() >= t) {
			int j = edges.size() - 2;
			int k = edges.size() - 1;
			Point A = points[edges[j].a];
			Point B = points[edges[j].b];
			Point C = points[edges[k].b];

			if (crossProduct(Point(B.x - A.x, B.y - A.y), Point(C.x - B.x, C.y - B.y)) > 0) {
				break;
			}

			edges.pop_back();
		}
		edges.push_back(Edge(i, edges.size()));
	}

	// Removing the duplicate edges from the hull
	edges.pop_back();

	// Creating the triangulation
	std::vector<Edge> result;
	for (int i = 0; i < edges.size(); i++) {
		int a = edges[i].a;
		int b = edges[i].b;
		Point A = points[a];
		Point B = points[b];
		bool flag = true;

		for (int j = 0; j < n; j++) {
			if (j == a || j == b) {
				continue;
			}
			Point P = points[j];
			if (insideCircle(A, B, P, points[(a + b) >> 1])) {
				flag = false;
				break;
			}
		}
		if (flag) {
			result.push_back(Edge(a, b));
		}
	}

	return result;
}

void GeomManager::Triangulation2D()
{
	for (int i = 0; i < m_segments.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_segments[i]);
	}
	m_segments.clear();
	std::vector<glm::vec2> point2D;
	for (int i = 0; i < m_points_clouds.size(); i++)
	{
		point2D.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
	}

	std::sort(point2D.begin(), point2D.end(), [](const glm::vec2& a, const glm::vec2& b)
		{
		return a.x < b.x || (a.x == b.x && a.y < b.y);
		});

	std::vector<glm::uvec3> triangulation;
	std::vector<glm::uvec2> edges;	
	edges.push_back(glm::uvec2(0, 1));
	edges.push_back(glm::uvec2(1, 2));
	edges.push_back(glm::uvec2(2, 0));
	triangulation.push_back(glm::uvec3(0,1,2));

	for (int i = 3; i < point2D.size(); i++)
	{
		std::vector<glm::uvec2> new_edges;
		for (int j = 0; j < i; j++)
		{
			bool isIntersect = false;
			for (int k = 0; k < edges.size(); k++)
			{
				if (j != edges[k].x && j != edges[k].y)
				{
					if (SegmentsIntersect(point2D[i]*100.0f, point2D[j]*100.0f, point2D[edges[k].x]*100.0f, point2D[edges[k].y]*100.0f))
					{
						isIntersect = true;
						break;
					}
				}
			}
			if (!isIntersect)
			{
				new_edges.push_back(glm::uvec2(i,j));
			}
		}		
		for (int j = 0; j < new_edges.size(); j++)
		{
			edges.push_back(new_edges[j]);
		}			
	}

	for (int i = 0; i < edges.size(); i++)
	{
		m_segments.push_back(createSegment(point2D[edges[i].x], point2D[edges[i].y]));
	}
}

void GeomManager::grahamScan()
{
	for (int i = 0; i < m_segments.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_segments[i]);
	}
	m_segments.clear();

	std::vector<glm::vec2> point2D;
	std::vector<int> buffer2D;
	std::vector<int> order;
	for (int i = 0; i < m_points_clouds.size(); i++)
	{
		point2D.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
		buffer2D.push_back(i);
	}
	
	glm::vec2 barycentre = glm::vec2(0);
	for (int i = 0; i < point2D.size(); i++)
	{
		barycentre += point2D[i];
	}
	barycentre /= point2D.size();
	glm::vec2 v = barycentre + glm::vec2(1, 0);
	float angleMin = getAngle(v, point2D[buffer2D[0]] - barycentre);
	float lenghtMin = glm::length(barycentre - point2D[buffer2D[0]]);
	int iNew = 0;
	for (int i = 1; i < buffer2D.size(); i++)
	{
		float angle = getAngle(v, point2D[buffer2D[i]] - barycentre);
		if (angleMin > angle || (angleMin == angle && lenghtMin > glm::length(barycentre - point2D[buffer2D[i]])))
		{
			angleMin = angle;
			lenghtMin = glm::length(barycentre - point2D[buffer2D[i]]);
			iNew = i;
		}
	}
	order.push_back(buffer2D[iNew]);
	buffer2D.erase(buffer2D.begin() + iNew);
	while (buffer2D.size() > 0)
	{
		angleMin = getAngle(v, point2D[buffer2D[0]] - barycentre);
		lenghtMin = glm::length(barycentre - point2D[buffer2D[0]]);
		iNew = 0;
		for (int i = 1; i < buffer2D.size(); i++)
		{

			float angle = getAngle(v, point2D[buffer2D[i]] - barycentre);
			if (angleMin > angle || (angleMin == angle && lenghtMin > glm::length(barycentre - point2D[buffer2D[i]])))
			{
				angleMin = angle;
				lenghtMin = glm::length(barycentre - point2D[buffer2D[i]]);
				iNew = i;
				
			}
		}
		order.push_back(buffer2D[iNew]);
		buffer2D.erase(buffer2D.begin() + iNew);
	}
	circular_doubly_linked_list cdll1;
	for (int i = 0; i < order.size(); i++)
	{
		cdll1.append_node(order[i]);
	}
	auto s = cdll1.head;
	auto pivot = s;
	bool next = false;
	do
	{
		float crossProduct =	(point2D[pivot->data].x - point2D[pivot->prev->data].x) * (point2D[pivot->next->data].y - point2D[pivot->data].y) - 
								(point2D[pivot->data].y - point2D[pivot->prev->data].y) * (point2D[pivot->next->data].x - point2D[pivot->data].x);
		if (crossProduct >= 0.0f)
		{
			pivot = pivot->next;
			next = true;
		}
		else
		{
			s = pivot->prev;
			cdll1.delete_node(pivot);
			pivot = s;
			next = false;
		}
		Debug::Log("%d", pivot->data);
	} while (pivot != s || next == false);
	Node* tmp = cdll1.head;
	while (tmp->next != cdll1.head)
	{
		m_segments.push_back(createSegment(point2D[tmp->data], point2D[tmp->next->data]));
		tmp = tmp->next;
		Debug::Log("%d", tmp->data);

	}
	m_segments.push_back(createSegment(point2D[tmp->data], point2D[tmp->next->data]));
	//delete(&cdll1);
}

void GeomManager::marcheJarvis()
{
	for (int i = 0; i < m_segments.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_segments[i]);
	}
	m_segments.clear();
	int min_id = 0;

	std::vector<glm::vec2> point2D;
	for (int i = 0; i < m_points_clouds.size(); i++)
	{
		point2D.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
	}
	glm::vec2 pos_id = point2D[min_id];
	for (int i = 1; i < point2D.size(); i++)
	{
		glm::vec2 pos = point2D[i];
		if (pos.x < pos_id.x || (pos.x == pos_id.x && pos.y < pos_id.y))
		{
			min_id = i;
			pos_id = point2D[min_id];
		}
	}
	glm::vec2 v = glm::vec2(0, -1);
	std::vector<int> order;
	int i = min_id;
	int j = 0;
	do
	{
		order.push_back(i);
		if (i == 0)
		{
			j = 1;
		}
		else
		{
			j = 0;
		}
		float angleMin = getAngle(v, point2D[j] - point2D[i]);
		float lenghtMax = glm::length(point2D[i] - point2D[j]);
		int iNew = j;

		for (j = iNew + 1; j < point2D.size(); j++)
		{
			if (j != i) {
				float angle = getAngle(v, point2D[j] - point2D[i]);
				if (angleMin > angle || (angleMin == angle && lenghtMax < glm::length(point2D[i] - point2D[j])))
				{
					angleMin = angle;
					lenghtMax = glm::length(point2D[i] - point2D[j]);
					iNew = j;
				}
			}
		}

		v = point2D[iNew] - point2D[i];
		i = iNew;
	} while (i != min_id);

	for (int i = 0; i < order.size() - 1; i++)
	{
		m_segments.push_back(createSegment(point2D[order[i]], point2D[order[i + 1]]));
	}
	m_segments.push_back(createSegment(point2D[order[order.size() - 1]], point2D[order[0]]));
}

void GeomManager::start()
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

	m_segmentMat = m_pc.materialManager->createMaterial();
	m_segmentMat->setColor(glm::vec3(1.0f, 0.0f, 0.0f));
	m_segmentMat->setMetallic(0.7f);
	m_segmentMat->setRoughness(0.15f);
	m_segmentMat->setPipeline(gp_unlit);
}

void GeomManager::fixedUpdate()
{

}

void GeomManager::update()
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
	if (m_convexHullChange)
	{
		m_convexHullChange = false;
		m_convexHull = !m_convexHull;
	}
	if (m_clearCloudPoint)
	{
		m_clearCloudPoint = false;
		for (int i = 0; i < m_points_clouds.size(); i++)
		{
			m_pc.modelManager->destroyModel(m_points_clouds[i]);
			//m_pc.lightManager->destroyLight(m_points_light_clouds[i]);
		}
		for (int i = 0; i < m_segments.size(); i++)
		{
			m_pc.modelManager->destroyModel(m_segments[i]);
		}
		m_segments.clear();
		m_points_clouds.clear();
		m_points_light_clouds.clear();
	}
	if (m_cloudPoint && !m_priority && !m_isMouseOverUI)
	{
		if (m_pc.inputManager->getMouse(0))
		{
			if (m_cloudButton)
			{
				Model* m = m_pc.modelManager->createModel(m_sb);
				m->setMaterial(m_pointMat);
				m->setScale(glm::vec3(m_size));
				glm::vec2 mp = glm::vec2(m_pc.inputManager->getMousePosX(), m_pc.inputManager->getMousePosY());
				glm::vec2 ss = glm::vec2(m_pc.settingManager->getWindowWidth(), m_pc.settingManager->getWindowHeight());
				glm::vec2 pos = m_cam2D->ScreenToSpace(mp, ss);
				m->setPosition(glm::vec3(pos.x, 0.0f, pos.y));
				m_points_clouds.push_back(m);
				m_cloudButton = false;
				//PointLight* pl = m_pc.lightManager->createPointLight(glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
				//pl->setRange(1.0f);
				//pl->setPosition(m->getPosition());
				//m_points_light_clouds.push_back(pl);
				if (m_points_clouds.size() > 2)
				{
					if (!m_convexHull) 
					{
						if (m_marcheJarvis)
						{
							marcheJarvis();
						}
						else
						{
							grahamScan();
						}
					}
					else
					{
						Triangulation2D();
					}
				}
			}
		}
		else
		{
			m_cloudButton = true;
		}
	}


}

void GeomManager::stop()
{

}

void GeomManager::preRender(VulkanMisc* vM)
{

}

void GeomManager::render(VulkanMisc* vM)
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
		if (m_cloudPoint)
		{
			if (ImGui::Button("Close"))
			{
				m_cloudPoint = false;
			}
		}
		else
		{
			if (ImGui::Button("Create Cloud Point"))
			{
				m_clearCloudPoint = m_points_clouds.size() > 0;
				m_cloudPoint = true;
				m_cloudButton = false;
			}
			if (ImGui::Checkbox("Create Convex Hull", &m_convexHullCheck))
			{
				m_convexHullChange = true;
			}

			if (ImGui::Button("Marche de Jarvis"))
			{
				m_marcheJarvis = true;
			}

			if (ImGui::Button("Graham-Scan"))
			{
				m_marcheJarvis = false;
			}
			if (ImGui::Button("Triangulation 2D"))
			{
				if (m_points_clouds.size() > 2)
				{
					Triangulation2D();
				}
			}
		}
		cnames.clear();
		cExtrusionNames.clear();
		valueS.clear();
		valueSExtrusion.clear();

	
		

		ImGui::Spacing();
	}
	ImGui::End();
}

void GeomManager::onGUI()
{

}