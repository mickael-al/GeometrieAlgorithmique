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
#include <unordered_set>
#include <functional>

#define ROUGE 0
#define VIOLET 1
#define BLEU 2

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

void GeomManager::updateSegment(glm::vec2 p1, glm::vec2 p2, Model* m)
{
	glm::vec2 pos = (p1 + p2) / 2.0f;
	float scale = glm::distance(p1, p2);
	m->setPosition(glm::vec3(pos.x, 0.0f, pos.y));
	m->setScale(glm::vec3(scale, m_size, m_size));
	glm::vec2 direction = p2 - p1;
	m->setEulerAngles(directionToRotation(glm::vec3(direction.x, 0.0f, direction.y)));	
}

Model* GeomManager::createSegment(glm::vec2 p1, glm::vec2 p2,Materials * mat)
{
	glm::vec2 pos = (p1 + p2) / 2.0f;
	float scale = glm::distance(p1, p2);
	Model* m = m_pc.modelManager->createModel(m_sb);
	m->setMaterial(mat == nullptr ? m_segmentMat : mat);
	m->setPosition(glm::vec3(pos.x, 0.0f, pos.y));
	m->setScale(glm::vec3(scale, m_size, m_size));
	glm::vec2 direction = p2 - p1;
	m->setEulerAngles(directionToRotation(glm::vec3(direction.x, 0.0f, direction.y)));
	return m;
}

Model* GeomManager::createSegment3D(glm::vec3 p1, glm::vec3 p2, Materials* mat)
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

double crossProduct(const glm::vec2& A, const glm::vec2& B) {
	return A.x * B.y - A.y * B.x;
}

bool isInCircumcircle(const glm::vec2& point, const std::vector<glm::vec2>& points, const glm::uvec3& triangle)
{
	glm::vec3 _A = glm::vec3(points[triangle.x],0);
	glm::vec3 _B = glm::vec3(points[triangle.y],0);
	glm::vec3 _C = glm::vec3(points[triangle.z],0);
	glm::vec3 _D = glm::vec3(point,0);


	glm::mat4 m({
		_A.x,_A.y,std::pow(_A.x,2) + std::pow(_A.y,2),1,
		_D.x,_D.y,std::pow(_D.x,2) + std::pow(_D.y,2),1,
		_C.x,_C.y,std::pow(_C.x,2) + std::pow(_C.y,2),1,
		_B.x,_B.y,std::pow(_B.x,2) + std::pow(_B.y,2),1
		});


	return glm::determinant(m) > 0;
}

void edgeFlip(std::vector<glm::uvec2>& edges, std::vector<glm::uvec3>& triangles, const std::vector<glm::vec2>& points) 
{
	std::set<unsigned int> uniqueIndices;
	for (int i = 0; i < triangles.size(); i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{			
			uniqueIndices.clear();
			uniqueIndices.insert(triangles[i].x);
			uniqueIndices.insert(triangles[i].y);
			uniqueIndices.insert(triangles[i].z);
			uniqueIndices.insert(triangles[j].x);
			uniqueIndices.insert(triangles[j].y);
			uniqueIndices.insert(triangles[j].z);
			if (uniqueIndices.size() == 4)
			{
				uniqueIndices.clear();
				uniqueIndices.insert(triangles[i].x);
				uniqueIndices.insert(triangles[i].y);
				uniqueIndices.insert(triangles[i].z);
				int fp = 0;				
				for (int k = 0; k < 3; k++)
				{
					if (std::find(uniqueIndices.begin(), uniqueIndices.end(), triangles[j][k]) == uniqueIndices.end())
					{
						fp = triangles[j][k];
						k = 3;
					}
				}
				if (isInCircumcircle(points[fp], points, triangles[i]))
				{
					uniqueIndices.clear();
					uniqueIndices.insert(triangles[j].x);
					uniqueIndices.insert(triangles[j].y);
					uniqueIndices.insert(triangles[j].z);
					int fp2 = 0;
					for (int k = 0; k < 3; k++)
					{
						if (std::find(uniqueIndices.begin(), uniqueIndices.end(), triangles[i][k]) == uniqueIndices.end())
						{
							fp2 = triangles[i][k];
							k = 3;
						}
					}
					uniqueIndices.erase(fp);					
					auto it = uniqueIndices.begin();
					triangles[i] = glm::uvec3(fp, fp2,*it);
					it++;
					triangles[j] = glm::uvec3(fp2, fp,*it);
				}
			}
		}
	}
}


void edgeFlip2(std::vector<glm::uvec2>& edges, std::vector<glm::uvec3>& triangles, const std::vector<glm::vec2>& points)
{
	std::set<unsigned int> uniqueIndices;
	for (int i = triangles.size()-1; i < triangles.size(); i++)
	{
		for (int j = 0; j < triangles.size(); j++)
		{
			uniqueIndices.clear();
			uniqueIndices.insert(triangles[i].x);
			uniqueIndices.insert(triangles[i].y);
			uniqueIndices.insert(triangles[i].z);
			uniqueIndices.insert(triangles[j].x);
			uniqueIndices.insert(triangles[j].y);
			uniqueIndices.insert(triangles[j].z);
			if (uniqueIndices.size() == 4)
			{
				uniqueIndices.clear();
				uniqueIndices.insert(triangles[i].x);
				uniqueIndices.insert(triangles[i].y);
				uniqueIndices.insert(triangles[i].z);
				int fp = 0;
				for (int k = 0; k < 3; k++)
				{
					if (std::find(uniqueIndices.begin(), uniqueIndices.end(), triangles[j][k]) == uniqueIndices.end())
					{
						fp = triangles[j][k];
						k = 3;
					}
				}
				if (isInCircumcircle(points[fp], points, triangles[i]))
				{
					uniqueIndices.clear();
					uniqueIndices.insert(triangles[j].x);
					uniqueIndices.insert(triangles[j].y);
					uniqueIndices.insert(triangles[j].z);
					int fp2 = 0;
					for (int k = 0; k < 3; k++)
					{
						if (std::find(uniqueIndices.begin(), uniqueIndices.end(), triangles[i][k]) == uniqueIndices.end())
						{
							fp2 = triangles[i][k];
							k = 3;
						}
					}
					uniqueIndices.erase(fp);
					auto it = uniqueIndices.begin();
					triangles[i] = glm::uvec3(fp, fp2, *it);
					it++;
					triangles[j] = glm::uvec3(fp2, fp, *it);
					break;
				}
			}
		}
	}
}

bool almost_equal(float x, float y, int ulp = 2)
{
	return fabsf(x - y) <= std::numeric_limits<float>::epsilon() * fabsf(x + y) * static_cast<float>(ulp)
		|| fabsf(x - y) < std::numeric_limits<float>::min();
}

bool almost_equal(const glm::vec2& v1, const glm::vec2& v2)
{
	return almost_equal(v1.x, v2.x) && almost_equal(v1.y, v2.y);
}

bool almost_equal(const glm::vec3& v1, const glm::vec3& v2)
{
	return almost_equal(v1.x, v2.x) && almost_equal(v1.y, v2.y) && almost_equal(v1.z, v2.z);
}

struct Triangle {
	glm::vec2 * x;
	glm::vec2 * y;
	glm::vec2 * z;
	bool isBad;	

	Triangle(glm::vec2* tx, glm::vec2* ty, glm::vec2* tz)
	{
		x = tx;
		y = ty;
		z = tz;
		isBad = false;
	}

	float norm2(const glm::vec2 &p) const
	{
		return p.x * p.x + p.y * p.y;
	}

	float dist2(const glm::vec2 &v1, const glm::vec2& v2) const
	{
		float dx = v1.x - v2.x;
		float dy = v1.y - v2.y;
		return dx * dx + dy * dy;
	}

	bool containsVertex(const glm::vec2& v) const
	{
		// return p1 == v || p2 == v || p3 == v;
		return almost_equal(*x, v) || almost_equal(*y, v) || almost_equal(*z, v);
	}

	bool circumCircleContains(const glm::vec2& v) const
	{
		float ab = norm2(*x);
		float cd = norm2(*y);
		float ef = norm2(*z);

		float ax = x->x;
		float ay = x->y;
		float bx = y->x;
		float by = y->y;
		float cx = z->x;
		float cy = z->y;

		float circum_x = (ab * (cy - by) + cd * (ay - cy) + ef * (by - ay)) / (ax * (cy - by) + bx * (ay - cy) + cx * (by - ay));
		float circum_y = (ab * (cx - bx) + cd * (ax - cx) + ef * (bx - ax)) / (ay * (cx - bx) + by * (ax - cx) + cy * (bx - ax));

		glm::vec2 circum(circum_x / 2, circum_y / 2);
		float circum_radius = dist2(*x,circum);
		float dist = dist2(v,circum);
		return dist <= circum_radius;
	}
};

struct Edge
{
	glm::vec2* x;
	glm::vec2* y;
	bool isBad;
	Edge(glm::vec2* ex, glm::vec2* ey)
	{
		x = ex;
		y = ey;
		isBad = false;
	}
};

struct Vec3PtrHash {
	size_t operator()(const glm::vec3* ptr) const {
		return reinterpret_cast<size_t>(ptr);
	}
};

struct Vec3PtrEqual {
	bool operator()(const glm::vec3* lhs, const glm::vec3* rhs) const {
		return almost_equal(*lhs, *rhs);
	}
};

struct Edge3D;
class Sommet3D
{
public:
	glm::vec3* s;
	int color = ROUGE;
	std::vector<Edge3D*> edges;
	Sommet3D(glm::vec3* es)
	{
		s = es;
	}
};

struct Face;
struct Edge3D
{
	std::vector<Sommet3D*> points;
	Sommet3D* x;
	Sommet3D* y;
	int color = ROUGE;
	std::vector<Face*> faces;

	Edge3D(Sommet3D* ex, Sommet3D* ey)
	{
		x = ex;
		y = ey;
		x->edges.push_back(this);
		y->edges.push_back(this);
		points.push_back(x);
		points.push_back(y);
	}
};

struct Face
{
	std::vector<Edge3D*> edges;
	Edge3D* x;
	Edge3D* y;
	Edge3D* z;
	std::vector<Sommet3D*> list_Point;
	glm::vec3 centre;
	bool inactive = false;
	int color = ROUGE;

	std::vector<Sommet3D*> GetPoints()
	{
		std::vector<Sommet3D*> lp;
		for (int k = 0; k < edges.size(); k++)
		{
			for (int l = 0; l < 2; l++)
			{
				bool contain = false;
				for (int i = 0; i < lp.size() && !contain; i++)
				{
					if (almost_equal(*edges[k]->points[l]->s, *lp[i]->s))
					{
						contain = true;
					}
				}
				if (!contain)
				{
					lp.push_back(edges[k]->points[l]);
				}
			}
		}
		return lp;
	}
	Face(Edge3D* ex, Edge3D* ey, Edge3D* ez)
	{
		x = ex;
		y = ey;
		z = ez;
		x->faces.push_back(this);
		y->faces.push_back(this);
		z->faces.push_back(this);
		edges.push_back(x);
		edges.push_back(y);
		edges.push_back(z);
		list_Point = GetPoints();
		std::cout << "Face(si > 3 error) =" << list_Point.size() << std::endl;
		centre = (*list_Point[0]->s + *list_Point[1]->s + *list_Point[2]->s) / 3.0f;
	}

	glm::vec3 Normal()
	{
		glm::vec3 edge1 = *list_Point[2]->s - *list_Point[0]->s;
		glm::vec3 edge2 = *list_Point[2]->s - *list_Point[1]->s;
		return glm::normalize(glm::cross(edge1, edge2));
	}

	glm::vec3 Normal(glm::vec3 bary)
	{
		glm::vec3 edge1 = *list_Point[2]->s - *list_Point[0]->s;
		glm::vec3 edge2 = *list_Point[2]->s - *list_Point[1]->s;
		glm::vec3 normal = glm::normalize(glm::cross(edge1, edge2));
		if (glm::dot(glm::normalize(centre-bary), normal) < 0)
		{
			normal = -normal;
		}
		return normal;
	}
	
};

struct Tetra
{
	std::vector<Face*> faces;

	glm::vec3 bary;
	Tetra(Face* ex, Face* ey, Face* ez, Face* ew)
	{
		faces.push_back(ex);
		faces.push_back(ey);
		faces.push_back(ez);
		faces.push_back(ew);

		std::vector<Sommet3D*> lp;
		for (int m = 0; m < faces.size(); m++)
		{
			for (int k = 0; k < faces[m]->list_Point.size(); k++)
			{

				bool contain = false;
				for (int i = 0; i < lp.size() && !contain; i++)
				{
					if (almost_equal(*faces[m]->list_Point[k]->s, *lp[i]->s))
					{
						contain = true;
					}
				}
				if (!contain)
				{
					lp.push_back(faces[m]->list_Point[k]);
				}
			}
		}
		bary = (*lp[0]->s + *lp[1]->s + *lp[2]->s + *lp[3]->s)/4.0f;

	}
};

bool almost_equal(const Edge& e1, const Edge& e2)
{
	return	(almost_equal(*e1.x, *e2.x) && almost_equal(*e1.y, *e2.y)) ||
		(almost_equal(*e1.x, *e2.y) && almost_equal(*e1.y, *e2.x));
}

std::vector<Triangle> calculateDelaunayTriangulation(std::vector<Triangle> triangles,std::vector<glm::vec2>& points)
{
	// Créer un super-triangle contenant tous les points
	float minX = points[0].x, minY = points[0].y;
	float maxX = points[0].x, maxY = points[0].y;

	for (size_t i = 1; i < points.size(); ++i) 
	{
		minX = std::min(minX, points[i].x);
		minY = std::min(minY, points[i].y);
		maxX = std::max(maxX, points[i].x);
		maxY = std::max(maxY, points[i].y);
	}

	float dx = maxX - minX;
	float dy = maxY - minY;
	float deltaMax = std::max(dx, dy);
	float midx = (minX + maxX) / 2.0f;
	float midy = (minY + maxY) / 2.0f;

	glm::vec2 p1(midx - 20.0f * deltaMax, midy - deltaMax);
	glm::vec2 p2(midx, midy + 20.0f * deltaMax);
	glm::vec2 p3(midx + 20.0f * deltaMax, midy - deltaMax);
	for (auto& t : triangles)
	{
		t.isBad = false;
	}
	if (triangles.size() == 0)
	{
		triangles.push_back(Triangle(&p1, &p2, &p3));
	}

	for (int i = 0; i < points.size(); i++)
	{
		std::vector<Edge> polygon;

		for (auto& t : triangles)
		{
			if (t.circumCircleContains(points[i]))
			{
				t.isBad = true;
				polygon.push_back(Edge(t.x, t.y));
				polygon.push_back(Edge(t.y, t.z));
				polygon.push_back(Edge(t.z, t.x));
			}
		}

		triangles.erase(std::remove_if(begin(triangles), end(triangles), [](Triangle& t) {
			return t.isBad;
			}), end(triangles));

		for (auto e1 = begin(polygon); e1 != end(polygon); ++e1)
		{
			for (auto e2 = e1 + 1; e2 != end(polygon); ++e2)
			{
				if (almost_equal(*e1, *e2))
				{
					e1->isBad = true;
					e2->isBad = true;
				}
			}
		}

		polygon.erase(std::remove_if(begin(polygon), end(polygon), [](Edge& e) {
			return e.isBad;
			}), end(polygon));

		for (const auto e : polygon)
		{
			triangles.push_back(Triangle(e.x, e.y, &points[i]));
		}
	}
	
	triangles.erase(std::remove_if(begin(triangles), end(triangles), [p1, p2, p3](Triangle& t) {
		return t.containsVertex(p1) || t.containsVertex(p2) || t.containsVertex(p3);
		}), end(triangles));

	return triangles;
}


void addPointToDelaunayTriangulation(glm::vec2& newPoint, std::vector<Triangle>& triangles) {
	std::vector<Edge> polygon;

	for (auto& t : triangles) {
		if (t.circumCircleContains(newPoint)) {
			t.isBad = true;
			polygon.push_back(Edge(t.x, t.y));
			polygon.push_back(Edge(t.y, t.z));
			polygon.push_back(Edge(t.z, t.x));
		}
	}

	triangles.erase(std::remove_if(begin(triangles), end(triangles), [](Triangle& t) {
		return t.isBad;
		}), end(triangles));

	for (auto e1 = begin(polygon); e1 != end(polygon); ++e1) {
		for (auto e2 = e1 + 1; e2 != end(polygon); ++e2) {
			if (almost_equal(*e1, *e2)) {
				e1->isBad = true;
				e2->isBad = true;
			}
		}
	}

	polygon.erase(std::remove_if(begin(polygon), end(polygon), [](Edge& e) {
		return e.isBad;
		}), end(polygon));

	for (auto& e : polygon) {
		triangles.push_back(Triangle(e.x, e.y, &newPoint));
	}

	triangles.erase(std::remove_if(begin(triangles), end(triangles), [&newPoint](Triangle& t) {
		return t.containsVertex(newPoint);
		}), end(triangles));
}

glm::vec2 calculateCircumcenter(const glm::vec2& A, const glm::vec2& B, const glm::vec2& C) {
	float D = 2 * (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
	glm::vec2 circumcenter;

	circumcenter.x = ((A.x * A.x + A.y * A.y) * (B.y - C.y) + (B.x * B.x + B.y * B.y) * (C.y - A.y) +
		(C.x * C.x + C.y * C.y) * (A.y - B.y)) /
		D;

	circumcenter.y = ((A.x * A.x + A.y * A.y) * (C.x - B.x) + (B.x * B.x + B.y * B.y) * (A.x - C.x) +
		(C.x * C.x + C.y * C.y) * (B.x - A.x)) /
		D;

	return circumcenter;
}

struct VoronoiEdge {
	glm::vec2 x;
	glm::vec2 y;
	VoronoiEdge(glm::vec2 ex, glm::vec2 ey)
	{		
		x = ex;
		y = ey;
	}
};

bool point_in_triangle(glm::vec2 p, glm::vec2 a, glm::vec2 b, glm::vec2 c) 
{
	double alpha = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) /
		((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	double beta = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) /
		((b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y));
	double gamma = 1.0 - alpha - beta;

	return (alpha > 0 && beta > 0 && gamma > 0);
}

bool sameEdge(glm::vec2 a, glm::vec2 b, glm::vec2 x, glm::vec2 y)
{
	return almost_equal(a, x) && almost_equal(b, y) || almost_equal(b, x) && almost_equal(a, y);
}

int GetIndexFromTriangle(const std::vector<Triangle>& delaunayTriangles,glm::vec2 * e1, glm::vec2* e2,int ignore)
{	
	for (int i = 0 ; i < delaunayTriangles.size();i++)
	{
		if (ignore == i)
		{
			continue;
		}
		if (
			delaunayTriangles[i].x == e1 && delaunayTriangles[i].y == e2
			|| delaunayTriangles[i].x == e2 && delaunayTriangles[i].y == e1

			|| delaunayTriangles[i].y == e1 && delaunayTriangles[i].z == e2
			|| delaunayTriangles[i].y == e2 && delaunayTriangles[i].z == e1

			|| delaunayTriangles[i].z == e1 && delaunayTriangles[i].x == e2
			|| delaunayTriangles[i].z == e2 && delaunayTriangles[i].x == e1

			)
		{
			return i;
		}
	}
	return -1;
}

std::vector<VoronoiEdge> GeomManager::calculateVoronoiDiagram(const std::vector<Triangle>& delaunayTriangles, std::vector<glm::vec2> points2d)
{	
	std::vector<VoronoiEdge> voronoiEdges;
	std::vector<glm::vec2> circumcenterTriangle;
	for (const auto& triangle : delaunayTriangles) 
	{		
		glm::vec2 circumcenter = calculateCircumcenter(*triangle.x, *triangle.y, *triangle.z);
		circumcenterTriangle.push_back(circumcenter);
	}
	int count = 0;
	for (const auto& triangle : delaunayTriangles)
	{
		glm::vec2 midpointXY = 0.5f * (*triangle.x + *triangle.y);
		glm::vec2 midpointYZ = 0.5f * (*triangle.y + *triangle.z);
		glm::vec2 midpointZX = 0.5f * (*triangle.z + *triangle.x);

		int v = GetIndexFromTriangle(delaunayTriangles, triangle.x, triangle.y, count);
		voronoiEdges.push_back(VoronoiEdge(circumcenterTriangle[count], v == -1 ? midpointXY : circumcenterTriangle[v]));
		if (v == -1)
		{
			glm::vec2 dir = glm::normalize(voronoiEdges[voronoiEdges.size()-1].x - voronoiEdges[voronoiEdges.size() - 1].y);
			glm::vec2 dir1 = glm::normalize(*triangle.z - *triangle.x);
			glm::vec2 dir2 = glm::normalize(*triangle.z - *triangle.y);
			if (glm::dot(dir1, dir2) > 0)
			{
				voronoiEdges[voronoiEdges.size() - 1].x = circumcenterTriangle[count];
				voronoiEdges[voronoiEdges.size() - 1].y = voronoiEdges[voronoiEdges.size() - 1].x - dir * 1000.0f;

			}
			else
			{
				voronoiEdges[voronoiEdges.size() - 1].y += dir * 1000.0f;
			}
		}
		v = GetIndexFromTriangle(delaunayTriangles, triangle.y, triangle.z, count);
		voronoiEdges.push_back(VoronoiEdge(circumcenterTriangle[count], v == -1 ? midpointYZ : circumcenterTriangle[v]));
		if (v == -1)
		{
			glm::vec2 dir = glm::normalize(voronoiEdges[voronoiEdges.size() - 1].x - voronoiEdges[voronoiEdges.size() - 1].y);
			glm::vec2 dir1 = glm::normalize(*triangle.x - *triangle.y);
			glm::vec2 dir2 = glm::normalize(*triangle.x - *triangle.z);
			if (glm::dot(dir1, dir2) > 0)
			{
				voronoiEdges[voronoiEdges.size() - 1].x = circumcenterTriangle[count];
				voronoiEdges[voronoiEdges.size() - 1].y = voronoiEdges[voronoiEdges.size() - 1].x - dir * 1000.0f;

			}
			else
			{
				voronoiEdges[voronoiEdges.size() - 1].y += dir * 1000.0f;
			}
		}
		v = GetIndexFromTriangle(delaunayTriangles, triangle.z, triangle.x, count);
		voronoiEdges.push_back(VoronoiEdge(circumcenterTriangle[count], v == -1 ? midpointZX : circumcenterTriangle[v]));
		if (v == -1)
		{
			glm::vec2 dir = glm::normalize(voronoiEdges[voronoiEdges.size() - 1].x - voronoiEdges[voronoiEdges.size() - 1].y);
			glm::vec2 dir1 = glm::normalize(*triangle.y - *triangle.z);
			glm::vec2 dir2 = glm::normalize(*triangle.y - *triangle.x);
			if (glm::dot(dir1, dir2) > 0)
			{
				voronoiEdges[voronoiEdges.size() - 1].x = circumcenterTriangle[count];
				voronoiEdges[voronoiEdges.size() - 1].y = voronoiEdges[voronoiEdges.size() - 1].x - dir * 1000.0f;

			}
			else
			{
				voronoiEdges[voronoiEdges.size() - 1].y += dir * 1000.0f;
			}
		}
		count++;
	}
	return voronoiEdges;
}

void GeomManager::DelaunayNoyauxTriangulation()
{
	for (int i = 0; i < m_segments.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_segments[i]);
	}
	m_segments.clear();


	if (m_points_clouds.size() == 0)
	{
		return;
	}

	if (m_points_clouds.size() == 1)
	{
		point2d.push_back(glm::vec2(m_points_clouds[0]->getPosition().x, m_points_clouds[0]->getPosition().z));
		return;
	}
	if (m_points_clouds.size() == 2)
	{
		point2d.push_back(glm::vec2(m_points_clouds[1]->getPosition().x, m_points_clouds[1]->getPosition().z));
		m_segments.push_back(createSegment(point2d[0], point2d[1]));
		return;
	}
	if (m_points_clouds.size() == 3)
	{
		point2d.push_back(glm::vec2(m_points_clouds[2]->getPosition().x, m_points_clouds[2]->getPosition().z));
		m_segments.push_back(createSegment(point2d[0], point2d[1]));
		m_segments.push_back(createSegment(point2d[1], point2d[2]));
		m_segments.push_back(createSegment(point2d[2], point2d[0]));
		triangles.push_back(Triangle(&point2d[0], &point2d[1], &point2d[2]));
	}

	glm::vec2 last_point = m_points_clouds[m_points_clouds.size() - 1]->getPosition();

	point2d.push_back(glm::vec2(m_points_clouds[0]->getPosition().x, m_points_clouds[0]->getPosition().z));
	point2d.push_back(glm::vec2(m_points_clouds[1]->getPosition().x, m_points_clouds[1]->getPosition().z));
	point2d.push_back(glm::vec2(m_points_clouds[2]->getPosition().x, m_points_clouds[2]->getPosition().z));
	triangles.push_back(Triangle(&point2d[0], &point2d[1], &point2d[2]));

	if (point_in_triangle(last_point, *triangles[0].x, *triangles[0].y, *triangles[0].z))
	{
		return;
	}
	return;
	for (int i = 0; i < triangles.size(); i++)
	{
		m_segments.push_back(createSegment(*triangles[i].x, *triangles[i].y));
		m_segments.push_back(createSegment(*triangles[i].y, *triangles[i].z));
		m_segments.push_back(createSegment(*triangles[i].z, *triangles[i].x));
	}

	return;

	//Old code
	if (currentTriangle.size() == 0)
	{				
		currentTriangle.clear();
		for (int i = 0; i < m_points_clouds.size(); i++)
		{
			point2d.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
		}
		triangles = calculateDelaunayTriangulation(triangles, point2d);
		currentTriangle.resize(triangles.size());
		for (int i = 0; i < triangles.size(); i++)
		{
			for (int j = 0; j < point2d.size(); j++)
			{				
				if (triangles[i].x == &point2d[j])
				{
					currentTriangle[i].x = j;
				}
				if (triangles[i].y == &point2d[j])
				{
					currentTriangle[i].y = j;
				}
				if (triangles[i].z == &point2d[j])
				{
					currentTriangle[i].z = j;
				}
			}
		}		
	}
	else
	{		
		for (int i = 0; i < m_points_clouds.size(); i++)
		{
			point2d.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
		}		
		for (int i = 0; i < currentTriangle.size(); i++)
		{
			triangles.push_back(Triangle(&point2d[currentTriangle[i].x], &point2d[currentTriangle[i].y], &point2d[currentTriangle[i].z]));
		}
		triangles = calculateDelaunayTriangulation(triangles,point2d);
	}
	for (int i = 0; i < triangles.size(); i++)
	{
		m_segments.push_back(createSegment(*triangles[i].x, *triangles[i].y));
		m_segments.push_back(createSegment(*triangles[i].y, *triangles[i].z));
		m_segments.push_back(createSegment(*triangles[i].z, *triangles[i].x));
	}
}

void GeomManager::DelaunayTriangulation2DTest()
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
					if (SegmentsIntersect(point2D[i] * 100.0f, point2D[j] * 100.0f, point2D[edges[k].x] * 100.0f, point2D[edges[k].y] * 100.0f))
					{
						isIntersect = true;
						break;
					}
				}
			}
			if (!isIntersect)
			{
				new_edges.push_back(glm::uvec2(i, j));
			}
		}
		for (int j = 0; j < new_edges.size(); j++)
		{
			edges.push_back(new_edges[j]);
		}
		
		triangulation.clear();
		std::set<unsigned int> uniqueIndices;
		for (int i = 0; i < edges.size(); i++)
		{
			for (int j = i + 1; j < edges.size(); j++)
			{
				for (int k = j + 1; k < edges.size(); k++)
				{
					uniqueIndices.clear();
					uniqueIndices.insert(edges[i].x);
					uniqueIndices.insert(edges[i].y);
					uniqueIndices.insert(edges[j].x);
					uniqueIndices.insert(edges[j].y);
					uniqueIndices.insert(edges[k].x);
					uniqueIndices.insert(edges[k].y);

					if (uniqueIndices.size() == 3)
					{
						auto it = uniqueIndices.begin();
						glm::uvec3 trig;
						trig.x = *it;
						it++;
						trig.y = *it;
						it++;
						trig.z = *it;
						it++;

						if (orientation(point2D[trig.x], point2D[trig.y], point2D[trig.z]) == 2)
						{
							unsigned int temp = trig.y;
							trig.y = trig.z;
							trig.z = temp;
						}
						triangulation.push_back(trig);
					}
				}
			}
		}				
		if (triangulation.size() >= 2)
		{
			std::cout << edges.size() << " " << triangulation.size() << std::endl;
			edgeFlip2(edges, triangulation, point2D);			
			std::vector<glm::uvec2> nedges;
			for (int r = 0; r < triangulation.size(); r++)
			{
				glm::uvec2 ind = glm::uvec2(triangulation[r].x >= triangulation[r].y ? triangulation[r].x : triangulation[r].y, triangulation[r].x < triangulation[r].y ? triangulation[r].x : triangulation[r].y);
				if (std::find(nedges.begin(), nedges.end(), ind) == nedges.end())
				{
					nedges.push_back(ind);
				}
				ind = glm::uvec2(triangulation[r].y >= triangulation[r].z ? triangulation[r].y : triangulation[r].z, triangulation[r].y < triangulation[r].z ? triangulation[r].y : triangulation[r].z);
				if (std::find(nedges.begin(), nedges.end(), ind) == nedges.end())
				{
					nedges.push_back(ind);
				}
				ind = glm::uvec2(triangulation[r].z >= triangulation[r].x ? triangulation[r].z : triangulation[r].x, triangulation[r].z < triangulation[r].x ? triangulation[r].z : triangulation[r].x);
				if (std::find(nedges.begin(), nedges.end(), ind) == nedges.end())
				{
					nedges.push_back(ind);
				}
				//std::cout << ind.x << " " << ind.y << std::endl;
			}
			if (nedges.size() == edges.size())
			{
				edges = nedges;				
			}
			std::cout << edges.size() << " " << triangulation.size() << std::endl;
		}
	}



	for (int i = 0; i < triangulation.size(); i++)
	{
		m_segments.push_back(createSegment(point2D[triangulation[i].x], point2D[triangulation[i].y]));
		m_segments.push_back(createSegment(point2D[triangulation[i].y], point2D[triangulation[i].z]));
		m_segments.push_back(createSegment(point2D[triangulation[i].z], point2D[triangulation[i].x]));
	}
}

void GeomManager::DelaunayTriangulation2D()
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
	std::vector<Triangle> t;
	std::vector<Triangle> triangles = calculateDelaunayTriangulation(t,point2D);
	for (int i = 0; i < triangles.size(); i++)
	{
		m_segments.push_back(createSegment(*triangles[i].x, *triangles[i].y));
		m_segments.push_back(createSegment(*triangles[i].y, *triangles[i].z));
		m_segments.push_back(createSegment(*triangles[i].z, *triangles[i].x));
	}
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

	std::vector<glm::uvec2> edges;	
	edges.push_back(glm::uvec2(0, 1));
	edges.push_back(glm::uvec2(1, 2));
	edges.push_back(glm::uvec2(2, 0));

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

std::vector<glm::uvec2> GeomManager::grahamScan(std::vector<glm::vec2> point2D)
{
	std::vector<glm::uvec2> result;
	std::vector<int> buffer2D;
	std::vector<int> order;
	for (int i = 0; i < point2D.size(); i++)
	{		
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
		float crossProduct = (point2D[pivot->data].x - point2D[pivot->prev->data].x) * (point2D[pivot->next->data].y - point2D[pivot->data].y) -
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
		result.push_back(glm::uvec2(tmp->data, tmp->next->data));
		tmp = tmp->next;		
	}
	result.push_back(glm::uvec2(tmp->data, tmp->next->data));
	return result;
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

void GeomManager::Voronoi()
{
	std::vector<Triangle> t;
	std::vector<glm::vec2> point2d;
	for (int i = 0; i < m_points_clouds.size(); i++)
	{
		point2d.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
	}
	for (int i = 0; i < currentTriangle.size(); i++)
	{
		t.push_back(Triangle(&point2d[currentTriangle[i].x], &point2d[currentTriangle[i].y], &point2d[currentTriangle[i].z]));
	}
	std::vector<VoronoiEdge> voronoipoint = calculateVoronoiDiagram(t, point2d);
	Debug::Log("voronoi %d", t.size());
	for (int i = 0; i < voronoipoint.size(); i++)
	{
		Model* m = createSegment(voronoipoint[i].x, voronoipoint[i].y, m_segmentMat2);
		m_segments.push_back(m);
		m_segmentsVoronoi.push_back(m);
	}
}

Edge3D* alreadyExist(const std::vector<Edge3D*>& edges, Edge3D* newEdge,bool * exist)
{
	for (int i = 0; i < edges.size(); i++)
	{
		if (almost_equal(*edges[i]->x->s, *newEdge->x->s) && almost_equal(*edges[i]->y->s, *newEdge->y->s) ||
			almost_equal(*edges[i]->y->s, *newEdge->x->s) && almost_equal(*edges[i]->x->s, *newEdge->y->s))
		{
			std::cout << "AAA" << std::endl;
			*exist = true;
			return edges[i];
		}
	}
	*exist = false;
	return newEdge;
}

Face* alreadyExistFace(const std::vector<Face*>& faces, Face* newFace, bool* exist)
{
	for (int i = 0; i < faces.size(); i++)
	{
		if (faces[i]->x == newFace->x && faces[i]->y == newFace->y && faces[i]->z == newFace->z || 
			faces[i]->x == newFace->x && faces[i]->y == newFace->z && faces[i]->z == newFace->y ||

			faces[i]->x == newFace->y && faces[i]->y == newFace->z && faces[i]->z == newFace->x ||
			faces[i]->x == newFace->y && faces[i]->y == newFace->x && faces[i]->z == newFace->z ||

			faces[i]->x == newFace->z && faces[i]->y == newFace->x && faces[i]->z == newFace->y ||
			faces[i]->x == newFace->z && faces[i]->y == newFace->y && faces[i]->z == newFace->x
			)
		{
			std::cout << "BBB" << std::endl;
			*exist = true;
			return faces[i];
		}
	}
	*exist = false;
	return newFace;
}

void GeomManager::ConvexHull3D()
{
	for (int i = 0; i < m_segments.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_segments[i]);
	}
	for (int i = 0; i < m_faceModel.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_faceModel[i]);
	}
	m_faceModel.clear();
	m_segments.clear();
	std::vector<Sommet3D*> sommets;
	std::vector<glm::vec3> fdp;
	for (int i = 0; i < m_points_clouds.size(); i++)
	{
		fdp.push_back(m_points_clouds[i]->getPosition());
	}
	for (int i = 0; i < fdp.size(); i++)
	{
		glm::vec3* pos = fdp.data()+i;
		Sommet3D* s = new Sommet3D(pos);
		sommets.push_back(s);
		std::cout << pos << std::endl;
	}
	std::vector<Edge3D*> edges;
	std::vector<Face*> faces;
	std::vector<Tetra*> tetra;

	if (m_points_clouds.size() > 3)
	{
		edges.push_back(new Edge3D(sommets[0], sommets[1]));
		edges.push_back(new Edge3D(sommets[1], sommets[2]));
		edges.push_back(new Edge3D(sommets[2], sommets[0]));
		faces.push_back(new Face(edges[0], edges[1], edges[2]));
		edges.push_back(new Edge3D(sommets[0], sommets[3]));
		edges.push_back(new Edge3D(sommets[1], sommets[3]));
		edges.push_back(new Edge3D(sommets[2], sommets[3]));
		faces.push_back(new Face(edges[0], edges[3], edges[4]));
		faces.push_back(new Face(edges[2], edges[3], edges[5]));
		faces.push_back(new Face(edges[1], edges[4], edges[5]));
		tetra.push_back(new Tetra(faces[0], faces[1], faces[2], faces[3]));
		
		for (int i = 4; i < m_points_clouds.size(); i++)
		{
			for (int j = 0; j < edges.size(); j++)
			{
				if (edges[j]->color != BLEU)
				{
					edges[j]->color = ROUGE;
				}
			}

			for (int j = 0; j < tetra.size(); j++)
			{
				for (int k = 0; k < tetra[j]->faces.size(); k++)
				{
					if (!tetra[j]->faces[k]->inactive)
					{
						glm::vec3 normal = tetra[j]->faces[k]->Normal(tetra[j]->bary);
						glm::vec3 dir = glm::normalize(tetra[j]->faces[k]->centre - *sommets[i]->s);
						if (glm::dot(normal, dir) < 0)
						{
							tetra[j]->faces[k]->color = BLEU;
						}
					}
				}
			}
			for (int j = 0; j < edges.size(); j++)
			{
				int count = 0;
				for (int k = 0; k < edges[j]->faces.size(); k++)
				{
					if (edges[j]->faces[k]->color == BLEU && !edges[j]->faces[k]->inactive)
					{
						count++;
						edges[j]->color = VIOLET;
						if (count >= 2)
						{
							edges[j]->color = BLEU;
						}
					}
				}
			}

			for (int j = 0; j < faces.size(); j++)
			{
				if (faces[j]->color == BLEU && !faces[j]->inactive)
				{
					bool exist1, exist2, exist3;
					std::vector<Edge3D*> nedges;
					std::vector<Face*> nfaces;
					Edge3D * e1 = alreadyExist(edges, new Edge3D(faces[j]->list_Point[0], sommets[i]), &exist1);
					Edge3D * e2 = alreadyExist(edges, new Edge3D(faces[j]->list_Point[1], sommets[i]), &exist2);
					Edge3D * e3 = alreadyExist(edges, new Edge3D(faces[j]->list_Point[2], sommets[i]), &exist3);
					nedges.push_back(e1);
					nedges.push_back(e2);
					nedges.push_back(e3);
					bool fexist1, fexist2, fexist3;
					for (int s = 0; s < 3; s++)//faces
					{
						for (int k = 0; k < 3; k++)
						{
							for (int l = 0; l < 3; l++)
							{
								if (k != l)
								{
									Face* currentface = new Face(faces[j]->edges[s], nedges[k], nedges[l]);
									if (currentface->list_Point.size() == 3)
									{
										nfaces.push_back(currentface);
										l = 3;
										k = 3;
									}
									else
									{
										delete currentface;
									}
								}
							}
						}
					}

					Face * f1 = alreadyExistFace(faces, nfaces[0], &fexist1);
					Face * f2 = alreadyExistFace(faces, nfaces[1], &fexist2);
					Face * f3 = alreadyExistFace(faces, nfaces[2], &fexist3);
					if (faces[j]->edges[0]->color >= BLEU)
					{
						f1->inactive = true;
					}
					if (faces[j]->edges[1]->color >= BLEU)
					{
						f2->inactive = true;
					}
					if (faces[j]->edges[2]->color >= BLEU)
					{
						f3->inactive = true;
					}
					if (!exist1) { edges.push_back(e1); }
					if (!exist2) { edges.push_back(e2); }
					if (!exist3) { edges.push_back(e3); }
					if (!fexist1) { faces.push_back(f1); }
					if (!fexist2) { faces.push_back(f2); }
					if (!fexist3) { faces.push_back(f3); }

					faces[j]->inactive = true;
					tetra.push_back(new Tetra(faces[j], f1, f2, f3));
					/*for (int s = 0; s < tetra[tetra.size() - 1]->faces.size(); s++)
					{
						glm::vec3 normal = tetra[tetra.size() - 1]->faces[s]->Normal(tetra[tetra.size() - 1]->bary);
						m_segments.push_back(createSegment3D(tetra[tetra.size() - 1]->faces[s]->centre, tetra[tetra.size() - 1]->faces[s]->centre+ normal*5.0f, m_segmentMat2));
					}*/
				}
			}

		}

		for (int j = 0; j < tetra.size(); j++)
		{
			for (int k = 0; k < tetra[j]->faces.size(); k++)
			{
				for (int l = 0; l < tetra[j]->faces[k]->edges.size(); l++)
				{					
					m_segments.push_back(createSegment3D(*tetra[j]->faces[k]->edges[l]->x->s, *tetra[j]->faces[k]->edges[l]->y->s, tetra[j]->faces[k]->edges[l]->color >= BLEU ? m_pointMat : tetra[j]->faces[k]->edges[l]->color == VIOLET ? m_segmentMat3 : m_segmentMat));
					/*if (l == 0)
					{
						m_faceModel.push_back(createTriangle((*tetra[j]->faces[k]->edges[l]->x->s+ *tetra[j]->faces[k]->edges[l]->y->s)/2.0f))
					}*/
				}				
			}
		}

	}

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
	m_faceShape = m_pc.modelManager->allocateBuffer("../Model/triangle.obj");
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

	m_segmentMat2 = m_pc.materialManager->createMaterial();
	m_segmentMat2->setColor(glm::vec3(0.0f, 1.0f, 0.0f));
	m_segmentMat2->setMetallic(0.7f);
	m_segmentMat2->setRoughness(0.15f);
	m_segmentMat2->setPipeline(gp_unlit);

	m_segmentMat3 = m_pc.materialManager->createMaterial();
	m_segmentMat3->setColor(glm::vec3(0.53f, 0.42f, 1.0f));
	m_segmentMat3->setMetallic(0.7f);
	m_segmentMat3->setRoughness(0.15f);
	m_segmentMat3->setPipeline(gp_unlit);
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

	if (m_convexHull3DChange)
	{
		m_convexHull3DChange = false;
		m_convexHull3D = !m_convexHull3D;
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
		m_segmentsVoronoi.clear();
		m_points_clouds.clear();
		m_points_light_clouds.clear();
		currentTriangle.clear();
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
				if (m_convexHull3D)
				{
					m->setPosition(glm::vec3(pos.x, m_height, pos.y));
				}
				else
				{
					m->setPosition(glm::vec3(pos.x, 0.0f, pos.y));
				}
				m_points_clouds.push_back(m);
				m_cloudButton = false;
				//PointLight* pl = m_pc.lightManager->createPointLight(glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, 1.0f));
				//pl->setRange(1.0f);
				//pl->setPosition(m->getPosition());
				//m_points_light_clouds.push_back(pl);
				if (m_points_clouds.size() > 2)
				{
					if (m_convexHull) 
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
					else if(m_convexHull3D)
					{
						ConvexHull3D();
					}
					else
					{
						Triangulation2D();
					}

					
				}
				if (m_delaunayNoyaux)
				{
					DelaunayNoyauxTriangulation();
				}
			}
		}
		else
		{
			m_cloudButton = true;
		}
	}

	if (m_moveVoronoi)
	{
		if (m_points_clouds.size() > 2)
		{
			std::vector<Triangle> t;
			std::vector<glm::vec2> point2d;
			for (int i = 0; i < m_points_clouds.size(); i++)
			{
				point2d.push_back(glm::vec2(m_points_clouds[i]->getPosition().x, m_points_clouds[i]->getPosition().z));
			}
			for (int i = 0; i < point2d.size(); i++)
			{
				point2d[i].x += glm::sin(m_pc.time->getTime()+ point2d[i].y)*0.5f;
				//point2d[i].y += glm::sin(m_pc.time->getTime());
			}
			for (int i = 0; i < currentTriangle.size(); i++)
			{
				t.push_back(Triangle(&point2d[currentTriangle[i].x], &point2d[currentTriangle[i].y], &point2d[currentTriangle[i].z]));
			}
			std::vector<VoronoiEdge> voronoipoint = calculateVoronoiDiagram(t, point2d);
			for (int i = 0; i < voronoipoint.size(); i++)
			{				
				updateSegment(voronoipoint[i].x, voronoipoint[i].y, m_segmentsVoronoi[i]);
			}
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
			if (m_convexHull3D)
			{
				ImGui::SliderFloat("Height", &m_height, 0.0f, 5.0f);
			}
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
				m_convexHull3DCheck = false;
				m_convexHull3DChange = false;
				m_convexHull3D = false;
				m_delaunayNoyaux = false;
				m_delaunayNoyauxCheck = false;
			}

			if (ImGui::Checkbox("Create 3D Convex Hull", &m_convexHull3DCheck))
			{
				m_convexHullChange = false;
				m_convexHullCheck = false;
				m_convexHull3DChange = true;
				m_convexHull = false;
				m_delaunayNoyaux = false;
				m_delaunayNoyauxCheck = false;
			}

			
			if (ImGui::Checkbox("Create delaunay Noyaux", &m_delaunayNoyauxCheck))
			{
				m_delaunayNoyaux = true;
				m_convexHullChange = false;
				m_convexHullCheck = false;
				m_convexHull = false;
				m_convexHull3DCheck = false;
				m_convexHull3DChange = false;
				m_convexHull3D = false;
			}			
			if (ImGui::Checkbox("Move Voronoi", &m_moveVoronoi))
			{

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
			if (ImGui::Button("Delaunay Triangulation 2D"))
			{
				if (m_points_clouds.size() > 2)
				{
					DelaunayTriangulation2D();
				}
			}
			/*if (ImGui::Button("Delaunay Triangulation New 2D"))
			{
				if (m_points_clouds.size() > 2)
				{
					DelaunayTriangulation2DTest();
				}
			}*/
			if (ImGui::Button("Delaunay Triangulation noyaux"))
			{
				if (/*m_points_clouds.size() > 2*/true)
				{
					DelaunayNoyauxTriangulation();
				}
			}

			if (ImGui::Button("Voronoi Diagram"))
			{
				if (m_points_clouds.size() > 2)
				{
					Voronoi();
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