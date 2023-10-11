#include "PathFindingScene.hpp"

PathFindingScene::PathFindingScene(glm::vec3 position, glm::vec3 zoneSize, glm::vec3 pointCount, std::string path, float liasonPercent)
{
	m_pc = GameEngine::getPtrClass();
	m_position = position;
	m_zoneSize = zoneSize;
	m_pointCount = pointCount;
	m_liasonPercent = liasonPercent;
	m_finder = nullptr;
	m_boundsMax = glm::vec3(0);
	m_boundsMin = glm::vec3(0);
	m_vec3Astars = nullptr;
	m_path = path;
}

PathFindingScene::~PathFindingScene()
{
	if (m_finder != nullptr)
	{
		delete m_finder;
		delete[] m_vec3Astars;
	}
	m_points.clear();
	m_neighbors.clear();	
}

uint64_t morton3D_64(uint64_t x, uint64_t y, uint64_t z) 
{
	x = (x | (x << 32)) & 0x1F00000000FFFFL;
	x = (x | (x << 16)) & 0x1F0000FF0000FFL;
	x = (x | (x << 8)) & 0x100F00F00F00F00FL;
	x = (x | (x << 4)) & 0x10C30C30C30C30C3L;
	x = (x | (x << 2)) & 0x1249249249249249L;

	y = (y | (y << 32)) & 0x1F00000000FFFFL;
	y = (y | (y << 16)) & 0x1F0000FF0000FFL;
	y = (y | (y << 8)) & 0x100F00F00F00F00FL;
	y = (y | (y << 4)) & 0x10C30C30C30C30C3L;
	y = (y | (y << 2)) & 0x1249249249249249L;

	z = (z | (z << 32)) & 0x1F00000000FFFFL;
	z = (z | (z << 16)) & 0x1F0000FF0000FFL;
	z = (z | (z << 8)) & 0x100F00F00F00F00FL;
	z = (z | (z << 4)) & 0x10C30C30C30C30C3L;
	z = (z | (z << 2)) & 0x1249249249249249L;

	return x | (y << 1) | (z << 2);
}

uint64_t scaleAndCast(float min, float max, float target)
{		
	if (min == max)
	{
		return 0;
	}
	//return static_cast<uint64_t>((target - min) * (2097151.0f / (max - min)));
	//float scaledValue = ((target - min) / (max - min)) * 64.0f;
	//return static_cast<uint64_t>(scaledValue);
	return static_cast<uint64_t>((target - min)*16.0f);
}

uint64_t vec3ToMorton_64(const glm::vec3& vec, const glm::vec3& boundsMin, const glm::vec3& boundsMax)
{
	return morton3D_64(scaleAndCast(boundsMin.x, boundsMax.x,vec.x), scaleAndCast(boundsMin.y, boundsMax.y, vec.y), scaleAndCast(boundsMin.z, boundsMax.z, vec.z));
}

// Fonction pour comparer deux glm::vec3 selon le code de Morton
bool mortonCompare(const glm::vec3& a, const glm::vec3& b, const glm::vec3& boundsMin, const glm::vec3& boundsMax)
{		
	return vec3ToMorton_64(a, boundsMin, boundsMax) < vec3ToMorton_64(b, boundsMin, boundsMax);
}

void calculateBounds(const std::vector<glm::vec3>& points, glm::vec3& boundsMin, glm::vec3& boundsMax) 
{
	if (points.empty()) 
	{
		boundsMin = glm::vec3(0.0f);
		boundsMax = glm::vec3(0.0f);
		return;
	}

	boundsMin = boundsMax = points[0];

	for (const glm::vec3& point : points) 
	{	
		boundsMin.x = std::min(boundsMin.x, point.x);
		boundsMin.y = std::min(boundsMin.y, point.y);
		boundsMin.z = std::min(boundsMin.z, point.z);

		boundsMax.x = std::max(boundsMax.x, point.x);
		boundsMax.y = std::max(boundsMax.y, point.y);
		boundsMax.z = std::max(boundsMax.z, point.z);
	}
}

void PathFindingScene::SortTree()
{
	calculateBounds(m_points, m_boundsMin, m_boundsMax);

	std::sort(m_points.begin(), m_points.end(), [&](const glm::vec3& a, const glm::vec3& b) 
		{
		return mortonCompare(a, b, m_boundsMin,m_boundsMax);
		});
	m_morton.resize(m_points.size());
	for (int i = 0; i < m_points.size(); i++)
	{
		m_morton[i] = vec3ToMorton_64(m_points[i],m_boundsMin,m_boundsMax);
	}		
}


void PathFindingScene::generateMap()
{
	glm::vec3 pointSpacing = m_zoneSize / m_pointCount;	
	float x, y, z;
	glm::vec3 startPosition, endPosition, hitPosition;	
	glm::vec3 zoneSizeDivid = m_zoneSize / 2.0f;	

	m_points.reserve((m_pointCount.x * m_pointCount.y * m_pointCount.z) / 4);
	Debug::Log("PathFinding generate");
	for (x = 0; x < m_zoneSize.x; x += pointSpacing.x)
	{
		for (y = 0; y < m_zoneSize.y; y += pointSpacing.y)
		{
			for (z = 0; z < m_zoneSize.z; z += pointSpacing.z)
			{
				endPosition.x = startPosition.x = x + m_position.x - zoneSizeDivid.x;
				endPosition.y = startPosition.y = y + m_position.y - zoneSizeDivid.y;
				endPosition.z = startPosition.z = z + m_position.z - zoneSizeDivid.z;
				endPosition.y -= 1000.0f;

				if (m_pc.physicsEngine->raycast(startPosition, endPosition, hitPosition))
				{
					m_points.push_back(hitPosition);
				}
			}
		}
	}
	Debug::Log("PathFinding Points: %d", m_points.size());
	filterClosePoints();
	Debug::Log("PathFinding FilterClosePoints : %d", m_points.size());	
	SortTree();
	generateNeighbors();
	FillGraph();
}

void PathFindingScene::debugPoint()
{
	ShapeBuffer* sb = m_pc.modelManager->allocateBuffer("../Model/Cube.obj");
	Materials* material = m_pc.materialManager->createMaterial();
	Materials* materialEdge = m_pc.materialManager->createMaterial();
	material->setColor(glm::vec3(0, 1, 0));
	material->setMetallic(0.0f);
	material->setRoughness(1.0f);

	materialEdge->setColor(glm::vec3(0, 0, 1));
	materialEdge->setMetallic(0.0f);
	materialEdge->setRoughness(1.0f);
	float scale_base = 0.1f;
	for (int i = 0; i < m_points.size(); i++)
	{
		Model* m = m_pc.modelManager->createModel(sb);
		m->setPosition(m_points[i]);
		m->setMaterial(material);
		m->setScale(glm::vec3(scale_base));
	}

	for (int i = 0; i < m_points.size(); i++) 
	{
		for (unsigned int neighborIndex : m_neighbors[i]) 
		{
			glm::vec3 startPos = m_points[i];
			glm::vec3 endPos = m_points[neighborIndex];

			glm::vec3 midPos = (startPos + endPos) * 0.5f;

			glm::vec3 direction = glm::normalize(endPos - startPos);
			float distance = glm::distance(startPos, endPos);
			glm::vec3 scale = glm::vec3(scale_base*0.8f, scale_base*0.8f, distance);

			Model* edgeModel = m_pc.modelManager->createModel(sb);
			edgeModel->setPosition(midPos);
			edgeModel->setRotation(glm::quatLookAt(direction, glm::vec3(0, 1, 0)));
			edgeModel->setScale(scale);
			edgeModel->setMaterial(materialEdge);
		}
	}
}


void PathFindingScene::saveToFile() 
{
	std::ofstream outFile(m_path, std::ios::binary);
	if (!outFile.is_open()) 
	{
		Debug::Error("Failed to open file for writing: %s", m_path.c_str());
		return;
	}

	// Write header
	outFile.write(reinterpret_cast<char*>(&m_position), sizeof(glm::vec3));
	outFile.write(reinterpret_cast<char*>(&m_zoneSize), sizeof(glm::vec3));
	outFile.write(reinterpret_cast<char*>(&m_pointCount), sizeof(glm::vec3));
	outFile.write(reinterpret_cast<char*>(&m_liasonPercent), sizeof(float));

	// Write number of points
	int size_points = m_points.size();
	outFile.write(reinterpret_cast<char*>(&size_points), sizeof(int));

	// Write points
	for (const glm::vec3& point : m_points) 
	{
		outFile.write(reinterpret_cast<const char*>(&point), sizeof(glm::vec3));
	}

	// Write neighbors
	for (const std::vector<unsigned int>& neighbors : m_neighbors) 
	{
		int size_neighbors = neighbors.size();
		outFile.write(reinterpret_cast<const char*>(&size_neighbors), sizeof(int));
		for (unsigned int neighbor : neighbors) 
		{
			outFile.write(reinterpret_cast<const char*>(&neighbor), sizeof(unsigned int));
		}
	}

	outFile.write(reinterpret_cast<char*>(&m_boundsMin), sizeof(glm::vec3));
	outFile.write(reinterpret_cast<char*>(&m_boundsMax), sizeof(glm::vec3));

	// Write number of morton
	int size_morton = m_morton.size();
	outFile.write(reinterpret_cast<char*>(&size_morton), sizeof(int));

	for (const uint64_t& morton : m_morton)
	{
		outFile.write(reinterpret_cast<const char*>(&morton), sizeof(uint64_t));
	}

	outFile.close();
}

void PathFindingScene::loadFromFile() 
{
	std::ifstream inFile(m_path, std::ios::binary);
	if (!inFile.is_open()) 
	{
		Debug::Error("Failed to open file for reading: %s", m_path.c_str());
		return;
	}

	// header
	inFile.read(reinterpret_cast<char*>(&m_position), sizeof(glm::vec3));
	inFile.read(reinterpret_cast<char*>(&m_zoneSize), sizeof(glm::vec3));
	inFile.read(reinterpret_cast<char*>(&m_pointCount), sizeof(glm::vec3));
	inFile.read(reinterpret_cast<char*>(&m_liasonPercent), sizeof(float));

	// Read number of points
	int size_points = 0;
	inFile.read(reinterpret_cast<char*>(&size_points), sizeof(int));

	// Read points
	m_points.resize(size_points);
	for (glm::vec3& point : m_points) 
	{
		inFile.read(reinterpret_cast<char*>(&point), sizeof(glm::vec3));
	}

	// Read neighbors
	m_neighbors.resize(size_points);
	for (std::vector<unsigned int>& neighbors : m_neighbors) 
	{
		int size_neighbors = 0;
		inFile.read(reinterpret_cast<char*>(&size_neighbors), sizeof(int));
		neighbors.resize(size_neighbors);
		for (unsigned int& neighbor : neighbors) 
		{
			inFile.read(reinterpret_cast<char*>(&neighbor), sizeof(unsigned int));
		}
	}

	inFile.read(reinterpret_cast<char*>(&m_boundsMin), sizeof(glm::vec3));
	inFile.read(reinterpret_cast<char*>(&m_boundsMax), sizeof(glm::vec3));

	// Write number of morton
	int size_morton = 0;
	inFile.read(reinterpret_cast<char*>(&size_morton), sizeof(int));

	m_morton.resize(size_morton);

	for (uint64_t& morton : m_morton)
	{
		inFile.read(reinterpret_cast<char*>(&morton), sizeof(uint64_t));
	}

	inFile.close();
	FillGraph();
}

void PathFindingScene::FillGraph()
{
	if (m_finder != nullptr)
	{
		delete m_finder;
		delete[] m_vec3Astars;
	}

	m_finder = new PathFinder<Vec3AStar>();
	m_vec3Astars = new Vec3AStar[m_points.size()];
	for (int i = 0; i < m_points.size(); i++)
	{
		m_vec3Astars[i].position = m_points[i];
	}
	for (int i = 0; i < m_points.size(); i++)
	{
		for (unsigned int neighborIndex : m_neighbors[i])
		{			
			Node * n = (Node*)(m_vec3Astars+i);		
			n->addChild((Node*)(m_vec3Astars+neighborIndex), glm::distance(m_vec3Astars[i].position, m_vec3Astars[neighborIndex].position));			
		}
	}
}

unsigned int PathFindingScene::nearPointIndex(glm::vec3 target)
{
	uint64_t targetMorton = vec3ToMorton_64(target, m_boundsMin, m_boundsMax);

	unsigned int left = 0;
	unsigned int right = static_cast<unsigned int>(m_morton.size()) - 1;
	unsigned int nearestIndex = 0;
	uint64_t nearestDistance = std::numeric_limits<uint64_t>::max();

	while (left <= right) 
	{
		unsigned int middle = (left + right) / 2;
		uint64_t currentMorton = m_morton[middle];

		uint64_t currentDistance = abs(static_cast<int>(currentMorton) - static_cast<int>(targetMorton));

		if (currentDistance < nearestDistance) 
		{
			nearestDistance = currentDistance;
			nearestIndex = middle;
		}

		if (currentMorton < targetMorton)
		{
			left = middle + 1;
		}
		else
		{
			right = middle - 1;
		}
	}

	return nearestIndex;

	/*float nearDistance = FLT_MAX;
	float minDist;
	unsigned int index = 0;
	for (int i = 0; i < m_points.size(); i++)
	{
		minDist = glm::distance(m_points[i], target);
		if (nearDistance > minDist)
		{
			nearDistance = minDist;
			index = i;
		}		
	}
	return index;*/
}


bool PathFindingScene::pathFinding(glm::vec3 * startPosition, glm::vec3 *endPosition, std::vector<Vec3AStar*>* path)
{
	unsigned int startIndex = nearPointIndex(*startPosition);
	unsigned int endIndex = nearPointIndex(*endPosition);

	*startPosition = m_vec3Astars[startIndex].position;
	*endPosition = m_vec3Astars[endIndex].position;

	m_finder->setStart(m_vec3Astars[startIndex]);
	m_finder->setGoal(m_vec3Astars[endIndex]);
	
	if (m_finder->findPath<AStar>((*path)))
	{
		AStar::getInstance().clear();
		return true;
	}
	AStar::getInstance().clear();
	return false;
}

void PathFindingScene::filterClosePoints()
{
	for (int i = 0; i < m_points.size(); i++) 
	{
		for (int j = m_points.size() - 1; j > i; j--) 
		{
			if (glm::distance(m_points[i], m_points[j]) < 0.1f) 
			{
				m_points.erase(m_points.begin() + j);				
			}
		}
	}
}

void PathFindingScene::generateNeighbors()
{
	float maxDistance = (m_zoneSize.x / m_pointCount.x);
	maxDistance += maxDistance * m_liasonPercent;

	m_neighbors.resize(m_points.size());

	for (std::size_t i = 0; i < m_points.size(); ++i) 
	{
		for (std::size_t j = 0; j < m_points.size(); ++j) 
		{
			if (i != j) 
			{
				float distance = glm::distance(m_points[i], m_points[j]);
				if (distance <= maxDistance) 
				{
					m_neighbors[i].push_back(j);
				}
			}
		}
	}
}


void PathFindingScene::start()
{

}

void PathFindingScene::fixedUpdate()
{

}

void PathFindingScene::update()
{
	if (m_pc.inputManager->getKeyDown(GLFW_KEY_P))
	{
		generateMap();
	}
	if (m_pc.inputManager->getKeyDown(GLFW_KEY_O))
	{
		debugPoint();
		Debug::Log("debugPoint");
	}
	if (m_pc.inputManager->getKeyDown(GLFW_KEY_L))
	{
		loadFromFile();
		Debug::Log("load");
	}
	if (m_pc.inputManager->getKeyDown(GLFW_KEY_K))
	{
		saveToFile();
		Debug::Log("Save");
	}
}

void PathFindingScene::stop()
{

}

void PathFindingScene::onGUI()
{

}