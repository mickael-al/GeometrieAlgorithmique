#include "BaseScene.hpp"

void BaseScene::load()
{
	m_pc = GameEngine::getPtrClass();
	ShapeBuffer* sb = m_pc.modelManager->allocateBuffer("../Model/plane.obj");
	Model* plane = m_pc.modelManager->createModel(sb, "Plane");
	m_em = new BoneManager(plane);
	m_pc.behaviourManager->addBehaviour(m_em);
	m_pc.hud->addBlockUI(m_em);


	Textures* t = m_pc.textureManager->createTexture("../Texture/damier.png", false);
	Materials* mat = m_pc.materialManager->createMaterial();
	mat->setNormalTexture(m_pc.textureManager->createTexture("../Texture/normal_damier.png"));
	mat->setMetallic(0.8f);
	mat->setRoughness(0.2f);
	mat->setAlbedoTexture(t);
	mat->setMetallicTexture(t);
	mat->setTilling(glm::vec2(5.0f, 5.0f));
	plane->setPosition(glm::vec3(0.0f, -1.0f, 0.0f));
	plane->setScale(glm::vec3(10.0f, 1.0f, 10.0f));
	plane->setMaterial(mat);	
}

void BaseScene::unload()
{
	m_pc.behaviourManager->removeBehaviour(m_em);
	m_pc.hud->removeBlockUI(m_em);
}
