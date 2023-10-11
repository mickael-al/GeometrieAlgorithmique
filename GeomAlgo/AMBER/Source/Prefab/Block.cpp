#include "Block.hpp"

Block::Block(ptrClass pc,LoadPrefab* lp, glm::vec3 position, glm::vec3 euler, glm::vec3 scale, bool collider)
{
	m_lp = lp;
	m_pc = pc;
	for (int i = 0; i < lp->sb.size(); i++)
	{
		m_model.push_back(m_pc.modelManager->createModel(lp->sb[i]));
		m_model[i]->setPosition(position);
		m_model[i]->setEulerAngles(euler);
		m_model[i]->setScale(scale* lp->meta.mesh_resize_scale);
		m_model[i]->setMaterial(lp->mat[i% lp->mat.size()]);
	}	
	for (int i = 0; i < lp->meta.collider_transform.size() && collider; i++)
	{
		m_collision.push_back(m_pc.physicsEngine->AllocateCollision(new BoxShape(lp->meta.collider_transform[i].collider_size* scale, 1.0f)));
		m_pc.physicsEngine->AddCollision(m_collision[i]);
	}
	UpdateTransform();
}

Model* Block::getMainModel() const
{
	return m_model[0];
}

std::vector<CollisionBody*> Block::getCollision() const
{
	return m_collision;
}

void Block::DebugCollider()
{
	ShapeBuffer * sb = m_pc.modelManager->allocateBuffer("../Model/Cube.obj");
	Materials* material = m_pc.materialManager->createMaterial();
	material->setColor(glm::vec3(0, 1, 0));
	material->setMetallic(0.0f);
	material->setRoughness(1.0f);
	float scale_base = 0.03f;
	for (int i = 0; i < m_collision.size(); i++)
	{
		ColliderTransform ct = m_lp->meta.collider_transform.at(i);
		std::vector<Model*> edge;
		for (int j = 0; j < 12; j++)
		{
			edge.push_back(m_pc.modelManager->createModel(sb));
			edge[j]->setMaterial(material);
			edge[j]->setScale(glm::vec3(j < 4 ? ct.collider_size.x*2 : scale_base, (j >= 4 && j < 8) ? ct.collider_size.y * 2 : scale_base,  (j >= 8) ? ct.collider_size.z * 2 : scale_base));

			glm::vec3 calclocalpos = glm::vec3(
				(j < 4) ? 0.0f : ((j%2==0) ? ct.collider_size.x : -ct.collider_size.x)
				, (j >= 4 && j < 8) ? 0.0f : ((j < 4) ? ((j%2==0) ? ct.collider_size.y : -ct.collider_size.y) : (j == 8 || j == 9) ? ct.collider_size.y : -ct.collider_size.y)
				, (j >= 8) ? 0.0f : (j%4 <= 1 ? ct.collider_size.z : -ct.collider_size.z)
			);
		
			Transform t = m_collision[i]->globalTransform(calclocalpos, glm::vec3(0.0f), glm::vec3(1.0f), ct.collider_pos, glm::quat(glm::radians(ct.collider_eul)), glm::vec3(1.0f));

			Transform t2 = m_collision[i]->globalTransform(t.position, t.rotation, glm::vec3(1.0f), m_model[0]->getPosition(), m_model[0]->getRotation(), m_model[0]->getScale() / m_lp->meta.mesh_resize_scale);
			edge[j]->setPosition(t2.position);
			edge[j]->setRotation(t2.rotation);
		}
	}

}


void Block::UpdateTransform()
{
	for (int i = 0; i < m_collision.size(); i++)
	{
		ColliderTransform ct = m_lp->meta.collider_transform.at(i);
		Transform t = m_collision[i]->globalTransform(ct.collider_pos, ct.collider_eul, glm::vec3(1.0f), m_model[0]->getPosition(), m_model[0]->getRotation(), m_model[0]->getScale() / m_lp->meta.mesh_resize_scale);
		m_collision[i]->setPosition(t.position);
		m_collision[i]->setRotation(t.rotation);		
	}
	for (int i = 1; i < m_model.size(); i++)
	{
		m_model[i]->setPosition(m_model[0]->getPosition());
		m_model[i]->setRotation(m_model[0]->getRotation());
		m_model[i]->setScale(m_model[0]->getScale());
	}
}

Block::~Block() 
{
	for (int i = 0; i < m_model.size(); i++)
	{
		m_pc.modelManager->destroyModel(m_model[i]);
	}
	m_model.clear();
	for (int i = 0; i < m_collision.size(); i++)
	{
		m_pc.physicsEngine->ReleaseCollision(m_collision[i]);
	}
	m_collision.clear();
}