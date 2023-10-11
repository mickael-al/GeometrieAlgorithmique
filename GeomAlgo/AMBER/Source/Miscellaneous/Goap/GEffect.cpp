#include "GEffect.hpp"
#include "EnumModifier.hpp"
#include "GNode.hpp"

GEffect::GEffect(const std::string & ressource, const Modifier mod,const float value)
{
	this->m_ressource = ressource;
	this->m_modifier = mod;
	this->m_value = value;
}

GEffect::GEffect(const GEffect* const effect)
{
	this->m_ressource = effect->m_ressource;
	this->m_modifier = effect->m_modifier;
	this->m_value = effect->m_value;
}

GEffect::~GEffect()
{
	this->m_ressource.clear();
}

const std::string & GEffect::GetRessource() const
{
	return this->m_ressource;
}

const Modifier & GEffect::GetModifier() const
{
	return this->m_modifier;
}

void GEffect::ExecuteEffect(GNode* const node) const
{	
	node->ressource[m_ressource] += (m_modifier == Modifier::ADD ? this->m_value : -this->m_value);
	this->Execute();
}

void GEffect::Execute() const
{ 

}