#include "GAction.hpp"
#include "GEffect.hpp"

GAction::GAction(const std::string & name, const unsigned short cost)
{
	this->m_name = name;
	this->m_cost = cost;
}
GAction::GAction(const GAction* const action)
{
	this->m_name = action->GetName();
	this->m_preConditions = action->GetPreConditions();
	this->m_effects = action->GetEffects();
}

GAction::~GAction()
{
	this->m_preConditions.clear();
	this->m_effects.clear();
}

void GAction::SetName(const std::string & name)
{
	this->m_name = name;
}

const std::string & GAction::GetName() const
{
	return this->m_name;
}
void GAction::AddPreCondition(GPreCondition* preCondition)
{
	this->m_preConditions.push_back(preCondition);
}
const std::vector<GPreCondition*> & GAction::GetPreConditions() const
{
	return this->m_preConditions;
}

void GAction::AddEffect(const GEffect* effect)
{
	this->m_effects.push_back(effect);
}

const std::vector<const GEffect*> & GAction::GetEffects() const
{
	return this->m_effects;
}

unsigned short GAction::GetCost() const
{
	return this->m_cost;
}

void GAction::SetCost(const unsigned short cost)
{
	this->m_cost = cost;
}