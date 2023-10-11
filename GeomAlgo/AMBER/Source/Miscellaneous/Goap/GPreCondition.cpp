#include "GPreCondition.hpp"
#include "GAction.hpp"
#include "GNode.hpp"

GPreCondition::GPreCondition(const std::string ressource, const Condition condition, const float value)
{
	m_ressource = ressource;
	m_condition = condition;
	m_value = value;
}

GPreCondition::GPreCondition(const GPreCondition* const preCondition)
{
	m_ressource = preCondition->m_ressource;
	m_condition = preCondition->m_condition;
	m_value = preCondition->m_value;
}

GPreCondition::~GPreCondition()
{
	m_ressource.clear();
}

void GPreCondition::Check() const
{

}

const std::string& GPreCondition::GetRessource() const
{
	return m_ressource;
}

const Condition& GPreCondition::GetCondition() const
{
	return m_condition;
}

bool GPreCondition::CheckPreCondition(GNode* const node) const
{
	Check();
	const float value = node->ressource[m_ressource];
	switch (m_condition)
	{
	case Condition::INF:
		return value < m_value;
	case Condition::SUP:
		return value > m_value;
	case Condition::INF_EQUALS:
		return value <= m_value;
	case Condition::SUP_EQUALS:
		return value >= m_value;
	default:
		return false;
	}
}

void GPreCondition::SetResolver(std::vector<const GAction*>& actions)
{
	m_resolver = actions;
}
const std::vector<const GAction*>& GPreCondition::GetResolver() const
{
	return m_resolver;
}

float GPreCondition::GetValue() const
{
	return m_value;
}