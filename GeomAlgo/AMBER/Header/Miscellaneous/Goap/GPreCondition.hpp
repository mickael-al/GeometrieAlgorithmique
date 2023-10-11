#ifndef __G_PRE_CONDITION__
#define __G_PRE_CONDITION__

#include <vector>
#include "EnumCondition.hpp"
#include <string>

class GAction;
class GNode;
enum Condition;

class GPreCondition
{
public:
	GPreCondition(const std::string ressource, const Condition condition, const float value);
	GPreCondition(const GPreCondition* const preCondition);
	~GPreCondition();
	virtual void Check() const;
	const std::string& GetRessource() const;
	const Condition& GetCondition() const;	
	bool CheckPreCondition(GNode* const node) const;
	void SetResolver(std::vector<const GAction*>& actions);
	const std::vector<const GAction*>& GetResolver() const;
	float GetValue() const;
private:
	std::string m_ressource;
	Condition m_condition;
	float m_value;
	std::vector<const GAction*> m_resolver;
};

#endif //!__G_PRE_CONDITION__