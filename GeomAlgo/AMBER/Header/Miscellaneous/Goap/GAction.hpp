#ifndef __ACTION__
#define __ACTION__

#include <iostream>
#include <vector>

class GEffect;
class GPreCondition;

class GAction final
{
public:
	GAction(const std::string & name,const unsigned short cost = 1);
	GAction(const GAction* const action);
	~GAction();
	void SetName(const std::string & name);
	const std::string & GetName() const;
	unsigned short GetCost() const;
	void SetCost(const unsigned short couts);
	void AddPreCondition(GPreCondition* preCondition);
	const std::vector<GPreCondition*> & GetPreConditions() const;
	void AddEffect(const GEffect* effect);
	const std::vector<const GEffect*> & GetEffects() const;
private:
	std::string m_name;
	std::vector<GPreCondition*> m_preConditions;
	std::vector<const GEffect*> m_effects;
	unsigned short m_cost;
};

#endif // !__ACTION__
