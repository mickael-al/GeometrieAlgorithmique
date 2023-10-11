#ifndef __GOAP_MANAGER__
#define __GOAP_MANAGER__

#include <vector>
#include <unordered_map>
#include <string>
class GAction;

class Goap final
{
public:
	Goap(std::vector<GAction*> & allActions, std::vector<GAction*> & objectifs);
	~Goap();
	std::vector<const GAction*> Resolve();
private:
	void CalculateResolver() const;
private:
	std::vector<GAction*> m_allActions;
	std::vector<GAction*> m_objectifs;
	std::unordered_map<std::string, float> m_ressource;
}; 

#endif //!__GOAP_MANAGER__
