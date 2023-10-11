#ifndef __G_GEffect__
#define __G_GEffect__

#include <string>

class GNode;
enum Modifier;

class GEffect
{
public:
	GEffect(const std::string & ressource,const Modifier mod,const float value);
	GEffect(const GEffect * const effect);
	~GEffect();	
	const std::string & GetRessource() const;
	const Modifier & GetModifier() const;
	void ExecuteEffect(GNode* const node) const;
	virtual void Execute() const;
private:
	std::string m_ressource;	
	Modifier m_modifier;
	float m_value;
};

#endif //!__G_GEffect__
