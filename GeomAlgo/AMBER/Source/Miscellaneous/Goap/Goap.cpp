#include "Goap.hpp"
#include "CompareCost.hpp"
#include "EnumModifier.hpp"
#include "GAction.hpp"
#include "GPreCondition.hpp"
#include "GEffect.hpp"
#include "Debug.hpp"
#include <queue>
#include <algorithm>

using namespace Ge;

Goap::Goap(std::vector<GAction*> & allActions, std::vector<GAction*> & objectifs)
{
	m_allActions = allActions;
	m_objectifs = objectifs;
	for (int i = 0; i < m_objectifs.size(); i++)
	{
		bool isExist = false;
		for (int j = 0; j < m_allActions.size() && !isExist; j++)
		{
			if (m_allActions[j] == m_objectifs[i])
			{
				isExist = true;
			}
		}
		if (!isExist)
		{
			Debug::Error("Objectif n'est pas dans action");
		}
	}	
	CalculateResolver();
}

Goap::~Goap()
{

}

void Goap::CalculateResolver() const
{
	std::vector<GPreCondition*> precondition;
	std::vector<GPreCondition*> local_precondition;
	for (int i = 0; i < m_allActions.size(); i++)
	{
		local_precondition = m_allActions[i]->GetPreConditions();
		for (int j = 0; j < local_precondition.size(); j++)
		{
			if (std::find(precondition.begin(), precondition.end(), local_precondition[j]) == precondition.end())
			{
				precondition.push_back(local_precondition[j]);
			}
		}
	}
	std::vector<const GEffect*> local_effects;
	Condition cond;
	for (int p = 0; p < precondition.size(); p++)
	{
		std::vector<const GAction*> effect_to_actions;		
		for (int i = 0; i < m_allActions.size(); i++)
		{
			local_effects = m_allActions[i]->GetEffects();
			for (int j = 0; j < local_effects.size(); j++)
			{
				if (precondition[p]->GetRessource() == local_effects[j]->GetRessource() && std::find(effect_to_actions.begin(), effect_to_actions.end(), m_allActions[i]) == effect_to_actions.end())
				{
					cond = precondition[p]->GetCondition();

					if (((cond == Condition::INF || cond == Condition::INF_EQUALS) && local_effects[j]->GetModifier() == Modifier::SUB) ||
						((cond == Condition::SUP || cond == Condition::SUP_EQUALS) && local_effects[j]->GetModifier() == Modifier::ADD))
					{
						effect_to_actions.push_back(m_allActions[i]);
					}

				}
			}
		}
		precondition[p]->SetResolver(effect_to_actions);
	}

	/*for (int i = 0; i < precondition.size(); i++)
	{
		std::string allAction;
		std::vector<const GAction*> resolver = precondition[i]->GetResolver();
		for (int j = 0; j < resolver.size(); j++)
		{
			allAction += resolver[j]->GetName() + " ";
		}
		std::cout << "Precondition : " << precondition[i]->GetRessource() << " -> " << allAction << std::endl;
	}*/
}

std::vector<const GAction*> Goap::Resolve()
{
	if (m_objectifs.empty()) 
	{
		Debug::Error("Aucun objectif à résoudre.");
		return std::vector<const GAction*>();
	}
	std::priority_queue<GNode*, std::vector<GNode*>, CompareCost> openGNode;
	const GAction* current = nullptr;
	GNode* currentGNode = nullptr;
	std::vector<GPreCondition*> conditions;
	std::vector<const GAction*> actions_resolve;
	std::vector<const GEffect*> action_effects;
	GNode* firstwa = new GNode();
	
	std::vector<GNode*> worlActionRessource;
	firstwa->openNode.push_back(m_objectifs[0]);	
	firstwa->action = m_objectifs[0];
	firstwa->ressource = m_ressource;
	firstwa->cost = 0.0f;
	
	std::vector<const GAction*> action_allready_add;
	openGNode.push(firstwa);
	GNode* finalAction = nullptr;
	GNode* wa = nullptr;
	bool popCurrentWorld = false;
	bool checkCond = false;
	bool isInOpenNode = false;
	int curindex, i, j, s;
	while (!openGNode.empty() && finalAction == nullptr)
	{		
		currentGNode = openGNode.top();
		popCurrentWorld = false;		
		action_allready_add.clear();
		for (curindex = 0; curindex < currentGNode->openNode.size(); curindex++)
		{
			current = currentGNode->openNode[curindex];
			conditions = current->GetPreConditions();
			checkCond = true;
			for (i = 0; i < conditions.size() && checkCond; i++)
			{
				if (!conditions[i]->CheckPreCondition(currentGNode))
				{
					checkCond = false;
					actions_resolve = conditions[i]->GetResolver();
					for (j = 0; j < actions_resolve.size(); j++)
					{
						isInOpenNode = false;
						for (s = 0; s < currentGNode->openNode.size() && !isInOpenNode; s++)
						{
							if (currentGNode->openNode[s] == actions_resolve[j])
							{
								isInOpenNode = true;
							}
						}
						if (!isInOpenNode)
						{
							currentGNode->openNode.push_back(actions_resolve[j]);
						}
					}							
				}
			}
			if (checkCond)
			{
				if (!popCurrentWorld)
				{
					openGNode.pop();
					popCurrentWorld = true;
				}
				currentGNode->openNode.erase(currentGNode->openNode.begin() + curindex);
				curindex--;
				
				isInOpenNode = false;
				for (s = 0; s < action_allready_add.size() && !isInOpenNode; s++)
				{
					if (action_allready_add[s] == current)
					{
						isInOpenNode = true;
					}
				}
				if (!isInOpenNode)
				{
					wa = new GNode();
					wa->action = current;
					wa->openNode = currentGNode->openNode;
					wa->parent = currentGNode;
					wa->ressource = currentGNode->ressource;
					wa->cost += current->GetCost() + currentGNode->cost;
					action_effects = current->GetEffects();
					action_allready_add.push_back(current);
					for (int i = 0; i < action_effects.size(); i++)
					{
						action_effects[i]->ExecuteEffect(wa);
					}
					openGNode.push(wa);	
					worlActionRessource.push_back(wa);
					if (firstwa->action == wa->action)
					{
						finalAction = wa;
					}
				}
			}		
		}		
	}

	if (finalAction != nullptr)
	{
		m_ressource = finalAction->ressource;
	}

	std::vector<const GAction*> path;
	while (finalAction != nullptr)
	{
		if (finalAction->parent == nullptr)
		{
			break;
		}
		path.push_back(finalAction->action);
		finalAction = finalAction->parent;
	}

	for (int i = 0; i < worlActionRessource.size(); i++)
	{
		delete worlActionRessource[i];
	}
	worlActionRessource.clear();
	
	if (path.empty())
	{
		Debug::Error("Aucun chemin trouvé.");
	}
	m_objectifs.erase(m_objectifs.begin());	
	return path;
}