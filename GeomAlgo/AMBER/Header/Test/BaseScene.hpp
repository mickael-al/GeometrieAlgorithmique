#ifndef __BASE_SCENE__
#define __BASE_SCENE__

#include "Scene.hpp"
#include "GameEngine.hpp"
#include "BoneManager.hpp"
#include <vector>

class BaseScene : public Scene
{
public:
	void load();
	void unload();
private:
	ptrClass m_pc;
	BoneManager* m_em;
};

#endif //!__BASE_SCENE__