#ifndef __BASE_SCENE__
#define __BASE_SCENE__

#include "Scene.hpp"
#include "GameEngine.hpp"
#include "GeomManager.hpp"
#include <vector>

class BaseScene : public Scene
{
public:
	void load();
	void unload();
private:
	ptrClass m_pc;
	GeomManager* m_em;
};

#endif //!__BASE_SCENE__