#include "InputManager.hpp"

namespace Ge
{
	bool InputManager::initialize(const VulkanMisc * vM)
	{
		Debug::Info("Initialisation de l'InputManager");
		m_window = vM->str_VulkanDeviceMisc->str_window;
		glfwGetCursorPos(m_window, &m_lastMousePos.x, &m_lastMousePos.y);
		return true;
	}

	void InputManager::updateAxis()
	{
		glfwGetCursorPos(m_window, &m_Mpos.x, &m_Mpos.y);
		m_axisMouse = m_Mpos - m_lastMousePos;	
		m_lastMousePos = m_Mpos;
	}

	void InputManager::release()
	{
		Debug::Info("Liberation de l'InputManager");
	}

	const char * InputManager::getGamepadName(int jid) const
	{
		return glfwGetGamepadName(jid);
	}

	bool InputManager::getGamepadState(int jid, int key) const
	{
		GLFWgamepadstate state;
		glfwGetGamepadState(jid, &state);
		if (key > 15)
		{
			return false;
		}
		return state.buttons[key] == GLFW_PRESS;
	}

	float InputManager::getGamepadAxis(int jid, int indice) const
	{
		GLFWgamepadstate state;
		glfwGetGamepadState(jid, &state);
		if (indice > 5)
		{
			Debug::Error("l'indice de GetGamepadAxis n'existe pas");
			return -2;
		}
		return state.axes[indice];
	}

	bool InputManager::getKey(int key) const
	{
		return glfwGetKey(m_window, key) == GLFW_PRESS;
	}

	bool InputManager::getKeyUp(int key)
	{
		if (mapReleasedInput.find(key) == mapReleasedInput.end())
		{
			mapReleasedInput[key] = false;
		}
		if (glfwGetKey(m_window, key) == GLFW_PRESS)
		{
			mapReleasedInput[key] = true;
			return false;
		}
		else
		{
			if (mapReleasedInput[key])
			{
				mapReleasedInput[key] = false;
				return true;
			}

		}
		return false;
	}

	bool InputManager::getKeyDown(int key)
	{
		if (mapPressedInput.find(key) == mapPressedInput.end())
		{
			mapPressedInput[key] = false;
		}
		if (glfwGetKey(m_window, key) == GLFW_PRESS)
		{
			if (!mapPressedInput[key])
			{
				mapPressedInput[key] = true;
				return true;
			}
		}
		else
		{
			mapPressedInput[key] = false;
		}
		return false;
	}

	bool InputManager::getMouse(int key) const
	{
		return glfwGetMouseButton(m_window, key) == GLFW_PRESS;
	}
	
	void InputManager::HideMouse(bool state)
	{
		glfwSetInputMode(m_window, GLFW_CURSOR, state ? GLFW_CURSOR_DISABLED : GLFW_CURSOR_NORMAL);
	}

	double InputManager::axisMouseX()
	{
		return m_axisMouse.x;
	}

	double InputManager::axisMouseY()
	{
		return m_axisMouse.y;
	}

	double InputManager::getMousePosX() const
	{
		return m_Mpos.x;
	}

	double InputManager::getMousePosY() const
	{
		return m_Mpos.y;
	}
}