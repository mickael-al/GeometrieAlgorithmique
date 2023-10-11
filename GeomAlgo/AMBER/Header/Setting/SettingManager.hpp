#ifndef __ENGINE_SETTING_MANAGER__
#define __ENGINE_SETTING_MANAGER__

#include "Initializer.hpp"
#include "Debug.hpp"
#include "SettingInfo.hpp"

namespace Ge
{
    class SettingManager
    {
    public:
        double getFps() const;
        void setFps(double fps);
        const char *getName() const;
        void setName(const char *name);
        double getWindowHeight() const;
        void setWindowHeight(double height);
        double getWindowWidth() const;
        void setWindowWidth(double Width);
        glm::vec3 getGravity() const;
        void setGravity(glm::vec3 gravity);
        void setVersion(Version v);
        Version getVersion() const;
        void setClearColor(glm::vec4 color);
		glm::vec4 getClearColor() const;
        void setGamma(float gamma);
        float getGamma() const;
		void setIconPath(const char * path);
		const char * getIconPath() const;

    private:
        SettingInfo m_settingInfo;
    };
}

#endif //__ENGINE_SETTING_MANAGER__