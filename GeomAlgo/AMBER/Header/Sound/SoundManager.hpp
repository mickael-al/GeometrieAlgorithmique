#ifndef __SOUND_MANAGER__
#define __SOUND_MANAGER__

#include "Debug.hpp"
#include "SoundTypeLoader.hpp"
#include "AL/al.h"
#include "AL/alc.h"
#include "SoundBuffer.hpp"
#include "AudioSource.hpp"

namespace Ge
{	
	class SoundManager
	{
	public:
		SoundBuffer * createBuffer(const char * filepath);
		void releaseBuffer(SoundBuffer * sb);		
		AudioSource* createSource(SoundBuffer* sb, std::string name = "AudioSource");
		void releaseSource(AudioSource* as);
		void setListenerPosition(glm::vec3 position);
		void setListenerVelocity(glm::vec3 velocity);
		void setListenerDirection(glm::vec3 direction);		
		void setDistanceMode(ALenum mode);	
		ALenum getDistanceMode() const;
	private:
		friend class GameEngine;
		bool Initialize();
		void Release();		
	private:
		std::vector<SoundBuffer*> m_buffers;
		std::vector<AudioSource*> m_audios;
		ALCdevice* m_pDevice;
		ALCcontext* m_pContext;
		float m_orientation[6];
		ALenum m_distanceMode = AL_NONE;
	};

};

#endif //!__SOUND_MANAGER__