#include "SoundManager.hpp"

namespace Ge
{
	bool SoundManager::Initialize()
	{
		m_pDevice = alcOpenDevice(nullptr);
		if (!m_pDevice)
		{
			Debug::Warn("Echec il y a aucun pheripherique audio!");
			return true;
		}
		m_pContext = alcCreateContext(m_pDevice, nullptr);
		alcMakeContextCurrent(m_pContext);
		if (!m_pContext)
		{
			Debug::Error("Echec de la creation du context audio!");
			return false;
		}
		Debug::INITSUCCESS("SoundManager");
		return true;
	}

	SoundBuffer * SoundManager::createBuffer(const char* filepath)
	{
		uint64_t size;
		uint32_t frequency;
		ALenum format;
		int8_t* buffer = SoundTypeLoader::LoadWavFormat(filepath,&size,&frequency,&format);
		
		if (buffer == nullptr)
		{
			Debug::Error("SoundBuffer nullptr");
			return nullptr;
		}

		SoundBuffer* sb = new SoundBuffer(size, frequency, format,buffer);

		delete[] buffer;
		buffer = nullptr;
		m_buffers.push_back(sb);
		return sb;
	}

	void SoundManager::releaseBuffer(SoundBuffer* sb)
	{
		auto it = std::find(m_buffers.begin(), m_buffers.end(), sb);
		if (it != m_buffers.end()) 
		{
			m_buffers.erase(it);
			delete sb; 
		}
	}

	AudioSource* SoundManager::createSource(SoundBuffer* sb, std::string name)
	{
		AudioSource* audio = new AudioSource(sb, name);
		m_audios.push_back(audio);
		return audio;
	}

	void SoundManager::releaseSource(AudioSource* as)
	{
		auto it = std::find(m_audios.begin(), m_audios.end(), as);
		if (it != m_audios.end())
		{
			m_audios.erase(it);
			delete as;
		}
	}

	void SoundManager::setListenerPosition(glm::vec3 position)
	{
		alListener3f(AL_POSITION, position.x, position.y, position.z);
	}

	void SoundManager::setListenerVelocity(glm::vec3 velocity)
	{
		alListener3f(AL_POSITION, velocity.x, velocity.y, velocity.z);
	}

	void SoundManager::setListenerDirection(glm::vec3 direction)
	{
		m_orientation[0] = direction.x;
		m_orientation[1] = direction.y;
		m_orientation[2] = direction.z;
		m_orientation[3] = 0.0;
		m_orientation[4] = -1.0;
		m_orientation[5] = 0.0;
		alListenerfv(AL_ORIENTATION, m_orientation);
	}

	void SoundManager::setDistanceMode(ALenum mode)
	{
		m_distanceMode = mode;
		alDistanceModel(mode);
	}

	ALenum SoundManager::getDistanceMode() const
	{
		return m_distanceMode;
	}

	void SoundManager::Release()
	{
		for (int i = 0; i < m_audios.size(); i++)
		{
			delete m_audios[i];
		}
		m_audios.clear();
		for (int i = 0; i < m_buffers.size(); i++)
		{
			delete m_buffers[i];
		}
		m_buffers.clear();
		alcMakeContextCurrent(NULL);
		alcDestroyContext(m_pContext);
		alcCloseDevice(m_pDevice);
		Debug::RELEASESUCCESS("SoundManager");
	}
}