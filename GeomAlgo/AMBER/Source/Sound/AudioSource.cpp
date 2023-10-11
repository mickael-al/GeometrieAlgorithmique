#include "AudioSource.hpp"

namespace Ge
{
	AudioSource::AudioSource(SoundBuffer* sb, std::string name) : GObject()
	{
		alGenSources(1, &m_sourceID);
		alSourcei(m_sourceID, AL_BUFFER, sb->getbufferID());
		alSourcef(m_sourceID, AL_PITCH, m_pitch);
		alSourcef(m_sourceID, AL_GAIN, m_gain);
		setName(name);
		m_velocity = glm::vec3(0, 0, 0);
	}

	void AudioSource::play()
	{
		alSourcePlay(m_sourceID);
	}

	void AudioSource::stop()
	{
		alSourceStop(m_sourceID);
	}

	void AudioSource::pause()
	{
		alSourcePause(m_sourceID);
	}

	void AudioSource::setPitch(float pitch)
	{
		m_pitch = pitch;
		alSourcef(m_sourceID, AL_PITCH, m_pitch);
	}

	void AudioSource::setGain(float gain)
	{
		m_gain = gain;
		alSourcef(m_sourceID, AL_GAIN, m_gain);
	}

	ALint AudioSource::SourceState()
	{
		ALint sourceState;
		alGetSourcei(m_sourceID, AL_SOURCE_STATE, &sourceState);
		//AL_PLAYING
		return sourceState;
	}

	void AudioSource::setPosition(glm::vec3 pos)
	{
		GObject::setPosition(pos);
		alSource3f(m_sourceID, AL_POSITION, pos.x, pos.y, pos.z);
	}

	void AudioSource::setRotation(glm::quat rot)
	{
		GObject::setRotation(rot);
		glm::vec3 dir = getDirection();		
		orientation[0] = dir.x;
		orientation[1] = dir.y;
		orientation[2] = dir.z;
		orientation[3] = 0;
		orientation[4] = -1.0f;
		orientation[5] = 0.0f;		
		alSourcefv(m_sourceID, AL_ORIENTATION, orientation);
	}

	void AudioSource::setEulerAngles(glm::vec3 eul)
	{
		GObject::setEulerAngles(eul);
		glm::vec3 dir = getDirection();
		orientation[0] = dir.x;
		orientation[1] = dir.y;
		orientation[2] = dir.z;
		orientation[3] = 0;
		orientation[4] = -1.0f;
		orientation[5] = 0.0f;
		alSourcefv(m_sourceID, AL_ORIENTATION, orientation);
	}

	void AudioSource::setVelocity(glm::vec3 vel)
	{
		m_velocity = vel;
		alSource3f(m_sourceID, AL_VELOCITY, m_velocity.x, m_velocity.y, m_velocity.z);
	}

	glm::vec3 AudioSource::getVelocity() const
	{
		return m_velocity;
	}

	void AudioSource::setLoop(bool state)
	{
		m_loop = state;
		alSourcei(m_sourceID, AL_LOOPING, m_loop ? AL_TRUE : AL_FALSE);
	}

	bool AudioSource::getLoop() const
	{
		return m_loop;
	}

	void AudioSource::setRolloffFactor(float rolloffFactor)
	{
		m_rolloffFactor = rolloffFactor;
		alSourcef(m_sourceID, AL_ROLLOFF_FACTOR, m_rolloffFactor);
	}

	float AudioSource::getRolloffFactor() const
	{
		return m_rolloffFactor;
	}

	void AudioSource::setMaxDistance(float maxDistance)
	{
		m_maxDistance = maxDistance;
		alSourcef(m_sourceID, AL_MAX_DISTANCE, m_maxDistance);
	}

	float AudioSource::getMaxDistance() const
	{
		return m_maxDistance;
	}

	void AudioSource::setRefDistance(float refDistance)
	{
		m_refDistance = refDistance;
		alSourcef(m_sourceID, AL_REFERENCE_DISTANCE, m_refDistance);
	}

	float AudioSource::getRefDistance() const
	{
		return m_refDistance;
	}

	void AudioSource::mapMemory(){}

	void AudioSource::onGUI()
	{
		GObject::onGUI();
		ImGui::TextColored(ImVec4(0.2f, 1, 0.2f, 1), "AudioSource\n");

		if (ImGui::DragFloat("Pitch", &m_pitch, 0.05f, 0.0f))
		{
			setPitch(m_pitch);
		}
		if (ImGui::DragFloat("Gain", &m_gain, 0.05f, 0.0f,1.0f))
		{
			setGain(m_gain);
		}
		if (ImGui::DragFloat("RolloffFactor", &m_rolloffFactor,0.05f,0.0f))
		{
			setRolloffFactor(m_rolloffFactor);
		}
		if (ImGui::DragFloat("MaxDistance", &m_maxDistance, 0.5f, 0.0f))
		{
			setMaxDistance(m_maxDistance);
		}
		if (ImGui::DragFloat("RefDistance", &m_refDistance, 0.5f, 0.0f))
		{
			setRefDistance(m_refDistance);
		}
		if (ImGui::DragFloat3("Velocity", (float*)&m_velocity))
		{
			setVelocity(m_velocity);
		}
	}

	AudioSource::~AudioSource()
	{
		alSourceStop(m_sourceID);
		alDeleteSources(1, &m_sourceID);
	}
}