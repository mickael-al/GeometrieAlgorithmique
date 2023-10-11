#include "SoundBuffer.hpp"

namespace Ge
{
	SoundBuffer::SoundBuffer(uint64_t size, uint32_t frequency, ALenum format, int8_t* buffer)
	{
		alGenBuffers(1, &m_bufferID);
		alBufferData(m_bufferID, format, buffer, size, frequency);
	}

	ALuint SoundBuffer::getbufferID() const
	{
		return m_bufferID;
	}

	SoundBuffer::~SoundBuffer()
	{
		alDeleteBuffers(1, &m_bufferID);
	}
}