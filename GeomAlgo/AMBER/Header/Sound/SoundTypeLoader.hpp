#ifndef __SOUND_TYPE_LOADER__
#define __SOUND_TYPE_LOADER__

#include <iostream>
#include <string>
#include <fstream>
#include "Debug.hpp"

using namespace std;
using std::string;
using std::fstream;

#include "Wav.hpp"
#include "AL/al.h"

namespace Ge
{
	class SoundTypeLoader
	{
	public:
		static int8_t* LoadWavFormat(const char* filePath, uint64_t* size, uint32_t* frequency, int32_t* format);
	private:
		static int getFileSize(FILE* inFile);
	};
}

#endif //!__SOUND_TYPE_LOADER__