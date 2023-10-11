#include "SoundTypeLoader.hpp"

namespace Ge
{
    int8_t* SoundTypeLoader::LoadWavFormat(const char* filePath, uint64_t* size, uint32_t* frequency, int32_t* format)
    {
        RIFF_Header riffHeader;
        WAVE_Format waveFormat;
        int8_t* buffer = nullptr;

        FILE* wavFile = fopen(filePath, "rb");
        if (wavFile == nullptr)
        {
            Debug::Error("Unable to open wave file: %s\n", filePath);
            return buffer;
        }

        fread(&riffHeader, 1, sizeof(RIFF_Header), wavFile);

        if (strncmp((char*)&riffHeader.chunkID, "RIFF", 4) != 0 || strncmp((char*)&riffHeader.format, "WAVE", 4) != 0)
        {            
            Debug::Error("invalid .wav file");            
            fclose(wavFile);
            return buffer;
        }

        fread(&waveFormat, sizeof(WAVE_Format), 1, wavFile);
        if (strncmp((char*)&waveFormat.subChunkID, "fmt ", 4) != 0)
        {
            Debug::Error("invalid fmt .wav file");
            fclose(wavFile);
            return buffer;
        }
        
        uint8_t * subChunkID = new uint8_t[4];
        fread(subChunkID, sizeof(char)*4, 1, wavFile);
        if (strncmp((char*)subChunkID, "LIST", 4) == 0)
        {
            uint32_t size;
            fread(&size, sizeof(uint32_t), 1, wavFile);
            fseek(wavFile, sizeof(uint8_t)* size, SEEK_CUR);
            fread(subChunkID, sizeof(char) * 4, 1, wavFile);            
        }

        if (strncmp((char*)subChunkID, "data", 4) == 0)
        {
            uint32_t data_size;
            fread(&data_size, sizeof(uint32_t), 1, wavFile);
            buffer = new int8_t[data_size];
            fread(buffer, 1, data_size, wavFile);

            *size = data_size;
            *frequency = waveFormat.sampleRate;

            if (waveFormat.bitsPerSample == 8)
            {
                *format = (waveFormat.numChannels == 1) ? AL_FORMAT_MONO8 : AL_FORMAT_STEREO8;
            }
            else if (waveFormat.bitsPerSample == 16)
            {
                *format = (waveFormat.numChannels == 1) ? AL_FORMAT_MONO16 : AL_FORMAT_STEREO16;
            }
        }
        delete subChunkID;
        fclose(wavFile);
        return buffer;
    }

    int SoundTypeLoader::getFileSize(FILE* inFile)
    {
        int fileSize = 0;
        fseek(inFile, 0, SEEK_END);

        fileSize = ftell(inFile);

        fseek(inFile, 0, SEEK_SET);
        return fileSize;
    }
}