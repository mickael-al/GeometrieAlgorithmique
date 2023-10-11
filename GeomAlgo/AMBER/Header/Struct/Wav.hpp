#ifndef __WAV__
#define __WAV__

#include <iostream>

struct RIFF_Header 
{
    uint8_t chunkID[4];
    uint32_t chunkSize;
    uint8_t format[4];
};

struct WAVE_Format 
{
    uint8_t subChunkID[4];
    uint32_t subChunkSize;
    uint16_t audioFormat;
    uint16_t numChannels;
    uint32_t sampleRate;
    uint32_t byteRate;
    uint16_t blockAlign;
    uint16_t bitsPerSample;
};

#endif //!__WAV__