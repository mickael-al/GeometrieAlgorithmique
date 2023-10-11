#ifndef __ENGINE_INITIALIZER__
#define __ENGINE_INITIALIZER__

namespace Ge
{
    class Initializer
    {
        public:
            virtual bool initialize() = 0;
            virtual void release() = 0;
    };
}

#endif // __ENGINE_INITIALIZER__