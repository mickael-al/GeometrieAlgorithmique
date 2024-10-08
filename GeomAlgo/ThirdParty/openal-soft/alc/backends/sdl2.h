#ifndef BACKENDS_SDL2_H
#define BACKENDS_SDL2_H

#include "backends/base.h"

struct SDL2BackendFactory final : public BackendFactory {
public:
    bool init() override;

    bool querySupport(BackendType type) override;

    std::string probe(BackendType type) override;

    BackendPtr createBackend(ALCdevice *device, BackendType type) override;

    static BackendFactory &getFactory();
};

#endif /* BACKENDS_SDL2_H */
