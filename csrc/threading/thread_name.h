// Thread names for `top -H` and profilers. Linux keeps at most 15 characters; errors are
// ignored.
#pragma once

#include <pthread.h>

#include <cstdio>

namespace fa { namespace cpu { namespace threading {

inline void name_current_thread(const char* name) {
#if defined(__APPLE__)
    (void)pthread_setname_np(name);
#elif defined(__linux__)
    (void)pthread_setname_np(pthread_self(), name);
#else
    (void)name;
#endif
}

// "prefix-N", e.g. "fa-map-17".
inline void name_current_thread(const char* prefix, int index) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "%s-%d", prefix, index);
    name_current_thread(buf);
}

}}}  // namespace fa::cpu::threading
