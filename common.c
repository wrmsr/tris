#include <time.h>

#include "common.h"

double old_t = 0;
double time_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    double t = ts.tv_sec + (ts.tv_nsec / 1.0e9);
    ASSERT(t > old_t);
    old_t = t;
    return t;
}
