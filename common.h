#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdlib.h>

#define MAX(a,b) ({ \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a > _b ? _a : _b; \
})

#define MIN(a,b) ({ \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b; \
})

#define ASSERT(cond) { \
    if (!(cond)) { \
        printf("ASSERT(%s) failed at %s L%i\n", #cond, __FILE__, __LINE__); \
        exit(1); \
    } \
}

#endif
