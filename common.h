#ifndef __COMMON_H__
#define __COMMON_H__

#include <stdlib.h>

#define ASSERT(cond) \
    { \
        if (!(cond)) { \
            printf("ASSERT(%s) failed at %s L%i\n", #cond, __FILE__, __LINE__); \
            exit(1); \
        } \
    }

#endif
