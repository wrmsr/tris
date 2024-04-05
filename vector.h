#ifndef __VECTOR_H__
#define __VECTOR_H__

typedef union {
    struct {
        double x, y, z;
    };
    double v3[3];
} v3_t;

void v3_add(v3_t *, v3_t *, v3_t *);
void v3_cross(v3_t *, v3_t *, v3_t *);
double v3_dot(v3_t *, v3_t *);
double v3_len(v3_t *);
void v3_sub(v3_t *, v3_t *, v3_t *);

#endif
