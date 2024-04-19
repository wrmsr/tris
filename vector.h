#ifndef __VECTOR_H__
#define __VECTOR_H__

typedef union {
    struct {
        double x, y, z;
    };
    double v3[3];  // TODO: xyz
} v3_t;

void v3_add(v3_t *, v3_t *, v3_t *);
void v3_cross(v3_t *, v3_t *, v3_t *);
double v3_dot(v3_t *, v3_t *);
double v3_len(v3_t *);
void v3_sub(v3_t *, v3_t *, v3_t *);
int v3_ray_plane(v3_t *, v3_t *, v3_t *, v3_t *, v3_t *);
void v3_rotate(v3_t *in, double, double, double, v3_t *);
void v3_rodrigues(v3_t *v, v3_t *k, double);
void v3_rotate_like_a_plane(v3_t *, double, v3_t *, double, v3_t *, double);

#endif
