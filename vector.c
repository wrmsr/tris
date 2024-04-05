#include <math.h>

#include "vector.h"

void v3_add(v3_t *a, v3_t *b, v3_t *c) {
    c->x = a->x + b->x;
    c->y = a->y + b->y;
    c->z = a->z + b->z;
}

void v3_cross(v3_t *a, v3_t *b, v3_t *c) {
    c->x = (a->y * b->z) - (a->z * b->y);
    c->y = (a->z * b->x) - (a->x * b->z);
    c->z = (a->x * b->y) - (a->y * b->x);
}

double v3_dot(v3_t *a, v3_t *b) {
    return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

double v3_len(v3_t *v) {
    return sqrt((v->x * v->x) + (v->y * v->y) + (v->z * v->z));
}

void v3_sub(v3_t *a, v3_t *b, v3_t *c) {
    c->x = a->x - b->x;
    c->y = a->y - b->y;
    c->z = a->z - b->z;
}
