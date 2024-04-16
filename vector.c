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

int v3_ray_plane(v3_t *plane_coord, v3_t *plane_n, v3_t *ray_coord, v3_t *ray_v, v3_t *intersection) {
    // See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection,
    // "Algebraic form".
    double denom = v3_dot(ray_v, plane_n);
    if ((denom > -1e-6) && (denom < 1e-6)) {
        // The ray and the plane are parallel.
        return 0;
    }
    v3_t val;
    v3_sub(plane_coord, ray_coord, &val);
    double d = v3_dot(&val, plane_n) / denom;
    intersection->x = ray_coord->x + (ray_v->x * d);
    intersection->y = ray_coord->y + (ray_v->y * d);
    intersection->z = ray_coord->z + (ray_v->z * d);
    return 1;
}

void v3_rotate(v3_t *in, double yaw, double pitch, double roll, v3_t *out) {
    // Decomposed intrinsic rotations, borrowed from:
    // https://danceswithcode.net/engineeringnotes/rotations_in_3d/rotations_in_3d_part1.html
    double x = in->x,
           y = in->y,
           z = in->z;
    double x0, y0, z0;
    // roll
    // TODO clean this shit up
    x0 = x; y0 = y; z0 = z;
    x = x0;
    y = (y0 * cos(roll)) - (z0 * sin(roll));
    z = (y0 * sin(roll)) + (z0 * cos(roll));
    // pitch
    x0 = x; y0 = y; z0 = z;
    x = (x0 * cos(pitch)) + (z0 * sin(pitch));
    y = y0;
    z = -(x0 * sin(pitch)) + (z0 * cos(pitch));
    // yaw
    x0 = x; y0 = y; z0 = z;
    x = (x0 * cos(yaw)) - (y0 * sin(yaw));
    y = (x0 * sin(yaw)) + (y0 * cos(yaw));
    z = z0;
    out->x = x;
    out->y = y;
    out->z = z;
}

void v3_rodrigues(v3_t *v, v3_t *k, double theta) {
    v3_t c;
    v3_cross(k, v, &c);
    v->x = (v->x * cos(theta)) + (c.x * sin(theta)) + (k->x * v3_dot(k, v) * (1 - cos(theta)));
    v->y = (v->y * cos(theta)) + (c.y * sin(theta)) + (k->y * v3_dot(k, v) * (1 - cos(theta)));
    v->z = (v->z * cos(theta)) + (c.z * sin(theta)) + (k->z * v3_dot(k, v) * (1 - cos(theta)));
}

void v3_rotate_like_a_plane(v3_t *camera_fwd, double roll, v3_t *camera_left, double pitch, v3_t *camera_up, double yaw) {
    v3_rodrigues(camera_left, camera_fwd, roll);
    v3_rodrigues(camera_up, camera_fwd, roll);
    v3_rodrigues(camera_fwd, camera_left, pitch);
    v3_rodrigues(camera_up, camera_left, pitch);
    v3_rodrigues(camera_fwd, camera_up, yaw);
    v3_rodrigues(camera_left, camera_up, yaw);
}
