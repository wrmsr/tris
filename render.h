#ifndef __RENDER_H__
#define __RENDER_H__

#include <stdint.h>

#include "vector.h"

#define FAKESCREEN_W 640
#define FAKESCREEN_H 480

typedef struct {
    double start_time;
    int frames_drawn;
    int pixels_rejected_by_z;
    int pixels_drawn;
} render_stats_t;

typedef struct {
    double start_time;

    int n_skipped_triangles;
    int n_skipped_triangles_by_backface_cull;
    int n_skipped_triangles_by_near_plane;
    int num_triangles;
    int num_onespan_triangles;
    int num_twospan_triangles;
    int num_degenerate_triangles;
    int num_pointup_spans;
    int num_pointdown_spans;
    int num_degenerate_spans;

    int pixels_rejected_by_z;
    int pixels_drawn;
} render_frame_stats_t;

typedef struct {
    uint8_t r, g, b, a;
} rgba_t;

typedef struct {
    rgba_t *texture;
    int w, h;
} material_t;

typedef struct {
    v3_t xyz;
    struct {
        double u, v;
    } uv;
} triangle_coord_t;

typedef struct {
    union {
        triangle_coord_t abc[3];
        struct {
            triangle_coord_t a, b, c;
        };
    };
    material_t *material;
} triangle_t;

// Triangles in screen space, the visible area of which is (0, 0) (screen
// width, screen height). The Z-value is used only to populate the Z-buffer.
typedef struct {
    int16_t x, y;
    double z;
} screen_vertex_t;

typedef struct {
    screen_vertex_t v[3];
    uint8_t r, g, b;
    triangle_t *parent;
} screen_triangle_t;

#define TRIANGLE_POOL_SIZE 20000

// Any triangle can be decomposed into one or two triangles with a flat top or
// bottom (a "span"). A span is simple and fast to draw.
#define NUM_SPANS 20000
typedef struct span_t {
    // screen space: needed for quickly figuring out which pixels on screen to draw
    int16_t y_lo, y_hi;
    screen_vertex_t ref; // this is either the lowest (down pointing span) or highest (up pointing span) point of the span.
    double dx_dy_lo;
    double dx_dy_hi;

    // object space: 
    double dz_dy_lo;
    double dz_dx_lo;

    // needed for materials, original triangle coords for texture mapping, etc
    triangle_t *triangle;

    // We'll insert spans in a linked list (the "y range table"). That table
    // indicates on which y-value spans begin and end, to speed up the raster
    // loop.
    struct span_t *next_span_y_entry;
    struct span_t *next_span_y_exit;

    // At raster time, we'll keep a list of active spans.
    struct span_t *next_active_span;
    struct span_t *prev_active_span;
} span_t;

void render_begin_frame(void);
triangle_t *render_add_triangle(void);
void render_draw_screen_triangle(screen_triangle_t *);
void render_draw_triangle(v3_t *, v3_t *, v3_t *, v3_t *, triangle_t *);
void render_end_frame(void);
void render_print_stats(void);

#endif
