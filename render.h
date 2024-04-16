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
    v3_t v[3];
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
} screen_triangle_t;

// Any triangle can be decomposed into one or two triangles with a flat top or
// bottom (a "span"). A span is simple and fast to draw.
#define NUM_SPANS 20000
typedef struct span_t {
    int16_t y_lo, y_hi;
    screen_vertex_t ref; // this is either the lowest (down pointing span) or highest (up pointing span) point of the span.
    double dx_dy_lo;
    double dx_dy_hi;
    double dz_dy_lo;
    double dz_dx_lo;
    //screen_triangle_t *parent; // the triangle this span comes from, so that we can do texture/attribute/etc lookups

    // We'll insert spans in a linked list (the "y range table"). That table
    // indicates on which y-value spans begin and end, to speed up the raster
    // loop.
    struct span_t *next_span_y_entry;
    struct span_t *next_span_y_exit;

    // At raster time, we'll keep a list of active spans.
    struct span_t *next_active_span;
    struct span_t *prev_active_span;
} span_t;

void begin_frame(void);
void draw_screen_triangle(screen_triangle_t *);
void draw(v3_t *, v3_t *, v3_t *, v3_t *, triangle_t *);
void end_frame(void);
void print_stats(void);

#endif
