#include <GL/gl.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdint.h>
#include <time.h>

#include "common.h"
#include "vector.h"

#define EPSILON 1e-6

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
span_t spans[NUM_SPANS];
int num_spans;

#define FAKESCREEN_W 640
#define FAKESCREEN_H 480

span_t *span_y_entry_table[FAKESCREEN_H];
span_t *span_y_exit_table[FAKESCREEN_H];

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

void add_span(screen_vertex_t *a, screen_vertex_t *b, screen_vertex_t *c, screen_vertex_t *hi, screen_vertex_t *lo) {
    ASSERT((hi == NULL) != (lo == NULL));
    screen_vertex_t *hi_or_lo = (hi != NULL ? hi : lo);
    // Of the other two vertices, which is on the left, and which is on the
    // right?
    span_t *span = &spans[num_spans];
    screen_vertex_t *other1 = (a == hi_or_lo ? b : a);
    screen_vertex_t *other2 = (other1 == b ? c : (b == hi_or_lo ? c : b));
    screen_vertex_t *x_lo_vert, *x_hi_vert;
    if (other1->x < other2->x) {
        x_lo_vert = other1;
        x_hi_vert = other2;
    } else if (other2->x < other1->x) {
        x_lo_vert = other2;
        x_hi_vert = other1;
    } else {
        // We may receive degenerate spans, just make a note.
        num_degenerate_spans++;
        return;
    }
    if (hi != NULL) {
        num_pointdown_spans++;
        span->y_lo = other1->y; // or other2[1], doesn't matter.
        span->y_hi = hi->y;
    } else {
        num_pointup_spans++;
        span->y_lo = lo->y;
        span->y_hi = other1->y; // or other2[1], doesn't matter.
    }
    memcpy(&(span->ref), hi_or_lo, sizeof(span->ref));
    span->dx_dy_lo = (span->ref.x - x_lo_vert->x) / (double)(span->ref.y - x_lo_vert->y);
    span->dx_dy_hi = (span->ref.x - x_hi_vert->x) / (double)(span->ref.y - x_hi_vert->y);
    span->dz_dy_lo = (span->ref.z - x_lo_vert->z) / (double)(span->ref.y - x_lo_vert->y);
    span->dz_dx_lo = (x_hi_vert->z - x_lo_vert->z) / (double)(x_hi_vert->x - x_lo_vert->x);
    ASSERT(++num_spans < NUM_SPANS);

    int clamped_y_lo = MAX(0, MIN(FAKESCREEN_H - 1, span->y_lo));
    int clamped_y_hi = MAX(0, MIN(FAKESCREEN_H - 1, span->y_hi));
    span->next_span_y_entry = span_y_entry_table[clamped_y_lo];
    span_y_entry_table[clamped_y_lo] = span;
    span->next_span_y_exit = span_y_exit_table[clamped_y_hi];
    span_y_exit_table[clamped_y_hi] = span;
}

int ray_plane(v3_t *plane_coord, v3_t *plane_n, v3_t *ray_coord, v3_t *ray_v, v3_t *intersection) {
    // See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection,
    // "Algebraic form".
    double denom = v3_dot(ray_v, plane_n);
    if ((denom > -EPSILON) && (denom < EPSILON)) {
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

void rotate(v3_t *in, double yaw, double pitch, double roll, v3_t *out) {
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

void rodrigues(v3_t *v, v3_t *k, double theta) {
    v3_t c;
    v3_cross(k, v, &c);
    v->x = (v->x * cos(theta)) + (c.x * sin(theta)) + (k->x * v3_dot(k, v) * (1 - cos(theta)));
    v->y = (v->y * cos(theta)) + (c.y * sin(theta)) + (k->y * v3_dot(k, v) * (1 - cos(theta)));
    v->z = (v->z * cos(theta)) + (c.z * sin(theta)) + (k->z * v3_dot(k, v) * (1 - cos(theta)));
}

void rotate_like_a_plane(v3_t *camera_fwd, double roll, v3_t *camera_left, double pitch, v3_t *camera_up, double yaw) {
    rodrigues(camera_left, camera_fwd, roll);
    rodrigues(camera_up, camera_fwd, roll);
    rodrigues(camera_fwd, camera_left, pitch);
    rodrigues(camera_up, camera_left, pitch);
    rodrigues(camera_fwd, camera_up, yaw);
    rodrigues(camera_left, camera_up, yaw);
}

void draw_screen_triangle(screen_triangle_t *);

void draw(v3_t *camera_pos, v3_t *camera_fwd, v3_t *camera_up, v3_t *camera_left, triangle_t *triangle) {
    // The screen exists on a plane "in front of" the camera.
    v3_t screen_center;
    v3_add(camera_pos, camera_fwd, &screen_center);

    // Do backface culling: a triangle facing the wrong way from the
            // camera (according to its normal) doesn't get rendered.
            v3_t v1, v2;
            v3_sub(&(triangle->v[0]), &(triangle->v[1]), &v1);
            v3_sub(&(triangle->v[0]), &(triangle->v[2]), &v2);
            v3_t normal;
            v3_cross(&v2, &v1, &normal);
            if (v3_dot(camera_fwd, &normal) < 0) {
                n_skipped_triangles++;
                n_skipped_triangles_by_backface_cull++;
                return;
            }

            // Turn each arbitrary triangle into triangles suitable for
            // rendering, by clipping away the parts that cannot be correctly
            // rendered.
            int n_inside = 0, inside[3],
                n_outside = 0, outside[3];
            for (int j = 0; j < 3; j++) {
                // Calculate signed distance to figure out which side of the
                // screen plane the vertex is on.
                v3_t screen_to_v;
                v3_sub(&(triangle->v[j]), &screen_center, &screen_to_v);
                double d = v3_dot(&screen_to_v, camera_fwd);
                if (d < 0) {
                    inside[n_inside++] = j;
                    ASSERT(n_inside <= 3);
                } else {
                    outside[n_outside++] = j;
                    ASSERT(n_outside <= 3);
                }
            }

            int n_clipped_triangles = 0;
            triangle_t clipped_triangles[2];
            if (n_inside == 0) {
                // This triangle is fine as-is.
                n_clipped_triangles = 1;
                memcpy(&clipped_triangles[0], triangle, sizeof(triangle_t));
            } else if (n_inside == 1) {
                ASSERT(n_outside == 2);
                // If one point is inside (v3), then when we clip the triangle
                // it becomes a quadrilateral. The quadrilateral has two of the
                // original vertices, and an additional 2 vertices (named u1,
                // u2) at the point where (v3, v1) and (v3, v2) intersect the
                // screen plane.
                v3_t *v3 = &(triangle->v[inside[0]]);
                v3_t *v1 = &(triangle->v[outside[0]]);
                v3_t v1_to_v3;
                v3_sub(v3, v1, &v1_to_v3);
                v3_t u1;
                // v1 and v3 are on opposite sides of the camera plane, so
                // v3-v1 and the camera plane cannot be parallel.
                ASSERT(ray_plane(&screen_center, camera_fwd, v1, &v1_to_v3, &u1));
                v3_t *v2 = &(triangle->v[outside[1]]);
                v3_t v2_to_v3;
                v3_sub(v3, v2, &v2_to_v3);
                v3_t u2;
                // Ditto, see above.
                ASSERT(ray_plane(&screen_center, camera_fwd, v2, &v2_to_v3, &u2));
                // We can't render a quadrilateral, but we can split it into
                // two triangles. Since we know that the perimeter of the
                // quadrilateral is formed by visiting v1, v2, u2, u1 in that
                // order, we know that we can form two triangles (v1, v2, u2)
                // and (u2, u1, v1).
                n_clipped_triangles = 2;
                memcpy(&clipped_triangles[0].v[0], v1, sizeof(v3_t));
                memcpy(&clipped_triangles[0].v[1], v2, sizeof(v3_t));
                memcpy(&clipped_triangles[0].v[2], &u2, sizeof(v3_t));
                memcpy(&clipped_triangles[1].v[0], &u2, sizeof(v3_t));
                memcpy(&clipped_triangles[1].v[1], &u1, sizeof(v3_t));
                memcpy(&clipped_triangles[1].v[2], v1, sizeof(v3_t));
            } else if (n_inside == 2) {
                ASSERT(n_outside == 1);
                // if two points inside (v1, v2), then make one triangle (v3, u1, u2)
                v3_t *v3 = &(triangle->v[outside[0]]);
                v3_t *v1 = &(triangle->v[inside[0]]);
                v3_t v1_to_v3;
                v3_sub(v3, v1, &v1_to_v3);
                v3_t u1;
                // v1 and v3 are on opposite sides of the camera plane, so
                // v3-v1 and the camera plane cannot be parallel.
                ASSERT(ray_plane(&screen_center, camera_fwd, v1, &v1_to_v3, &u1));
                v3_t *v2 = &(triangle->v[inside[1]]);
                v3_t v2_to_v3;
                v3_sub(v3, v2, &v2_to_v3);
                v3_t u2;
                // Ditto, see above.
                ASSERT(ray_plane(&screen_center, camera_fwd, v2, &v2_to_v3, &u2));
                n_clipped_triangles = 1;
                memcpy(&clipped_triangles[0].v[0], v3, sizeof(v3_t));
                memcpy(&clipped_triangles[0].v[1], &u1, sizeof(v3_t));
                memcpy(&clipped_triangles[0].v[2], &u2, sizeof(v3_t));
            } else if (n_inside == 3) {
                // Just don't render this triangle.
                ASSERT(n_outside == 0);
                n_skipped_triangles++;
                n_skipped_triangles_by_near_plane++;
            } else {
                ASSERT(0);
            }
            ASSERT((n_clipped_triangles >= 0) && (n_clipped_triangles <= 2));

                // Anything on the same side of the screen plane is the camera is "inside the camera".
                //We can't render (partially or
                // completely) inside objects correctly, and we have to do something about it. Example:
                
                    /*
                              v1 . . . . . . . v2
                               .               .
                                .             .
                    _*___________.___________.____ camera plane
                                  .         .
                                   .       .
                           \cam/    .     .
                                     .   .
                                      . .
                                       v3

                    in this case, cam->v3 does intersect the plane, but the
                    intersection point is to the left (see *), which results in
                    an unintuitive and weird render, given that the other
                    vertices of the triangle will be projected onto the right of the screen.

                    however, not rendering any part of the triangle would also
                    look bad. to make this triangle renderable, we can remove the part of
                    the triangle that won't render properly, which in this case
                    results in a quadrilateral. We can split it into two
                    nonoverlapping triangles.

                              v1 . . . . . . . v2
                               .'  ,           .
                                .     ' .     .
                    _____________. . . . . .'.____ camera plane


                           \cam/

                    there are three cases in total for how to do this culling,
                    depending on how many vertices (in this case, 1) are
                    "inside".
                    */

                // side_of_plane(screen_center, camera_fwd, point) is negative
                // if point is between the camera and the screen plane
                // ("inside" the camera). 

            // Project each clipped triangle into screen space.
            for (int j = 0; j < n_clipped_triangles; j++) {
            screen_triangle_t screen_triangle;
            for (int k = 0; k < 3; k++) {
                v3_t *v = &clipped_triangles[j].v[k];

                // Project this vertex into screen space: draw a ray from the
                // vertex to the camera, intersecting with a plane (the
                // screen).
                v3_t camera_pos_to_v;
                v3_sub(v, camera_pos, &camera_pos_to_v);
                v3_t poi;
                if (ray_plane(&screen_center, camera_fwd, camera_pos, &camera_pos_to_v, &poi)) {
                    // Find the intersection, relative to the plane center.
                    v3_t poi_rel;
                    v3_sub(&poi, &screen_center, &poi_rel);

                    // Put this vertex into screenspace.
                    double half_screen_w = FAKESCREEN_W / 2,
                           half_screen_h = FAKESCREEN_H / 2;
                    screen_triangle.v[k].x = half_screen_w - (v3_dot(&poi_rel, camera_left) * half_screen_w);
                    screen_triangle.v[k].y = half_screen_h + (v3_dot(&poi_rel, camera_up) * half_screen_h);
                    // If the vertex is behind the camera, then the distance
                    // should be negative.  i dont think this is right TODO maybe remove this assert
                    double z = v3_dot(&camera_pos_to_v, camera_fwd);
                    ASSERT(z >= 0);
                    screen_triangle.v[k].z = z;

                    // For debugging purposes, give every triangle a different color.
                    screen_triangle.r = 255; //(i * 13) % 0xff;
                    screen_triangle.g = 0; //(i * 101) % 0xff;
                    screen_triangle.b = 0; //(i * 211) % 0xff;
                } else {
                    // If the ray doesn't project onto the screen, it's because
                    // the ray is parallel to the screen, so it will be
                    // perpendicular to the screen normal.
                    ASSERT(0);
                }
            }
            draw_screen_triangle(&screen_triangle);
            }
}

void draw_screen_triangle(screen_triangle_t *screen_triangle) {
    // Pick a vertex that is above all other vertices.
    screen_vertex_t *a = &(screen_triangle->v[0]);
    screen_vertex_t *b = &(screen_triangle->v[1]);
    screen_vertex_t *c = &(screen_triangle->v[2]);
    screen_vertex_t *hi = NULL;
    if ((a->y > b->y) && (a->y > c->y)) {
        hi = a;
    } else if ((b->y > a->y) && (b->y > c->y)) {
        hi = b;
    } else if ((c->y > a->y) && (c->y > b->y)) {
        hi = c;
    }

    // Pick a vertex that is below all other vertices.
    screen_vertex_t *lo = NULL;
    if ((a->y < b->y) && (a->y < c->y)) {
        lo = a;
    } else if ((b->y < a->y) && (b->y < c->y)) {
        lo = b;
    } else if ((c->y < a->y) && (c->y < b->y)) {
        lo = c;
    }

    // If there's neither a hi nor lo vertex, then the triangle has no area and
    // cannot be rendered.
    if ((hi == NULL) && (lo == NULL)) {
        num_degenerate_triangles++;
        return;
    }

    // If there's only a hi or lo vertex, then there is only one span, e.g.
    // this can be rendered as one span:
    //   *******
    //    *   *
    //     * *
    //      *
    if ((hi == NULL) != (lo == NULL)) {
        if (hi != NULL) {
            add_span(a, b, c, hi, NULL);
        } else if (lo != NULL) {
            add_span(a, b, c, NULL, lo);
        } else ASSERT(0);
        num_onespan_triangles++;
    }

    // If there is both a hi and lo vertex, then we need to draw two spans.
    if ((hi != NULL) && (lo != NULL)) {
        // First, we need to find the vertex ('mid') which isn't the hi or lo
        // vertex.
        screen_vertex_t *mid = NULL;
        if ((a != hi) && (a != lo)) mid = a;
        else if ((b != hi) && (b != lo)) mid = b;
        else mid = c;
        ASSERT(mid != NULL);
        // Find the point on the edge linking hi and lo which is at the same
        // y-coordinate as 'mid'.
        screen_vertex_t split = {
            .x = lo->x + (((hi->x - lo->x) / (double)(hi->y - lo->y)) * (mid->y - lo->y)),
            .y = mid->y,
            .z = lo->z + (((hi->z - lo->z) / (double)(hi->y - lo->y)) * (mid->y - lo->y))
        };
        // Create two spans!
        add_span(hi, mid, &split, hi, NULL);
        add_span(lo, mid, &split, NULL, lo);
        num_twospan_triangles++;
    }
}

double old_t = 0;
double time_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    double t = ts.tv_sec + (ts.tv_nsec / 1.0e9);
    ASSERT(t > old_t);
    old_t = t;
    return t;
}

#define MAX_VERTS 10000
v3_t verts[MAX_VERTS];
int n_verts = 0;
#define MAX_TRIS 10000
int tris[MAX_TRIS][3];
int n_tris = 0;

void load_obj(char *fn, double scale, double yaw, double pitch, double roll) {
    FILE *fp = fopen(fn, "r");
    ASSERT(fp != NULL);
    char line[128];
    int n_verts_at_begin = n_verts;
    while (fgets(line, sizeof(line), fp) != NULL) {
        ASSERT(strlen(line) > 0);
        if (line[strlen(line) - 1] == '\n') {
            line[strlen(line) - 1] = '\0';
        }
        char *token = strtok(line, " ");
        char type = '\0';
        int off = 0;
        while (token != NULL) {
            switch (type) {
            case '\0':
                ASSERT(strlen(token) == 1);
                type = token[0];
                break;
            case 'v':
                ASSERT(off < 4);
                ASSERT(sscanf(token, "%lf", &verts[n_verts].v3[off]) == 1);
                verts[n_verts].v3[off] *= scale;
                if (off == 2) {
                    rotate(&verts[n_verts], yaw, pitch, roll, &verts[n_verts]);
                }
                off++;
                break;
            case 'f':
                ASSERT(off < 3);
                int i;
                ASSERT(sscanf(token, "%i", &i) == 1);
                tris[n_tris][off++] = i - 1 + n_verts_at_begin;
                break;
            // Intentionally no other cases: ignore other
            // directives/comments/etc.
            }
            token = strtok(NULL, " ");
        }
        switch (type) {
        case 'v':
            ASSERT(++n_verts < MAX_VERTS);
            break;
        case 'f':
            ASSERT(++n_tris < MAX_TRIS);
            break;
        }
    }
    ASSERT(feof(fp) != 0);
    ASSERT(ferror(fp) == 0);
}

typedef struct {
    GLubyte r, g, b;
} gl_rgb_t;

int main(void) {
    load_obj("cow-nonormals.obj", 1, 0, 0, -M_PI / 2);

    ASSERT(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) == 0);
    SDL_Window *sdl_window = SDL_CreateWindow("aspng", 0, 0, 100, 100, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    ASSERT(sdl_window != NULL);
    ASSERT(SDL_GL_CreateContext(sdl_window) != NULL);
    GLuint screen_texture;
    glGenTextures(1, &screen_texture);

    double *depth_buffer = malloc(sizeof(double) * FAKESCREEN_W * FAKESCREEN_H);
    gl_rgb_t *texture = malloc(sizeof(gl_rgb_t) * FAKESCREEN_W * FAKESCREEN_H);

    v3_t camera_pos = { .x = -10, .y = 0, .z = 0 };
    v3_t camera_fwd_reset = { .x = 1, .y = 0, .z = 0 };
    v3_t camera_fwd;
    memcpy(&camera_fwd, &camera_fwd_reset, sizeof(v3_t));
    v3_t camera_left_reset = { .x = 0, .y = -1, .z = 0 };
    v3_t camera_left;
    memcpy(&camera_left, &camera_left_reset, sizeof(v3_t)); // TODO necessary? we know it from the cross product
    v3_t camera_up_reset = { .x = 0, .y = 0, .z = 1 };
    v3_t camera_up;
    memcpy(&camera_up, &camera_up_reset, sizeof(v3_t));
    struct {
        double start_time;
        int frames_drawn;
        int pixels_rejected_by_z;
        int pixels_drawn;
    } render_stats = {
        .start_time = time_now(),
        .frames_drawn = 0,
        .pixels_drawn = 0,
        .pixels_rejected_by_z = 0,
    };
    struct {
        double start_time;
        int pixels_rejected_by_z;
        int pixels_drawn;
    } frame_stats;
    int do_quit = 0;
    while (!do_quit && render_stats.frames_drawn < 200) {
        frame_stats.start_time = time_now();
        frame_stats.pixels_rejected_by_z = 0;
        frame_stats.pixels_drawn = 0;

        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            switch (e.type) {
            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    do_quit = 1;
                    break;
                }
                break;
            }
        }

        // Rotate the camera around the origin.
        double now = render_stats.frames_drawn / 10.;
        double theta = 0.25 * now;
        double scale = 3 + (.5 * cos(2 * now));
        camera_pos.x = (2. * cos(theta)) * scale;
        camera_pos.y = (2. * sin(theta)) * scale;
        camera_pos.z = 0;

        // Point the camera at the origin.
        camera_fwd.x = 0 - camera_pos.x;
        camera_fwd.y = 0 - camera_pos.y;
        camera_fwd.z = 0 - camera_pos.z;
        double l = v3_len(&camera_fwd);
        camera_fwd.x = camera_fwd.x / l;
        camera_fwd.y = camera_fwd.y / l;
        camera_fwd.z = camera_fwd.z / l;
        camera_up.x = 0;
        camera_up.y = 0;
        camera_up.z = 1;
        v3_cross(&camera_fwd, &camera_up, &camera_left);

        // TODO clear the screen: shouldn't have to do this out of band
        for (int y = 0; y < FAKESCREEN_H; y++) {
            for (int x = 0; x < FAKESCREEN_W; x++) {
                int off = (y * FAKESCREEN_W) + x;
                depth_buffer[off] = DBL_MAX;
                texture[off].r = texture[off].g = texture[off].b = 0;
            }
        }

        for (int y = 0; y < FAKESCREEN_H; y++) {
            span_y_entry_table[y] = NULL;
            span_y_exit_table[y] = NULL;
        }

        n_skipped_triangles = 0;
        n_skipped_triangles_by_backface_cull = 0;
        n_skipped_triangles_by_near_plane = 0;
        num_pointup_spans = 0;
        num_pointdown_spans = 0;
        num_degenerate_spans = 0;
        num_onespan_triangles = 0;
        num_twospan_triangles = 0;
        num_degenerate_triangles = 0;

        num_spans = 0;
        for (int i = 0; i < n_tris; i++) {
            // Build the triangle to be rendered.
            triangle_t triangle;
            for (int j = 0; j < 3; j++) {
                memcpy(&triangle.v[j], &verts[tris[i][j]], sizeof(v3_t));
            }
            draw(&camera_pos, &camera_fwd, &camera_up, &camera_left, &triangle);
        }

        // Draw all spans to the screen, respecting the z-buffer.
        for (int i = 0; i < num_spans; i++) {
            spans[i].next_active_span = NULL;
            spans[i].prev_active_span = NULL;
        }
        span_t *active_span_table = NULL;
        double min_z = DBL_MAX;
        double max_z = -DBL_MAX;
        for (int y = 0; y < FAKESCREEN_H; y++) {
            // Add spans which are starting.
            for (span_t *span = span_y_entry_table[y]; span != NULL; span = span->next_span_y_entry) {
                span->next_active_span = active_span_table;
                span->prev_active_span = NULL;
                if (active_span_table != NULL)
                    active_span_table->prev_active_span = span;
                active_span_table = span;
            }

            // Remove spans which are ending.
            for (span_t *span = span_y_exit_table[y]; span != NULL; span = span->next_span_y_exit) {
                if (span->prev_active_span != NULL)
                    span->prev_active_span->next_active_span = span->next_active_span;
                else
                    active_span_table = span->next_active_span;
                if (span->next_active_span != NULL)
                    span->next_active_span->prev_active_span = span->prev_active_span;
            }

            // Render every active span.
            for (span_t *span = active_span_table; span != NULL; span = span->next_active_span) {
                int16_t x_fill_lo = span->ref.x + (span->dx_dy_lo * (y - span->ref.y));
                int16_t x_fill_hi = span->ref.x + (span->dx_dy_hi * (y - span->ref.y));
                ASSERT(x_fill_lo <= x_fill_hi);
                if (x_fill_lo < 0)
                    x_fill_lo = 0;
                if (x_fill_hi > FAKESCREEN_W - 1)
                    x_fill_hi = FAKESCREEN_W - 1;
                double z_lo = span->ref.z + (span->dz_dy_lo * (y - span->ref.y));
                for (int16_t x = x_fill_lo; x <= x_fill_hi; x++) {
                    // Do a z-check before we draw the pixel.
                    int off = (y * FAKESCREEN_W) + x;
                    double z = z_lo + (span->dz_dx_lo * (x - x_fill_lo));
                    if ((z < depth_buffer[off]) && (z >= 0)) {
                        depth_buffer[off] = z;
                        if (z > max_z) {
                            max_z = z;
                        }
                        if (z < min_z) {
                            min_z = z;
                        }
                        texture[off].r = 0; //span->parent->r;
                        texture[off].g = 0; //span->parent->g;
                        texture[off].b = 255; //span->parent->b;
                        render_stats.pixels_drawn++;
                        frame_stats.pixels_drawn++;
                    } else {
                        render_stats.pixels_rejected_by_z++;
                        frame_stats.pixels_rejected_by_z++;
                    }
                }
            }
        }

        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, screen_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, FAKESCREEN_W, FAKESCREEN_H, 0, GL_RGB, GL_UNSIGNED_BYTE, texture);

        // Figure out how big of a "screen" we can draw while not exceeding the
        // window's width or height and maintaining the aspect ratio.
        double ratio = FAKESCREEN_W / (double)FAKESCREEN_H;
        SDL_Surface *window_surface = SDL_GetWindowSurface(sdl_window);
        double
          sizing_1_w = window_surface->h * ratio,
          sizing_1_h = sizing_1_w / ratio;
        bool sizing_1_fits = (sizing_1_w <= window_surface->w) && (sizing_1_h <= window_surface->h);
        double
          sizing_2_h = window_surface->w / ratio,
          sizing_2_w = sizing_2_h * ratio;
        bool sizing_2_fits = (sizing_2_w <= window_surface->w) && (sizing_2_h <= window_surface->h);
        double screen_w, screen_h;
        if (sizing_1_fits) {
            screen_w = sizing_1_w;
            screen_h = sizing_1_h;
        } else if (sizing_2_fits) {
            screen_w = sizing_2_w;
            screen_h = sizing_2_h;
        } else ASSERT(0);

        // Draw the "screen" using OpenGL to scale it up to the window.
        glClearColor(1, 0, 1, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, window_surface->w, window_surface->h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(0, window_surface->w, 0, window_surface->h, -1, 1);
        glMatrixMode(GL_MODELVIEW);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, screen_texture);
        double
           left   = (window_surface->w / 2.) - (screen_w / 2),
           right  = (window_surface->w / 2.) + (screen_w / 2),
           top    = (window_surface->h / 2.) + (screen_h / 2),
           bottom = (window_surface->h / 2.) - (screen_h / 2);
        glBegin(GL_QUADS);
        glTexCoord2f(0, 1); glVertex3f(left, bottom, 0);
        glTexCoord2f(1, 1); glVertex3f(right, bottom, 0);
        glTexCoord2f(1, 0); glVertex3f(right, top, 0);
        glTexCoord2f(0, 0); glVertex3f(left, top, 0);
        glEnd();

        SDL_GL_SwapWindow(sdl_window);

        printf(
            "frame %i took %f seconds:\n"
            "    (%i triangles skipped: %i by backface culling, %i by near plane)\n"
            "    %i triangles onscreen (%i onespans, %i twospans, %i degenerate)\n"
            "    %i spans (%i pointups, %i pointdowns, %i degenerates)\n"
            "    (%e pixels rejected by z)\n"
            "    %e pixels drawn\n",
            render_stats.frames_drawn, time_now() - frame_stats.start_time,
            n_skipped_triangles, n_skipped_triangles_by_backface_cull, n_skipped_triangles_by_near_plane,
            num_triangles, num_onespan_triangles, num_twospan_triangles, num_degenerate_triangles,
            num_spans, num_pointup_spans, num_pointdown_spans, num_degenerate_spans,
            (double)frame_stats.pixels_rejected_by_z,
            (double)frame_stats.pixels_drawn
        );
        render_stats.frames_drawn++;
    }

    double elapsed = time_now() - render_stats.start_time;
    printf("%i frames, %e pixels in %f seconds: %f fps, %e pixels/s\n",
        render_stats.frames_drawn,
        (double)render_stats.pixels_drawn,
        elapsed,
        render_stats.frames_drawn / elapsed,
        render_stats.pixels_drawn / elapsed
    );
}
