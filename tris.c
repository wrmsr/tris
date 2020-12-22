#include <math.h>
#include <SDL2/SDL.h>
#include <stdint.h>
#include <time.h>

#include "common.h"

#define EPSILON 1e-6
#define DOT(u, v) ((u[0] * v[0]) + (u[1] * v[1]) + (u[2] * v[2]))
#define MAGNITUDE(v) sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2))

// Triangles in screen space, the visible area of which is (0, 0) (screen
// width, screen height). The Z-value is used only to populate the Z-buffer.
typedef struct {
    int16_t x, y;
    double z;
} screen_vertex;

#define NUM_TRIANGLES 10000
typedef struct {
    screen_vertex v[3];
    uint8_t r, g, b;
} triangle_t;
triangle_t triangles[NUM_TRIANGLES];
int num_triangles = 0;

// Any triangle can be decomposed into one or two triangles with a flat top or
// bottom (a "span").
#define NUM_SPANS 20000
typedef struct {
    int16_t y_lo, y_hi;
    screen_vertex ref;
    double dx_dy_lo;
    double dx_dy_hi;
    double dz_dy_lo;
    double dz_dx_lo;
    triangle_t *parent;
} span_t;
span_t spans[NUM_SPANS];
int num_spans;

int num_pointup_spans;
int num_pointdown_spans;
int num_degenerate_spans;
int num_onespan_triangles;
int num_twospan_triangles;
int num_degenerate_triangles;

void add_span(screen_vertex *a, screen_vertex *b, screen_vertex *c, screen_vertex *hi, screen_vertex *lo, triangle_t *parent) {
    ASSERT((hi == NULL) != (lo == NULL));
    screen_vertex *hi_or_lo = (hi != NULL ? hi : lo);
    // Of the other two vertices, which is on the left, and which is on the
    // right?
    span_t *span = &spans[num_spans];
    screen_vertex *other1 = (a == hi_or_lo ? b : a);
    screen_vertex *other2 = (other1 == b ? c : (b == hi_or_lo ? c : b));
    screen_vertex *x_lo_vert, *x_hi_vert;
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
    span->parent = parent;
    ASSERT(++num_spans < NUM_SPANS);
}

int ray_plane(double *plane_coord, double *plane_n, double *ray_coord,
    double *ray_v, double *intersection) {
    // See https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection,
    // "Algebraic form".
    double denom = DOT(ray_v, plane_n);
    if (denom < EPSILON) {
        // The ray and the plane are parallel.
        return 0;
    }
    double val[3] = {
        plane_coord[0] - ray_coord[0],
        plane_coord[1] - ray_coord[1],
        plane_coord[2] - ray_coord[2]
    };
    double d = DOT(val, plane_n) / denom;
    intersection[0] = ray_coord[0] + (ray_v[0] * d);
    intersection[1] = ray_coord[1] + (ray_v[1] * d);
    intersection[2] = ray_coord[2] + (ray_v[2] * d);
    return 1;
}

void rotate(double *in, double yaw, double pitch, double roll, double *out) {
    // 3x3 rotation matrix for yaw, pitch, and roll, yoinked from:
    // https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions
    out[0] =
        (in[0] * cos(yaw) * cos(roll))
      + (in[1] * ((cos(roll) * sin(pitch) * sin(yaw)) - (sin(roll) * cos(pitch))))
      + (in[2] * ((cos(pitch) * sin(yaw) * cos(roll)) + (sin(pitch) * sin(roll))));
    out[1] =
        (in[0] * sin(roll) * cos(yaw))
      + (in[1] * ((sin(pitch) * sin(yaw) * sin(roll)) + (cos(pitch) * cos(roll))))
      + (in[2] * ((sin(yaw) * sin(roll) * cos(pitch)) - (cos(roll) * sin(pitch))));
    out[2] =
        (in[0] * -sin(yaw))
      + (in[1] * cos(yaw) * sin(pitch))
      + (in[2] * cos(pitch) * cos(yaw));
}

double time_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + (ts.tv_nsec / 1.0e9);
}

#define MAX_VERTS 10000
double verts[MAX_VERTS][3];
int n_verts = 0;
#define MAX_TRIS 10000
int tris[MAX_TRIS][3];
int n_tris = 0;

void load_obj(char *fn, double scale) {
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
                ASSERT(sscanf(token, "%lf", &verts[n_verts][off]) == 1);
                verts[n_verts][off] *= scale;
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

int main(void) {
    load_obj("cow-nonormals.obj", 1);

    ASSERT(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) == 0);
    SDL_Window *sdl_window = SDL_CreateWindow("aspng", 0, 0, 100, 100, SDL_WINDOW_RESIZABLE);
    ASSERT(sdl_window != NULL);
    SDL_Renderer *sdl_renderer = SDL_CreateRenderer(sdl_window, -1, SDL_RENDERER_SOFTWARE);
    ASSERT(sdl_renderer != NULL);

    int did_print = 0;
    double *depth_buffer = NULL;
    int show_depth = 0;
    int free_flight = 0;

    double camera_pos[3] = {0, 0, -10};
    double camera_yaw = 0;
    double camera_pitch = 0;
    double camera_roll = 0;
    double since = time_now();
    int frames_since = 0;

    int do_quit = 0;
    while (!do_quit) {
        double now = time_now();
        if (now - since >= 1) {
            printf("fps: %f\n", frames_since / (now - since));
            since = now;
            frames_since = 0;
        }
        frames_since++;

        double forward[3] = {0, 0, 1};
        double camera_fwd[3];
        rotate(forward, camera_yaw, camera_pitch, camera_roll, camera_fwd);
        double up[3] = {0, 1, 0};
        double camera_up[3];
        rotate(up, camera_yaw, camera_pitch, camera_roll, camera_up);
        double right[3] = {1, 0, 0};
        double camera_right[3];
        rotate(right, camera_yaw, camera_pitch, camera_roll, camera_right);

        SDL_Surface *window_surface = SDL_GetWindowSurface(sdl_window);
        if (depth_buffer == NULL) {
            depth_buffer = malloc(sizeof(double) * window_surface->w * window_surface->h);
        }

        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            switch (e.type) {
            case SDL_KEYDOWN:
                switch (e.key.keysym.sym) {
                case SDLK_ESCAPE:
                    do_quit = 1;
                    break;
                case SDLK_d:
                    show_depth = 1 - show_depth;
                    break;
                case SDLK_f:
                    free_flight = 1 - free_flight;
                    break;
                }
                break;
            case SDL_WINDOWEVENT:
                switch (e.window.event) {
                case SDL_WINDOWEVENT_RESIZED:
                    free(depth_buffer);
                    window_surface = SDL_GetWindowSurface(sdl_window);
                    depth_buffer = malloc(sizeof(double) * window_surface->w * window_surface->h);
                    break;
                }
            }
            did_print = 0;
        }

        const uint8_t *keystate = SDL_GetKeyboardState(NULL);
        if (free_flight) {
            if (keystate[SDL_SCANCODE_LEFT])
                camera_yaw -= 0.01;
            if (keystate[SDL_SCANCODE_RIGHT])
                camera_yaw += 0.01;
            if (keystate[SDL_SCANCODE_UP])
                camera_pitch += 0.01;
            if (keystate[SDL_SCANCODE_DOWN])
                camera_pitch -= 0.01;
            if (keystate[SDL_SCANCODE_Z])
                camera_roll -= 0.01;
            if (keystate[SDL_SCANCODE_X])
                camera_roll += 0.01;
            if (keystate[SDL_SCANCODE_A]) {
                camera_pos[0] -= camera_fwd[0] * 0.1;
                camera_pos[1] -= camera_fwd[1] * 0.1;
                camera_pos[2] -= camera_fwd[2] * 0.1;
            }
            if (keystate[SDL_SCANCODE_S]) {
                camera_pos[0] += camera_fwd[0] * 0.1;
                camera_pos[1] += camera_fwd[1] * 0.1;
                camera_pos[2] += camera_fwd[2] * 0.1;
            }
        } else {
            // Rotate the camera around the origin.
            double now = time_now();
            camera_pos[0] = 10 * sin(5 * now / (2 * M_PI));
            camera_pos[1] = 0;
            camera_pos[2] = 10 * cos(5 * now / (2 * M_PI));

            // Point the camera at the origin.
            camera_yaw = (2 * M_PI) - fmod(M_PI + (5 * now / (2 * M_PI)), 2 * M_PI);
            camera_pitch = 0;
            camera_roll = M_PI;
        }

        // TODO clear the screen: shouldn't have to do this out of band
        for (int x = 0; x < window_surface->w; x++) {
            for (int y = 0; y < window_surface->h; y++) {
                uint32_t pixel = SDL_MapRGBA(window_surface->format, 0, 0, 0, 0xff);
                int off = (y * window_surface->w) + x;
                ((uint32_t *)window_surface->pixels)[off] = pixel;
                depth_buffer[off] = DBL_MAX;
            }
        }

        num_triangles = 0;
        for (int i = 0; i < n_tris; i++) {
            int projectible = 1;
            for (int j = 0; j < 3; j++) {
                double v[3] = {
                    verts[tris[i][j]][0],
                    verts[tris[i][j]][1],
                    verts[tris[i][j]][2]
                };

                // Project into screen space: draw a ray from the vertex to the
                // camera, intersecting with a plane (the screen).
                double screen_center[3] = {
                    camera_pos[0] + camera_fwd[0],
                    camera_pos[1] + camera_fwd[1],
                    camera_pos[2] + camera_fwd[2]
                };
                double camera_pos_to_v[3] = {
                    v[0] - camera_pos[0],
                    v[1] - camera_pos[1],
                    v[2] - camera_pos[2],
                };
                double poi[3];
                if (ray_plane(screen_center, camera_fwd, camera_pos, camera_pos_to_v, poi)) {
                    // Find the intersection, relative to the plane center.
                    double poi_rel[3] = {
                        poi[0] - screen_center[0],
                        poi[1] - screen_center[1],
                        poi[2] - screen_center[2]
                    };
                    // Put this vertex into screenspace.
                    //double plane_w = 2 * 1 /* camera_fwd is unit */ * tan(camera_fov / 2);
                    double screen_w = window_surface->w, screen_h = window_surface->h;
                    //double plane_h = plane_w * (screen_h / screen_w);
                    triangles[num_triangles].v[j].x =
                        (screen_w / 2) + (DOT(poi_rel, camera_right) * screen_w);
                    triangles[num_triangles].v[j].y =
                        (screen_h / 2) + (DOT(poi_rel, camera_up) * screen_h);
                    // If the vertex is behind the camera, then the distance
                    // should be negative.
                    triangles[num_triangles].v[j].z = 
                          DOT(camera_pos_to_v, camera_fwd) >= 0
                        ? MAGNITUDE(camera_pos_to_v)
                        : -MAGNITUDE(camera_pos_to_v);
                } else {
                    projectible = 0;
                }
            }
            if (projectible) {
                triangles[num_triangles].r = (i * 13) % 0xff;
                triangles[num_triangles].g = (i * 101) % 0xff;
                triangles[num_triangles].b = (i * 211) % 0xff;
                ASSERT(++num_triangles < NUM_TRIANGLES);
            }
        }

        num_spans = 0;
        num_pointup_spans = 0;
        num_pointdown_spans = 0;
        num_degenerate_spans = 0;
        num_onespan_triangles = 0;
        num_twospan_triangles = 0;
        num_degenerate_triangles = 0;
        for (int i = 0; i < num_triangles; i++) {
            screen_vertex *a = &triangles[i].v[0];
            screen_vertex *b = &triangles[i].v[1];
            screen_vertex *c = &triangles[i].v[2];
            // Look for (one) top ("hi") vertex.
            screen_vertex *hi = NULL;
            if ((a->y > b->y) && (a->y > c->y)) {
                hi = a;
            } else if ((b->y > a->y) && (b->y > c->y)) {
                hi = b;
            } else if ((c->y > a->y) && (c->y > b->y)) {
                hi = c;
            }
            // Look for (one) bottom ("lo") vertex.
            screen_vertex *lo = NULL;
            if ((a->y < b->y) && (a->y < c->y)) {
                lo = a;
            } else if ((b->y < a->y) && (b->y < c->y)) {
                lo = b;
            } else if ((c->y < a->y) && (c->y < b->y)) {
                lo = c;
            }
            // If there's neither a hi nor lo vertex, then it's
            // degenerate.
            if ((hi == NULL) && (lo == NULL)) {
                num_degenerate_triangles++;
                continue;
            }
            // If there's only a hi or lo vertex, then there is only one
            // span.
            if ((hi == NULL) != (lo == NULL)) {
                if (hi != NULL) {
                    add_span(a, b, c, hi, NULL, &triangles[i]);
                } else if (lo != NULL) {
                    add_span(a, b, c, NULL, lo, &triangles[i]);
                } else ASSERT(0);
                num_onespan_triangles++;
            }
            // If there is both a hi and lo vertex, then we need to draw
            // two spans.
            if ((hi != NULL) && (lo != NULL)) {
                // First, we need to find the vertex ('mid') which isn't the hi
                // or lo vertex.
                screen_vertex *mid = NULL;
                if ((a != hi) && (a != lo)) mid = a;
                else if ((b != hi) && (b != lo)) mid = b;
                else mid = c;
                ASSERT(mid != NULL);
                // Find the point on the edge linking hi and lo which is at the
                // same y-coordinate as 'mid'.
                screen_vertex split = {
                    .x = lo->x + (((hi->x - lo->x) / (double)(hi->y - lo->y)) * (mid->y - lo->y)),
                    .y = mid->y,
                    .z = lo->z + (((hi->z - lo->z) / (double)(hi->y - lo->y)) * (mid->y - lo->y))
                };
                // Create two spans!
                add_span(hi, mid, &split, hi, NULL, &triangles[i]);
                add_span(lo, mid, &split, NULL, lo, &triangles[i]);
                num_twospan_triangles++;
            }
        }

        if (!did_print) {
            did_print = 1;
            printf("******\n");
            printf("yaw = %f, pitch = %f, roll = %f\n", camera_yaw, camera_pitch, camera_roll);
            printf("camera: %f, %f, %f\n", camera_pos[0], camera_pos[1], camera_pos[2]);
            printf("forward: %f, %f, %f .. %f\n", camera_fwd[0], camera_fwd[1], camera_fwd[2], MAGNITUDE(camera_fwd));
            printf("up: %f, %f, %f .. %f\n", camera_up[0], camera_up[1], camera_up[2], MAGNITUDE(camera_up));
            printf("right: %f, %f, %f .. %f\n", camera_right[0], camera_right[1], camera_right[2], MAGNITUDE(camera_right));
            printf("%i triangles (%i onespans, %i twospans, %i degenerates)\n",
                num_triangles, num_onespan_triangles, num_twospan_triangles, num_degenerate_triangles);
            printf("%i spans (%i pointups, %i pointdowns, %i degenerates)\n",
                num_spans, num_pointup_spans, num_pointdown_spans, num_degenerate_spans);
        }

        double min_z = DBL_MAX;
        double max_z = -DBL_MAX;
        for (int y = 0; y < window_surface->h; y++) {
            for (int i = 0; i < num_spans; i++) {
                span_t *span = &spans[i];
                if ((span->y_lo <= y) && (y <= span->y_hi)) {
                    int16_t x_fill_lo = span->ref.x + (span->dx_dy_lo * (y - span->ref.y));
                    int16_t x_fill_hi = span->ref.x + (span->dx_dy_hi * (y - span->ref.y));
                    ASSERT(x_fill_lo <= x_fill_hi);
                    if (x_fill_lo < 0)
                        x_fill_lo = 0;
                    if (x_fill_hi > window_surface->w - 1)
                        x_fill_hi = window_surface->w - 1;
                    double z_lo = span->ref.z + (span->dz_dy_lo * (y - span->ref.y));
                    for (int16_t x = x_fill_lo; x <= x_fill_hi; x++) {
                        // Do a z-check before we draw the pixel.
                        int off = (y * window_surface->w) + x;
                        double z = z_lo + (span->dz_dx_lo * (x - x_fill_lo));
                        if ((z < depth_buffer[off]) && (z >= 0)) {
                            depth_buffer[off] = z;
                            if (z > max_z) {
                                max_z = z;
                            }
                            if (z < min_z) {
                                min_z = z;
                            }
                            uint32_t pixel = SDL_MapRGBA(
                                window_surface->format,
                                span->parent->r,
                                span->parent->g,
                                span->parent->b,
                                0xff
                            );
                            ((uint32_t *)window_surface->pixels)[off] = pixel;
                        }
                    }
                }
            }
        }

        if (show_depth) {
            for (int y = 0; y < window_surface->h; y++) {
                for (int x = 0; x < window_surface->w; x++) {
                    double z = depth_buffer[(y * window_surface->w) + x];
                    double ramped = 255 * (1 - ((z - min_z) / (max_z - min_z)));
                    uint32_t pixel = SDL_MapRGBA(window_surface->format, 0, ramped, 0, 0xff);
                    ((uint32_t *)window_surface->pixels)[(y * window_surface->w) + x] = pixel;
                }
            }
        }

        SDL_UpdateWindowSurface(sdl_window);
    }
}
