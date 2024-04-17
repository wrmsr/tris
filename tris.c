#include <GL/gl.h>
#include <SDL2/SDL.h>
#include <stdbool.h>
#include <stdint.h>

#include "common.h"
#include "render.h"
#include "vector.h"

extern render_stats_t render_stats;
extern render_frame_stats_t render_frame_stats;

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
                    v3_rotate(&verts[n_verts], yaw, pitch, roll, &verts[n_verts]);
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

// https://gamedev.stackexchange.com/a/23745
void barycentric(v3_t *p, triangle_t *triangle, double *b1, double *b2, double *b3) {
    v3_t v0;
    /*printf("point: {%f, %f, %f}\n",
        p->x, p->y, p->z);
    printf("triangle: {%f, %f, %f} , {%f, %f, %f} , {%f, %f, %f}\n",
        triangle->b.xyz.x, triangle->b.xyz.y, triangle->b.xyz.z,
        triangle->c.xyz.x, triangle->c.xyz.y, triangle->c.xyz.z,
        triangle->a.xyz.x, triangle->a.xyz.y, triangle->a.xyz.z );*/
    v3_sub(&(triangle->b.xyz), &(triangle->a.xyz), &v0);
    //printf("v0 = {%f, %f, %f}\n", v0.x, v0.y, v0.z);
    v3_t v1;
    v3_sub(&(triangle->c.xyz), &(triangle->a.xyz), &v1);
    v3_t v2;
    v3_sub(p, &(triangle->a.xyz), &v2);
    double d00 = v3_dot(&v0, &v0);
    double d01 = v3_dot(&v0, &v1);
    double d11 = v3_dot(&v1, &v1);
    double d20 = v3_dot(&v2, &v0);
    double d21 = v3_dot(&v2, &v1);
    //printf("d00 d01 d11 %f %f %f\n", d00, d01, d11);
    double denom = (d00 * d11) - (d01 * d01);
    *b1 = ((d11 * d20) - (d01 * d21)) / denom;
    //printf("%f, %f\n", denom, *b1);
    //ASSERT((*b1 >= 0) && (*b1 <= 1));
    *b2 = ((d00 * d21) - (d01 * d20)) / denom;
    //ASSERT((*b2 >= 0) && (*b2 <= 1));
    *b3 = 1. - *b1 - *b2;
    //ASSERT((*b3 >= 0) && (*b3 <= 1));
    //ASSERT(*b1 + *b2 + *b3 == 1);
}

double *depth_buffer = NULL;
gl_rgb_t *texture = NULL;

int main(void) {
    load_obj("cow-nonormals.obj", 1, 0, 0, -M_PI / 2);
    material_t material = {
        .w = 64,
        .h = 64,
        .texture = malloc(sizeof(rgba_t) * 64 * 64),
    };
    for (int y = 0; y < 64; y++) {
        for (int x = 0; x < 64; x++) {
            rgba_t *rgba = &(material.texture[(y * 64) + x]);
            rgba->r = (x % 16) * 16;
            rgba->g = (y % 16) * 16;
            rgba->b = 0;
            rgba->a = 0;
        }
    }

    ASSERT(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) == 0);
    SDL_Window *sdl_window = SDL_CreateWindow("aspng", 0, 0, 100, 100, SDL_WINDOW_RESIZABLE | SDL_WINDOW_OPENGL);
    ASSERT(sdl_window != NULL);
    ASSERT(SDL_GL_CreateContext(sdl_window) != NULL);
    GLuint screen_texture;
    glGenTextures(1, &screen_texture);

    depth_buffer = malloc(sizeof(double) * FAKESCREEN_W * FAKESCREEN_H);
    texture = malloc(sizeof(gl_rgb_t) * FAKESCREEN_W * FAKESCREEN_H);

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
    int do_quit = 0;
    while (!do_quit && render_stats.frames_drawn < 200) {
        render_frame_stats.start_time = time_now();
        render_frame_stats.pixels_rejected_by_z = 0;
        render_frame_stats.pixels_drawn = 0;

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

        render_begin_frame();
        for (int i = 0; i < n_tris; i++) {
            // Build the triangle to be rendered.
            triangle_t *triangle = render_add_triangle();
            for (int j = 0; j < 3; j++) {
                memcpy(&(triangle->abc[j].xyz), &verts[tris[i][j]], sizeof(v3_t));
                triangle->abc[j].uv.u = 0;
                triangle->abc[j].uv.v = 0;
                triangle->material = &material;
            }
            render_draw_triangle(&camera_pos, &camera_fwd, &camera_up, &camera_left, triangle);
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

        render_end_frame();
    }

    render_print_stats();
}
