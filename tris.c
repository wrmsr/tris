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
            rgba->r = ((x / 16) % 2) * 255;
            rgba->g = 0;
            rgba->b = ((y / 16) % 2) * 255;
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

    v3_t camera_pos = { .x = -10, .y = 0, .z = 0 }; // TODO this initialization doesnt matter
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
    while (!do_quit) { // && render_stats.frames_drawn < 400) {
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

        // TODO clear the screen: shouldn't have to do this out of band
        for (int y = 0; y < FAKESCREEN_H; y++) {
            for (int x = 0; x < FAKESCREEN_W; x++) {
                int off = (y * FAKESCREEN_W) + x;
                depth_buffer[off] = DBL_MAX;
                texture[off].r = texture[off].g = texture[off].b = 0;
            }
        }

        enum {
            DO_WHAT_COW,
            DO_WHAT_CUBE,
        } do_what;
        if (((int)time_now() / 5) % 2 == 0) {
            do_what = DO_WHAT_CUBE;
        } else {
            do_what = DO_WHAT_COW;
        }

        if (do_what == DO_WHAT_COW) {
            // Rotate the camera around the origin.
            double now = render_stats.frames_drawn / 10.;
            double theta = 0.25 * now;
            double scale = 3 + (.5 * cos(2 * now));
            camera_pos.x = 2. * cos(theta) * scale;
            camera_pos.y = 2. * sin(theta) * scale;
            camera_pos.z = 2. * sin(theta) * scale;
        } else {
            double now = time_now() * 2; //render_stats.frames_drawn / 10.;
            double theta = 0.25 * now;
            double scale = 4; //3 + (.5 * cos(2 * now));
            camera_pos.x = 2. * cos(theta) * scale;
            camera_pos.y = 0; //2. * sin(theta) * scale;
            camera_pos.z = 2. * sin(theta) * scale;
        }

        // Point the camera at the origin.
        camera_fwd.x = 0 - camera_pos.x;
        camera_fwd.y = 0 - camera_pos.y;
        camera_fwd.z = 0 - camera_pos.z;
        double l = v3_len(&camera_fwd);
        camera_fwd.x = camera_fwd.x / l;
        camera_fwd.y = camera_fwd.y / l;
        camera_fwd.z = camera_fwd.z / l;
        camera_up.x = 0;
        camera_up.y = 1;
        camera_up.z = 0;
        v3_cross(&camera_fwd, &camera_up, &camera_left);


        render_begin_frame();
        if (do_what == DO_WHAT_CUBE) {
            v3_t coords[] = {
                { .x = -2, .y = -2, .z = -2 },
                { .x =  2, .y = -2, .z = -2 },
                { .x =  2, .y =  2, .z = -2 },
                { .x = -2, .y =  2, .z = -2 },
                { .x = -2, .y = -2, .z =  2 },
                { .x =  2, .y = -2, .z =  2 },
                { .x =  2, .y =  2, .z =  2 },
                { .x = -2, .y =  2, .z =  2 },
            };
            triangle_t triangles[] = {
                {
                    .a = { .xyz = coords[2], .uv = { 0, 0 } },
                    .b = { .xyz = coords[1], .uv = { 1, 0 } },
                    .c = { .xyz = coords[0], .uv = { 0, 1 } },
                    .material = &material
                }/*, {
                    .a = { .xyz = coords[0], .uv = { 0, 0 } },
                    .b = { .xyz = coords[3], .uv = { 0, 1 } },
                    .c = { .xyz = coords[2], .uv = { 1, 0 } },
                    .material = &material
                },*/
                /*
                { .a = coords[1], .b = coords[2], .c = coords[6], .material = &material },
                { .a = coords[6], .b = coords[5], .c = coords[1], .material = &material },
                { .a = coords[2], .b = coords[3], .c = coords[7], .material = &material },
                { .a = coords[7], .b = coords[6], .c = coords[2], .material = &material },
                { .a = coords[3], .b = coords[0], .c = coords[4], .material = &material },
                { .a = coords[4], .b = coords[7], .c = coords[3], .material = &material },
                { .a = coords[0], .b = coords[1], .c = coords[5], .material = &material },
                { .a = coords[5], .b = coords[4], .c = coords[0], .material = &material },
                { .a = coords[4], .b = coords[5], .c = coords[6], .material = &material },
                { .a = coords[6], .b = coords[7], .c = coords[4], .material = &material },*/
            };
            for (int i = 0; i < sizeof(triangles) / sizeof(triangles[0]); i++) {
                render_draw_triangle(&camera_pos, &camera_fwd, &camera_up, &camera_left, &triangles[i]);
            }
        } else {
            for (int i = 0; i < n_tris; i++) {
                // Build the triangle to be rendered.
                triangle_t triangle;
                for (int j = 0; j < 3; j++) {
                    memcpy(&(triangle.abc[j].xyz), &verts[tris[i][j]], sizeof(v3_t));
                    triangle.abc[j].u = 0;
                    triangle.abc[j].v = 0;
                    triangle.material = NULL;
                }
                render_draw_triangle(&camera_pos, &camera_fwd, &camera_up, &camera_left, &triangle);
            }
        }

        for (int y = 0; y < 64; y++) {
            for (int x = 0; x < 64; x++) {
                gl_rgb_t *pixel = &texture[(y * FAKESCREEN_W) + x];
                rgba_t *texel = &(material.texture[(y * 64) + x]);
                pixel->r = texel->r;
                pixel->g = texel->g;
                pixel->b = texel->b;
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

        render_end_frame();
    }

    render_print_stats();
}
