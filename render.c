#include <float.h>
#include <stddef.h>
#include <string.h>

#include "common.h"
#include "render.h"

render_stats_t render_stats = {
    .start_time = -1,
    .frames_drawn = 0,
    .pixels_drawn = 0,
    .pixels_rejected_by_z = 0,
};

render_frame_stats_t render_frame_stats;
int color_index;

void render_begin_frame(void) {
    color_index = 0;

    if (render_stats.start_time < 0) {
        render_stats.start_time = time_now();
    }

    render_frame_stats.n_skipped_triangles = 0;
    render_frame_stats.n_skipped_triangles_by_backface_cull = 0;
    render_frame_stats.n_skipped_triangles_by_near_plane = 0;
    render_frame_stats.num_spans = 0;
    render_frame_stats.num_pointup_spans = 0;
    render_frame_stats.num_pointdown_spans = 0;
    render_frame_stats.num_degenerate_spans = 0;
    render_frame_stats.num_triangles = 0;
    render_frame_stats.num_onespan_triangles = 0;
    render_frame_stats.num_twospan_triangles = 0;
    render_frame_stats.num_degenerate_triangles = 0;
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
    v3_t v1;
    v3_sub(&(triangle->c.xyz), &(triangle->a.xyz), &v1);
    v3_t v2;
    v3_sub(p, &(triangle->a.xyz), &v2);
    double d00 = v3_dot(&v0, &v0);
    double d01 = v3_dot(&v0, &v1);
    double d11 = v3_dot(&v1, &v1);
    double d20 = v3_dot(&v2, &v0);
    double d21 = v3_dot(&v2, &v1);
    double denom = (d00 * d11) - (d01 * d01);
    *b1 = ((d11 * d20) - (d01 * d21)) / denom;
    *b2 = ((d00 * d21) - (d01 * d20)) / denom;
    *b3 = 1. - *b1 - *b2;
    //printf("b1 = %f, b2 = %f, b3 = %f\n", *b1, *b2, *b3);
    //ASSERT((*b1 >= 0) && (*b1 <= 1));
    //ASSERT((*b2 >= 0) && (*b2 <= 1));
    //ASSERT((*b3 >= 0) && (*b3 <= 1)); // need to solve precision issues first
    //ASSERT(*b1 + *b2 + *b3 == 1);
}

extern double *depth_buffer;
extern gl_rgb_t *texture;

void draw_span(screen_vertex_t *a, screen_vertex_t *b, screen_vertex_t *c, screen_vertex_t *hi, screen_vertex_t *lo, screen_triangle_t *screen_triangle) {
    render_frame_stats.num_spans++;
    ASSERT((hi == NULL) != (lo == NULL));
    screen_vertex_t *hi_or_lo = (hi != NULL ? hi : lo);
    // Of the other two vertices, which is on the left, and which is on the
    // right?
    span_t span;
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
        render_frame_stats.num_degenerate_spans++;
        return;
    }
    if (hi != NULL) {
        render_frame_stats.num_pointdown_spans++;
        span.y_lo = MAX(0, MIN(FAKESCREEN_H, other1->y)); // or other2[1], doesn't matter.
        span.y_hi = MAX(0, MIN(FAKESCREEN_H, hi->y));
    } else {
        render_frame_stats.num_pointup_spans++;
        span.y_lo = MAX(0, MIN(FAKESCREEN_H, lo->y));
        span.y_hi = MAX(0, MIN(FAKESCREEN_H, other1->y)); // or other2[1], doesn't matter.
    }
    memcpy(&(span.ref), hi_or_lo, sizeof(span.ref));
    // Of the other two vertices, which is more distant in x from the reference
    // point?
    screen_vertex_t *x_diff_vert;
    if (abs(span.ref.x - other1->x) > abs(span.ref.x - other2->x)) {
        x_diff_vert = other1;
    } else {
        x_diff_vert = other2;
    }

    // The dx/dy slope for drawing the left side of the triangle:
    float dsx_dsy_lo = (span.ref.x - x_lo_vert->x) / (double)(span.ref.y - x_lo_vert->y);
    // The dx/dy slope for drawing the right side of the triangle:
    float dsx_dsy_hi = (span.ref.x - x_hi_vert->x) / (double)(span.ref.y - x_hi_vert->y);
    // Slopes for conversion of screenspace x and y to objectspace x, y, and z.
    /*float dox_dsx = (span.ref.object.x - x_diff_vert->object.x) / (double)(span.ref.x - x_diff_vert->x);
    float dox_dsy = (span.ref.object.x - x_diff_vert->object.x) / (double)(span.ref.y - x_diff_vert->y);
    float doy_dsx = (span.ref.object.y - x_diff_vert->object.y) / (double)(span.ref.x - x_diff_vert->x);
    float doy_dsy = (span.ref.object.y - x_diff_vert->object.y) / (double)(span.ref.y - x_diff_vert->y);
    */
    //float dsz_dsx = (span.ref.object.z - x_diff_vert->object.z) / (double)(span.ref.x - x_diff_vert->x);
    //float dsz_dsy = (span.ref.object.z - x_diff_vert->object.z) / (double)(span.ref.y - x_diff_vert->y);
    float dsz_dsx = (span.ref.z - x_diff_vert->z) / (double)(span.ref.x - x_diff_vert->x);
    float dsz_dsy = (span.ref.z - x_diff_vert->z) / (double)(span.ref.y - x_diff_vert->y);
    span.triangle = screen_triangle->parent;

    // Draw the spans to the screen, respecting the z-buffer.
    double min_z = DBL_MAX;
    double max_z = -DBL_MAX;
    for (int y = span.y_lo; y < span.y_hi; y++) {
                int x_fill_lo = span.ref.x + (dsx_dsy_lo * (y - span.ref.y)); // TODO; this can just be +1 since we're going row by row, no mult needed
                int x_fill_hi = span.ref.x + (dsx_dsy_hi * (y - span.ref.y));
                ASSERT(x_fill_lo <= x_fill_hi);
                if (x_fill_lo < 0)
                    x_fill_lo = 0;
                if (x_fill_hi > FAKESCREEN_W - 1)
                    x_fill_hi = FAKESCREEN_W - 1;
                // This is a lot faster than doing the barymetric transform
                // now, given there's a chance this pixel gets z-rejected.
                double screen_z_lo = span.ref.z + (dsz_dsy * (y - span.ref.y));
                for (int x = x_fill_lo; x <= x_fill_hi; x++) {

                    // Do a z-check before we draw the pixel.
                    int off = (y * FAKESCREEN_W) + x;
                    double screen_z = screen_z_lo + (dsz_dsx * (x - x_fill_lo));
                    if ((screen_z < depth_buffer[off]) && (screen_z >= 0)) {
                        depth_buffer[off] = screen_z;
                        if (screen_z > max_z) {
                            max_z = screen_z;
                        }
                        if (screen_z < min_z) {
                            min_z = screen_z;
                        }
                        v3_t p = { .x = x, .y = y, .z = 0 };
                        double b1, b2, b3;
                        triangle_t triangle = {
                            .a = { .xyz = {.x = screen_triangle->v[0].x, .y = screen_triangle->v[0].y, .z = screen_triangle->v[0].z} },
                            .b = { .xyz = {.x = screen_triangle->v[1].x, .y = screen_triangle->v[1].y, .z = screen_triangle->v[1].z} },
                            .c = { .xyz = {.x = screen_triangle->v[2].x, .y = screen_triangle->v[2].y, .z = screen_triangle->v[2].z} },
                        };
                        barycentric(&p, &triangle, &b1, &b2, &b3);
                        //double object_x = (b1 * span.ref.object.x) + (b2 * span.ref.object.y) + (b3 * span.ref.object.z);
                    //double object_y = (b1 * span.ref.object.y) + (b2 * span.ref.object.y) + (b3 * span.ref.object.y);
                    //double screen_z = (b1 * x_lo_vert->z) + (b2 * x_hi_vert->z) + (b3 * hi_or_lo->z);
                    //printf("%f, %f, %f\n", span.triangle->a.u, span.triangle->b.u, b3 * span.triangle->c.u);
                    double u = (b1 * span.triangle->a.u) + (b2 * span.triangle->b.u) + (b3 * span.triangle->c.u);
                    //printf("%f, %f, %f\n", span.triangle->a.v, span.triangle->b.v, b3 * span.triangle->c.v);
                    double v = (b1 * span.triangle->a.v) + (b2 * span.triangle->b.v) + (b3 * span.triangle->c.v);

                        /*printf("dox_dsx = %f\n", dox_dsx);
                        printf("dox_dsy = %f\n", dox_dsy);
                        printf("doy_dsx = %f\n", doy_dsx);
                        printf("doy_dsy = %f\n", doy_dsy);
                        printf("dsz_dsx = %f\n", dsz_dsx);
                        printf("dsz_dsy = %f\n", dsz_dsy);
                        printf("triangle in screen: (%i, %i) (%i, %i) (%i, %i)\n",
                            x_lo_vert->x, x_lo_vert->y, hi_or_lo->x, hi_or_lo->y, x_hi_vert->x, x_hi_vert->y);
                        printf("current screen x, y = %i, %i\n", x, y);
                        printf("x_fill lo = %i, hi = %i\n", x_fill_lo, x_fill_hi);
                        */
                                                //double object_x = (b1 * span.ref.object.x) + (b2 * span.ref.object.y) + (b3 * span.ref.object.z);
                        //double object_y = (b1 * span.ref.object.y) + (b2 * span.ref.object.y) + (b3 * span.ref.object.y);
                        //double object_z = (b1 * span.ref.object.z) + (b2 * span.ref.object.z) + (b3 * span.ref.object.z);
                        //double u = (b1 * span.triangle->a.uv.u) + (b2 * span.triangle->b.uv.u) + (b3 * span.triangle->c.uv.u);
                        // TODO assert?
                        if (span.triangle->material != NULL) {
                            u = MAX(0, MIN(1, u));
                            int x_tex = u * span.triangle->material->w;
                            //double v = (b1 * span.triangle->a.uv.v) + (b2 * span.triangle->b.uv.v) + (b3 * span.triangle->c.uv.v);
                            v = MAX(0, MIN(1, v));
                            //printf("u, v = %f, %f\n", u, v);
                            int y_tex = v * span.triangle->material->h;
                            int tex_off = x_tex + (y_tex * span.triangle->material->w);
                            texture[off].r = span.triangle->material->texture[tex_off].r;
                            texture[off].g = span.triangle->material->texture[tex_off].g;
                            texture[off].b = span.triangle->material->texture[tex_off].b;
                        
                        } else {
                            texture[off].r = (color_index * 13) % 0xff;
                            texture[off].g = (color_index * 101) % 0xff;
                            texture[off].b = (color_index * 211) % 0xff;
                        }

                        render_stats.pixels_drawn++;
                        render_frame_stats.pixels_drawn++;
                    } else {
                        render_stats.pixels_rejected_by_z++;
                        render_frame_stats.pixels_rejected_by_z++;
                    }
                }
    }
}

void render_draw_triangle(v3_t *camera_pos, v3_t *camera_fwd, v3_t *camera_up, v3_t *camera_left, triangle_t *triangle) {
    color_index++;

    // The screen exists on a plane "in front of" the camera.
    v3_t screen_center;
    v3_add(camera_pos, camera_fwd, &screen_center);

    // Do backface culling: a triangle facing the wrong way from the camera
    // (according to its normal) doesn't get rendered.
    v3_t v1, v2;
    v3_sub(&(triangle->abc[0].xyz), &(triangle->abc[1].xyz), &v1);
    v3_sub(&(triangle->abc[0].xyz), &(triangle->abc[2].xyz), &v2);
    v3_t normal;
    v3_cross(&v2, &v1, &normal);
    v3_t cam_to_triangle;
    v3_sub(&(triangle->abc[0].xyz), camera_pos, &cam_to_triangle);
    if (v3_dot(&cam_to_triangle, &normal) < 0) {
        render_frame_stats.n_skipped_triangles++;
        render_frame_stats.n_skipped_triangles_by_backface_cull++;
        return;
    }

    // Turn each arbitrary triangle into triangles suitable for rendering, by
    // clipping away the parts that cannot be correctly rendered.
    int n_inside = 0, inside[3], n_outside = 0, outside[3];
    for (int j = 0; j < 3; j++) {
        // Calculate signed distance to figure out which side of the screen
        // plane the vertex is on.
        v3_t screen_to_v;
        v3_sub(&(triangle->abc[j].xyz), &screen_center, &screen_to_v);
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
        // If one point is inside (v3), then when we clip the triangle it
        // becomes a quadrilateral. The quadrilateral has two of the original
        // vertices, and an additional 2 vertices (named u1, u2) at the point
        // where (v3, v1) and (v3, v2) intersect the screen plane.
        v3_t *v3 = &(triangle->abc[inside[0]].xyz);
        v3_t *v1 = &(triangle->abc[outside[0]].xyz);
        v3_t v1_to_v3;
        v3_sub(v3, v1, &v1_to_v3);
        v3_t u1;
        // v1 and v3 are on opposite sides of the camera plane, so v3-v1 and
        // the camera plane cannot be parallel.
        ASSERT(v3_ray_plane(&screen_center, camera_fwd, v1, &v1_to_v3, &u1));
        v3_t *v2 = &(triangle->abc[outside[1]].xyz);
        v3_t v2_to_v3;
        v3_sub(v3, v2, &v2_to_v3);
        v3_t u2;
        // Ditto, see above.
        ASSERT(v3_ray_plane(&screen_center, camera_fwd, v2, &v2_to_v3, &u2));
        // We can't render a quadrilateral, but we can split it into
        // two triangles. Since we know that the perimeter of the
        // quadrilateral is formed by visiting v1, v2, u2, u1 in that
        // order, we know that we can form two triangles (v1, v2, u2)
        // and (u2, u1, v1).
        n_clipped_triangles = 2;
        memcpy(&(clipped_triangles[0].abc[0].xyz), v1, sizeof(v3_t));
        memcpy(&(clipped_triangles[0].abc[1].xyz), v2, sizeof(v3_t));
        memcpy(&(clipped_triangles[0].abc[2].xyz), &u2, sizeof(v3_t));
        memcpy(&(clipped_triangles[1].abc[0].xyz), &u2, sizeof(v3_t));
        memcpy(&(clipped_triangles[1].abc[1].xyz), &u1, sizeof(v3_t));
        memcpy(&(clipped_triangles[1].abc[2].xyz), v1, sizeof(v3_t));
    } else if (n_inside == 2) {
        ASSERT(n_outside == 1);
        // if two points inside (v1, v2), then make one triangle (v3, u1, u2)
        v3_t *v3 = &(triangle->abc[outside[0]].xyz);
        v3_t *v1 = &(triangle->abc[inside[0]].xyz);
        v3_t v1_to_v3;
        v3_sub(v3, v1, &v1_to_v3);
        v3_t u1;
        // v1 and v3 are on opposite sides of the camera plane, so v3-v1 and
        // the camera plane cannot be parallel.
        ASSERT(v3_ray_plane(&screen_center, camera_fwd, v1, &v1_to_v3, &u1));
        v3_t *v2 = &(triangle->abc[inside[1]].xyz);
        v3_t v2_to_v3;
        v3_sub(v3, v2, &v2_to_v3);
        v3_t u2;
        // Ditto, see above.
        ASSERT(v3_ray_plane(&screen_center, camera_fwd, v2, &v2_to_v3, &u2));
        n_clipped_triangles = 1;
        memcpy(&(clipped_triangles[0].abc[0].xyz), v3, sizeof(v3_t));
        memcpy(&(clipped_triangles[0].abc[1].xyz), &u1, sizeof(v3_t));
        memcpy(&(clipped_triangles[0].abc[2].xyz), &u2, sizeof(v3_t));
    } else if (n_inside == 3) {
        // Just don't render this triangle.
        ASSERT(n_outside == 0);
        render_frame_stats.n_skipped_triangles++;
        render_frame_stats.n_skipped_triangles_by_near_plane++;
    } else {
        ASSERT(0);
    }
    ASSERT((n_clipped_triangles >= 0) && (n_clipped_triangles <= 2));

    // Project each clipped triangle into screen space.
    for (int j = 0; j < n_clipped_triangles; j++) {
        screen_triangle_t screen_triangle;
        screen_triangle.parent = triangle;
        for (int k = 0; k < 3; k++) {
            v3_t *v = &(clipped_triangles[j].abc[k].xyz);

            // Project this vertex into screen space: draw a ray from the
            // vertex to the camera, intersecting with a plane (the screen).
            v3_t camera_pos_to_v;
            v3_sub(v, camera_pos, &camera_pos_to_v);
            v3_t poi;
            if (v3_ray_plane(&screen_center, camera_fwd, camera_pos, &camera_pos_to_v, &poi)) {
                // Find the intersection, relative to the plane center.
                v3_t poi_rel;
                v3_sub(&poi, &screen_center, &poi_rel);

                // Put this vertex into screenspace.
                double half_screen_w = FAKESCREEN_W / 2,
                       half_screen_h = FAKESCREEN_H / 2;
                screen_triangle.v[k].x = half_screen_w - (v3_dot(&poi_rel, camera_left) * half_screen_w);
                screen_triangle.v[k].y = half_screen_h + (v3_dot(&poi_rel, camera_up) * half_screen_h);
                // If the vertex is behind the camera, then the distance should
                // be negative.  i dont think this is right TODO maybe remove
                // this assert
                double z = v3_dot(&camera_pos_to_v, camera_fwd);
                ASSERT(z >= 0);
                screen_triangle.v[k].z = z;
                memcpy(&screen_triangle.v[k].object, v, sizeof(v3_t));
            } else {
                // If the ray doesn't project onto the screen, it's because the
                // ray is parallel to the screen, so it will be perpendicular
                // to the screen normal.
                ASSERT(0);
            }
        }
        render_draw_screen_triangle(&screen_triangle);
    }
}

void render_draw_screen_triangle(screen_triangle_t *screen_triangle) {
    render_frame_stats.num_triangles++;

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
        render_frame_stats.num_degenerate_triangles++;
        return;
    }

    // If there's only a hi or lo vertex, then there is only one span, e.g.:
    //   *******
    //    *   *
    //     * *
    //      *
    if ((hi == NULL) != (lo == NULL)) {
        if (hi != NULL) {
            draw_span(a, b, c, hi, NULL, screen_triangle);
        } else if (lo != NULL) {
            draw_span(a, b, c, NULL, lo, screen_triangle);
        } else ASSERT(0);
        render_frame_stats.num_onespan_triangles++;
    }

    // If there is both a hi and lo vertex, then we need to draw two spans:
    //       *
    //     * *
    //   *   *
    // *-----*
    //    *  *
    //       *
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
        draw_span(hi, mid, &split, hi, NULL, screen_triangle);
        draw_span(lo, mid, &split, NULL, lo, screen_triangle);
        render_frame_stats.num_twospan_triangles++;
    }
}

void render_end_frame(void) {
    printf(
        "frame %i took %f seconds:\n"
        "    (%i triangles skipped: %i by backface culling, %i by near plane)\n"
        "    %i triangles onscreen (%i onespans, %i twospans, %i degenerate)\n"
        "    %i spans (%i pointups, %i pointdowns, %i degenerates)\n"
        "    (%e pixels rejected by z)\n"
        "    %e pixels drawn\n",
        render_stats.frames_drawn, time_now() - render_frame_stats.start_time,
        render_frame_stats.n_skipped_triangles, render_frame_stats.n_skipped_triangles_by_backface_cull, render_frame_stats.n_skipped_triangles_by_near_plane,
        render_frame_stats.num_triangles, render_frame_stats.num_onespan_triangles, render_frame_stats.num_twospan_triangles, render_frame_stats.num_degenerate_triangles,
        render_frame_stats.num_spans, render_frame_stats.num_pointup_spans, render_frame_stats.num_pointdown_spans, render_frame_stats.num_degenerate_spans,
        (double)render_frame_stats.pixels_rejected_by_z,
        (double)render_frame_stats.pixels_drawn
    );
    render_stats.frames_drawn++;
}

void render_print_stats(void) {
    double elapsed = time_now() - render_stats.start_time;
    printf("%i frames, %e pixels in %f seconds: %f fps, %e pixels/s\n",
        render_stats.frames_drawn,
        (double)render_stats.pixels_drawn,
        elapsed,
        render_stats.frames_drawn / elapsed,
        render_stats.pixels_drawn / elapsed
    );
}
