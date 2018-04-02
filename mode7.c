/*
Software "Mode-7" emulation.
Based almost entirely on this book:
http://www.coranac.com/tonc/text/mode7ex.htm
*/
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "bmp.h"
#include "mode7.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define FOG_RATIO 6/16

/* ===========================================================================
Vectors an matrices
 ===========================================================================*/

typedef struct mat3 {
    /* Column vectors */
    Vector3 u, v, w;
} Mat3;

static Vector3 v3(double x, double y, double z) {
    Vector3 v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

static Vector3 vsub3(Vector3 v1, Vector3 v2) {
    return v3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

static double vdot3(Vector3 v1, Vector3 v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

/*
static Vector3 mat3_vec3(Mat3 m, Vector3 v) {
    Vector3 p;
    p.x = v.x * m.u.x + v.y * m.v.x + v.z * m.w.x;
    p.y = v.x * m.u.y + v.y * m.v.y + v.z * m.w.y;
    p.z = v.x * m.u.z + v.y * m.v.z + v.z * m.w.z;
    return p;
}

static Vector3 mat3_vec3_T(Mat3 m, Vector3 v) {
    Vector3 p;
    p.x = v.x * m.u.x + v.y * m.u.y + v.z * m.u.z;
    p.y = v.x * m.v.x + v.y * m.v.y + v.z * m.v.z;
    p.z = v.x * m.w.x + v.y * m.w.y + v.z * m.w.z;
    return p;
}
*/

static Mat3 make_camera(double phi, double theta) {
    Mat3 m;
    m.u.x = cos(phi); m.v.x = sin(phi)*sin(theta);  m.w.x = -sin(phi)*cos(theta);
    m.u.y = 0;        m.v.y = cos(theta);           m.w.y = sin(theta);
    m.u.z = sin(phi); m.v.z = -cos(phi)*sin(theta); m.w.z = cos(phi)*cos(theta);
    return m;
}

/* ===========================================================================
The actual Mode 7 code
 ===========================================================================*/

static struct {

    /* Position and Dimensions of the view plane */
    int X, Y, W, H;

    /* Distance between the eye and the view plane */
    double D;

    /* Left, right, top and bottom of the frustum */
    double L, R, T, B;
    /* Near/far of the frustum */
    double N, F;

    /* a_cw - Vector placing the camera in the world */
    Vector3 acw;

    /* Orientation of the camera */
    double phi;
    double theta;

    /* Camera matrix, computed from phi and theta */
    Mat3 C;

    int fog_enabled;
    unsigned int fog_color;

    /* Z-buffer */
    double *zbuf;

} mode7;

void m7_init(int x, int y, int w, int h) {

    mode7.X = x;
    mode7.Y = y;
    mode7.W = w;
    mode7.H = h;

    mode7.D = MODE7_D;
    /* Gotcha: L is negative and T is positive */
    mode7.L = -w/2;
    mode7.R = +w/2;
    mode7.T = h/2;
    mode7.B = -h/2;
    mode7.N = MODE7_N;
    mode7.F = MODE7_F;
    mode7.zbuf = malloc(w * h * sizeof *mode7.zbuf);

    m7_set_camera_pos(0, 0, 0);
    m7_set_camera_ang(0, 0);

    mode7.fog_enabled = 0;
    mode7.fog_color = 0;
}

void m7_deinit() {
    free(mode7.zbuf);
}

void m7_enable_fog(unsigned int color) {
    mode7.fog_enabled = 1;
    mode7.fog_color = color;
}
void m7_disable_fog() {
    mode7.fog_enabled = 0;
}

#define ZBUF_SET(x,y,v) mode7.zbuf[(y) * mode7.W + (x)] = v
#define ZBUF_GET(x,y) mode7.zbuf[(y) * mode7.W + (x)]

void m7_clear_zbuf() {
    int i;
    for(i=0;i < mode7.W * mode7.H; i++)
        mode7.zbuf[i] = 1000.0;
}

void m7_set_camera_pos(double x, double y, double z) {
    mode7.acw.x = x;
    mode7.acw.y = y;
    mode7.acw.z = z;
}

void m7_get_camera_pos(double *x, double *y, double *z) {
    *x = mode7.acw.x;
    *y = mode7.acw.y;
    *z = mode7.acw.z;
}

void m7_set_camera_ang(double phi, double theta) {
    while(phi < 0) phi += 2.0 * M_PI;
    while(phi >= 2 * M_PI) phi -= 2.0 * M_PI;
    mode7.phi = phi;

    while(theta < 0) theta += 2.0 * M_PI;
    while(theta >= 2 * M_PI) theta -= 2.0 * M_PI;
    mode7.theta = theta;

    mode7.C = make_camera(mode7.phi, mode7.theta);
}

void m7_get_camera_ang(double *phi, double *theta) {
    *phi = mode7.phi;
    *theta = mode7.theta;
}

void m7_draw_floor(Bitmap *dst, m7_plotfun plotfun, void *data) {
    int py,px;
    double l;

    /* Equation 21.9c: we don't need to draw anything above the horizon */
    int yph = mode7.D / mode7.F * (mode7.F * mode7.C.w.y - mode7.acw.y) / mode7.C.v.y;
    int ysh = mode7.T - yph;
    if(ysh < 0) ysh = 0;

    /* Refer to Eq 21.14 of the abovementioned Mode 7 book. */
    for(py = ysh; py < mode7.H; py++) {
        if((mode7.C.v.y * (py - mode7.T) + mode7.C.w.y * mode7.D) == 0.0)
            continue;
        l = mode7.acw.y / ((py - mode7.T) * mode7.C.v.y + mode7.D * mode7.C.w.y);
        if(l < 0)
            continue;

        for(px = 0; px < mode7.W; px++) {
            double pwx = mode7.acw.x + l * ((px + mode7.L) * mode7.C.u.x + (mode7.T - py) * mode7.C.v.x - mode7.D * mode7.C.w.x);
            double pwz = mode7.acw.z + l * ((px + mode7.L) * mode7.C.u.z + (mode7.T - py) * mode7.C.v.z - mode7.D * mode7.C.w.z);

            unsigned int c = plotfun(data, pwx, pwz);

            if(mode7.fog_enabled) /* [tonc] section 21.5.2 */
                c = bm_lerp(c, mode7.fog_color, l * FOG_RATIO);

            bm_set(dst, px + mode7.X, py + mode7.Y, c);
            ZBUF_SET(px, py, l);
        }
    }
}

void m7_draw_skybox(Bitmap *dst, Bitmap *src, int height, unsigned int background) {
    /* Equation 21.9c: */
    int yph = mode7.D / mode7.F * (mode7.F * mode7.C.w.y - mode7.acw.y) / mode7.C.v.y;
    int ysh = mode7.T - yph + 1;

    if(height <= 0) height = src->h;

    if(ysh >= 0) {
        int x;

        /* Field of view of our camera */
        double fov = atan2(mode7.W/2.0, mode7.D);

        double col = mode7.phi / (2 * M_PI);
        double dsx = fov/mode7.W;

        for(x = 0; x < mode7.W; x++) {
            int sx = (int)(col * src->w) % src->w;

            int dsrcy, ddsty, err, e2, c;

            /* A bit of bresenham magic to scale along the Y */
            int srcy1 = src->h, srcy0 = 0;
            int dsty1 = ysh, dsty0 = ysh - height, y;

            if(dsty0 > mode7.H)
                dsty0 = mode7.H;

            for(y = 0; y < dsty0; y++)
                bm_set(dst, x + mode7.X, y + mode7.Y, background);
            if(dsty0 >= mode7.H)
                continue;

            if(dsty0 < 0) {
                srcy0 = (src->h - 1) * (-dsty0) / height;
                assert(srcy0 >= 0 && srcy0 < src->h);
                dsty0 = 0;
            }
            if(dsty0 >= mode7.H) return;

            dsrcy = srcy1 - srcy0;
            ddsty = dsty1 - dsty0;

            err = dsrcy - ddsty;
            for(;;) {
                if(srcy0 >= src->h || dsty0 >= mode7.H) break;
                c = bm_get(src, sx, srcy0);
                bm_set(dst, x + mode7.X, dsty0 + mode7.Y, c);

                if(srcy0 == srcy1 && dsty0 == dsty1) break;

                e2 = 2 * err;

                if(e2 > -ddsty) {
                    err -= ddsty;
                    srcy0 += 1;
                }
                if(e2 < dsrcy) {
                    err += dsrcy;
                    dsty0 += 1;
                }
            }

            col += dsx;
        }
    }
}

#if 0
#define BLIT_CALLBACK
static int sprite_blit_fun(Bitmap *dst, int dx, int dy, Bitmap *src, int sx, int sy, int mask, void *data) {
    double l = *(double*)data, z;
    int zx = dx - mode7.X;
    int zy = dy - mode7.Y;

    if(zx < 0 || zx > mode7.W || zy < 0 || zy >= mode7.H) return 1;
    z = ZBUF_GET(zx, zy);
    int c = bm_get(src, sx, sy);
    if((c & 0xFFFFFF) != mask && z > l) {
        bm_set(dst, dx, dy, c);
        ZBUF_SET(zx, zy, l);
    }
    return 1;
}
#endif

/* Blits a scaled sprite, taking the distance l and the Z-buffer into account */
static void blit_zbuf(Bitmap *dst, int dx, int dy, int dw, int dh, double l, Bitmap *src, int sx, int sy, int sw, int sh) {
#ifdef BLIT_CALLBACK
    bm_blit_ex_fun(dst, dx, dy, dw, dh, src, sx, sy, sw, sh, sprite_blit_fun, &l);
#else
    /* This particular version might be slightly faster because of how
    it deals with clipping, and there's no function call overhead for the callback */
    int x, y, ssx, sdx;
    int ynum = 0;
    int xnum = 0;
    int mask, c;

    if(sw <= 0 || sh <= 0 || dw <= 0 || dh <= 0)
        return;

    /* Clip on the Y */
    y = dy;
    while(y + mode7.Y < dst->clip.y0 || y < 0 || sy < 0) {
        ynum += sh;
        while(ynum > dh) {
            ynum -= dh;
            sy++;
        }
        y++;
    }
    if(dy + mode7.Y >= dst->clip.y1 || dy > mode7.H || dy + dh + mode7.Y < dst->clip.y0 || sy >= src->h)
        return;

    /* Clip on the X */
    x = dx;
    while(x + mode7.X < dst->clip.x0 || x < 0 || sx < 0) {
        xnum += sw;
        while(xnum > dw) {
            xnum -= dw;
            sx++;
        }
        x++;
    }

    if(dx + mode7.X >= dst->clip.x1 || dx > mode7.W || dx + dw + mode7.X < dst->clip.x0 || sx >= src->w)
        return;

    mask = bm_get_color(src) & 0xFFFFFF;

    ssx = sx; /* Save sx for the next row */
    sdx = x;
    for(; y < dy + dh; y++){
        if(sy >= src->h || y + mode7.Y >= dst->clip.y1 || y >= mode7.H)
            break;
        xnum = 0;
        sx = ssx;

        assert(y + mode7.Y >= dst->clip.y0 && sy >= 0);
        for(x = sdx; x < dx + dw; x++) {
            if(sx >= src->w || x + mode7.X >= dst->clip.x1 || x >= mode7.W)
                break;

            if(ZBUF_GET(x, y) > l) {
                assert(x + mode7.X >= dst->clip.x0 && sx >= 0);

                c = bm_get(src, sx, sy);
                if((c & 0xFFFFFF) != mask) {
                    if(mode7.fog_enabled) /* [tonc] section 21.5.2 */
                        c = bm_lerp(c, mode7.fog_color, l * FOG_RATIO);
                    bm_set(dst, x + mode7.X, y + mode7.Y, c);
                    ZBUF_SET(x, y, l);
                }
            }

            xnum += sw;
            while(xnum > dw) {
                xnum -= dw;
                sx++;
            }
        }
        ynum += sh;
        while(ynum > dh) {
            ynum -= dh;
            sy++;
        }
    }
#endif
}

void m7_draw_sprite(Bitmap *dst, double wx, double wy, double wz, Bitmap *src, int sx, int sy, int sw, int sh) {

    Vector3 xw, asp, r;
    double zc, l;

    Vector3 xp, xs;

    int dw, dh;

    xw = v3(wx, wy, wz);
    asp = v3(mode7.L, mode7.T, -mode7.D);

    r = vsub3(xw, mode7.acw);

    zc = vdot3(mode7.C.w, r);
    l = -zc/mode7.D;

    xp.x = vdot3(mode7.C.u, r)/l;
    xp.y = vdot3(mode7.C.v, r)/l;
    xp.z = -mode7.D;

    xs = vsub3(xp, asp);
    xs.y = -xs.y;

    if(-zc < mode7.N || -zc > mode7.F) return; /* near/far clipping plane */

    dw = (double)sw/l;
    dh = (double)sh/l;

#if MODE7_ANCHOR_BOTTOM
    /* To anchor the sprites at the bottom: */
    blit_zbuf(dst, xs.x - dw/2, xs.y - dh, dw, dh, l, src, sx, sy, sw, sh);
#else
    blit_zbuf(dst, xs.x - dw/2, xs.y - dh/2, dw, dh, l, src, sx, sy, sw, sh);
#endif
}

static void horiz_line(Bitmap *bmp, double sx, double ex, int y, double sz, double ez, int col) {
    int x;
    double t;
    if(sx > ex) {
        t = sx; sx = ex; ex = t;
        t = sz; sz = ez; ez = t;
    }

    double mz = (ez - sz) / (ex - sx);
    double z = sz;
    for(x = sx; x <= ex; x++, z += mz) {
        if(x < 0) continue;
        else if(x >= mode7.W) break;
        if(ZBUF_GET(x, y) > z) {

            unsigned int c = col;
            if(mode7.fog_enabled) /* [tonc] section 21.5.2 */
                c = bm_lerp(c, mode7.fog_color, z * FOG_RATIO);

            bm_set(bmp, x + mode7.X, y + mode7.Y, c);
            ZBUF_SET(x, y, z);
        }
    }
}

static void fill_tri(Bitmap *bmp, Vector3 p0, Vector3 p1, Vector3 p2) {
    /* http://www-users.mat.uni.torun.pl/~wrona/3d_tutor/tri_fillers.html */

    /* This rounding is necessary to remove some artefacts */
    p0.x = (int)p0.x; p0.y = (int)p0.y;
    p1.x = (int)p1.x; p1.y = (int)p1.y;
    p2.x = (int)p2.x; p2.y = (int)p2.y;

    if(p1.y < p0.y) {
        Vector3 t = p0;
        p0 = p1;
        p1 = t;
    }
    if(p2.y < p0.y) {
        Vector3 t = p0;
        p0 = p2;
        p2 = t;
    }
    if(p2.y < p1.y) {
        Vector3 t = p1;
        p1 = p2;
        p2 = t;
    }

    double Sx = p0.x, Sz = p0.z;
    double Ex = p0.x, Ez = p0.z;
    double dx1, dx2, dx3;
    double dz1, dz2, dz3;

    unsigned int color = bm_get_color(bmp);

    if((p1.y - p0.y) > 0) {
        dx1 = (p1.x - p0.x) / (p1.y - p0.y);
        dz1 = (p1.z - p0.z) / (p1.y - p0.y);
    } else {
        dx1 = 0;
        dz1 = 0;
    }
    if((p2.y - p0.y) > 0) {
        dx2 = (p2.x - p0.x) / (p2.y - p0.y);
        dz2 = (p2.z - p0.z) / (p2.y - p0.y);
    } else {
        dx2 = 0;
        dz2 = 0;
    }
    if((p2.y - p1.y) > 0) {
        dx3 = (p2.x - p1.x) / (p2.y - p1.y);
        dz3 = (p2.z - p1.z) / (p2.y - p1.y);
    } else {
        dx3 = 0;
        dz3 = 0;
    }

    int y = p0.y;
    if(dx1 > dx2) {
        for(;y < p1.y; y++, Sx += dx2, Ex += dx1, Sz += dz2, Ez += dz1) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, Sx, Ex, y, Sz, Ez, color);
        }
        Ex = p1.x; Ez = p1.z;
        for(;y <= p2.y; y++, Sx += dx2, Ex += dx3, Sz += dz2, Ez += dz3) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, Sx, Ex, y, Sz, Ez, color);
        }
    } else {
        for(;y < p1.y; y++, Sx += dx1, Ex += dx2, Sz += dz1, Ez += dz2) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, Sx, Ex, y, Sz, Ez, color);
        }
        Sx = p1.x; Sz = p1.z;
        for(;y <= p2.y; y++, Sx += dx3, Ex += dx2, Sz += dz3, Ez += dz2) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, Sx, Ex, y, Sz, Ez, color);
        }
    }
}

void m7_draw_tri(Bitmap *bmp, Vector3 tri[3]) {
    Vector3 proj[3];
    int i;

    Vector3 asp = v3(mode7.L, mode7.T, -mode7.D);

    for(i = 0; i < 3; i++) {
        Vector3 xp;
        Vector3 r = vsub3(tri[i], mode7.acw);

        double zc = vdot3(mode7.C.w,r);
        double l = -zc/mode7.D;
        double x_p = vdot3(mode7.C.u,r)/l;
        double y_p = vdot3(mode7.C.v,r)/l;

        if(-zc < mode7.N || -zc > mode7.F) return; /* near/far clipping plane */

        xp = v3(x_p, y_p, -mode7.D);
        proj[i] = vsub3(xp, asp);
        proj[i].y = -proj[i].y;
        proj[i].z = l;
    }
    fill_tri(bmp, proj[0], proj[1], proj[2]);
}

static void draw_line(Bitmap *b, int x0, int y0, double z0, int x1, int y1, double z1) {
    int dx = x1 - x0;
    int dy = y1 - y0;
    int sx, sy;
    int err, e2;

    if(dx < 0) dx = -dx;
    if(dy < 0) dy = -dy;

    if(x0 < x1)
        sx = 1;
    else
        sx = -1;
    if(y0 < y1)
        sy = 1;
    else
        sy = -1;

    err = dx - dy;

    int x = x0, y = y0, steps = 0;

    /* Run the loop once, without plotting anything so that we
    can figure out how we'd need to adjust the z value each step.
    I'm sure there's a better way to do this, but I'm stumped. */
    for(;; steps++) {
        if(x == x1 && y == y1) break;
        e2 = 2 * err;
        if(e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if(e2 < dx) {
            err += dx;
            y += sy;
        }
    }
    if(!steps) return;

    double dz = (z1 - z0)/steps;
    unsigned int c = b->color;

    x = x0, y = y0;
    for(;;) {
        if(x >= 0 && x < mode7.W && y >= 0 && y < mode7.H){// && z0 >= 0 && z0 < 1.0) {
            if(ZBUF_GET(x, y) > z0 + DBL_EPSILON) {
                unsigned int c1 = c;
                if(mode7.fog_enabled) /* [tonc] section 21.5.2 */
                    c1 = bm_lerp(c, mode7.fog_color, z0 * FOG_RATIO);
                ZBUF_SET(x, y, z0);
                bm_set(b, x + mode7.X, y + mode7.Y, c1);
            }
        }
        if(x == x1 && y == y1) break;
        e2 = 2 * err;
        if(e2 > -dy) {
            err -= dy;
            x += sx;
        }
        if(e2 < dx) {
            err += dx;
            y += sy;
        }
        z0 += dz;
    }
}

void m7_line(Bitmap *bmp, Vector3 p[2]) {
    Vector3 proj[2];
    int i;
    Vector3 asp = v3(mode7.L, mode7.T, -mode7.D);

    for(i = 0; i < 2; i++) {
        Vector3 xp;
        Vector3 r = vsub3(p[i], mode7.acw);

        double zc = vdot3(mode7.C.w,r);
        double l = -zc/mode7.D;
        double x_p = vdot3(mode7.C.u,r)/l;
        double y_p = vdot3(mode7.C.v,r)/l;

        if(-zc < mode7.N || -zc > mode7.F) return; /* near/far clipping plane */

        xp = v3(x_p, y_p, -mode7.D);
        proj[i] = vsub3(xp, asp);
        proj[i].y = -proj[i].y;
        proj[i].z = l;
    }
    draw_line(bmp, proj[0].x, proj[0].y, proj[0].z, proj[1].x, proj[1].y, proj[1].z);
}

void m7_lookat(double x, double y, double z) {
    double dx =  x - mode7.acw.x;
    double dy =  y - mode7.acw.y;
    double dz =  z - mode7.acw.z;
    double phi = atan2(dx, -dz);
    double theta = -atan2(dy, sqrt(dx * dx + dz * dz));
    m7_set_camera_ang(phi, theta);
}

double m7_rel_angle(double phi_o) {
    double alpha = 0; /* = atan2(xc, -zc); // Ignored */
    double w = phi_o - mode7.phi - alpha;
    while (w < 0.0) w += 2 * M_PI;
    while (w >= 2 * M_PI) w -= 2 * M_PI;
    return w;
}
