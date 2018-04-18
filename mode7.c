/*
Software "Mode-7" emulation.
Based almost entirely on this book:
http://www.coranac.com/tonc/text/mode7ex.htm
*/
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "bmp.h"
#include "obj.h"
#include "mode7.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define FOG_RATIO 6/16

/* These are defined here, because it is unlikely they will require changing */
/* value for D, the draw plane distance */
#define MODE7_D 256
/* value for N, the near clipping pane */
#define MODE7_N +24
/* value for F, the far clipping pane */
#define MODE7_F +1024

/* For performance reasons you might want
to disable the stencil buffer here */
#ifndef M7_STENCIL
#  define M7_STENCIL 1
#endif

/* Do we support the `l` command in OBJ files? */
#ifndef M7_OBJ_LINES
#  define M7_OBJ_LINES 1
#endif

/* If there is z-fighting between a line and another element, let the line win. */
#ifndef LINE_ZBUF_FACTOR
#  define LINE_ZBUF_FACTOR 0.005
#endif

/* ===========================================================================
Vectors an matrices
 ===========================================================================*/

typedef struct mat3 {
    /* Column vectors */
    Vector3 u, v, w;
} Mat3;

Vector3 v3(double x, double y, double z) {
    Vector3 v;
    v.x = x; v.y = y; v.z = z;
    return v;
}

Vector3 v3_add(Vector3 v1, Vector3 v2) {
    return v3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
}

Vector3 v3_sub(Vector3 v1, Vector3 v2) {
    return v3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
}

int v3_equal(Vector3 v1, Vector3 v2) {
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

double v3_dot(Vector3 v1, Vector3 v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Vector3 v3_cross(Vector3 u, Vector3 v) {
    return v3(u.y*v.z - u.z*v.y, u.z*v.x - u.x*v.z, u.x*v.y - u.y*v.x);
}

Vector3 v3_normalize(Vector3 v) {
    double _len = 1.0/sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    return v3(v.x*_len, v.y*_len, v.z*_len);
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

typedef float ZBUF_TYPE;

#define ZBUF_SET(x,y,v) mode7.zbuf[(y) * mode7.W + (x)] = (v)
#define ZBUF_GET(x,y) mode7.zbuf[(y) * mode7.W + (x)]

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

    Vector3 asp;

    /* Orientation of the camera */
    double phi;
    double theta;

    /* Camera matrix, computed from phi and theta */
    Mat3 C;

    int fog_enabled;
    unsigned int fog_color;

    int backface_enabled;

    enum M7_ANCHOR anchor;

    /* Z-buffer */
    ZBUF_TYPE *zbuf;

    int stencil_enabled;
    Bitmap *stencil;

} mode7;

void m7_init(int x, int y, int w, int h) {

    memset(&mode7, 0, sizeof mode7);

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

    mode7.backface_enabled = 1;

    mode7.anchor = M7_ANCHOR_CENTER;

    mode7.asp = v3(mode7.L, mode7.T, -mode7.D);

    mode7.stencil_enabled = 1;
    mode7.stencil = bm_create(w,h);
}

void m7_dims(int *X, int *Y, int *W, int *H) {
    if(X) *X = mode7.X;
    if(Y) *Y = mode7.Y;
    if(W) *W = mode7.W;
    if(H) *H = mode7.H;
}

void m7_deinit() {
    free(mode7.zbuf);
    bm_free(mode7.stencil);
}

void m7_enable_fog(unsigned int color) {
    mode7.fog_enabled = 1;
    mode7.fog_color = color;
}
void m7_disable_fog() {
    mode7.fog_enabled = 0;
}
void m7_backface(int enable) {
    mode7.backface_enabled = enable;
}
void m7_anchor_mode(enum M7_ANCHOR mode) {
    mode7.anchor = mode;
}

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
            if(mode7.stencil_enabled)
                bm_putpixel(mode7.stencil, px, py);
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

/* Blits a scaled sprite, taking the distance l and the Z-buffer into account */
static void blit_zbuf(Bitmap *dst, int dx, int dy, int dw, int dh, double l, Bitmap *src, int sx, int sy, int sw, int sh) {

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
                    if(mode7.stencil_enabled)
                        bm_putpixel(mode7.stencil, x, y);
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
}

void m7_draw_sprite(Bitmap *dst, double wx, double wy, double wz, Bitmap *src, int sx, int sy, int sw, int sh) {

    Vector3 xw, r;
    double zc, l;

    Vector3 xp, xs;

    int dw, dh;

    xw = v3(wx, wy, wz);

    r = v3_sub(xw, mode7.acw);

    zc = v3_dot(mode7.C.w, r);
    l = -zc/mode7.D;

    xp.x = v3_dot(mode7.C.u, r)/l;
    xp.y = v3_dot(mode7.C.v, r)/l;
    xp.z = -mode7.D;

    xs = v3_sub(xp, mode7.asp);
    xs.y = -xs.y;

    if(-zc < mode7.N || -zc > mode7.F) return; /* near/far clipping plane */

    dw = (double)sw/l;
    dh = (double)sh/l;

    if(mode7.anchor == M7_ANCHOR_BOTTOM) {
        /* To anchor the sprites at the bottom: */
        blit_zbuf(dst, xs.x - dw/2, xs.y - dh, dw, dh, l, src, sx, sy, sw, sh);
    } else {
        blit_zbuf(dst, xs.x - dw/2, xs.y - dh/2, dw, dh, l, src, sx, sy, sw, sh);
    }
}

static void horiz_line(Bitmap *bmp, int sx, int ex, int y, double sz, double ez, int col) {
    int x;
    double t;
    if(sx > ex) {
        t = sx; sx = ex; ex = t;
        t = sz; sz = ez; ez = t;
    }

    double mz = ex > sx ? (ez - sz) / (ex - sx) : 0;

    assert(ex >= sx);
    for(x = sx; x <= ex; x++, sz += mz) {
        if(x < 0) continue;
        else if(x >= mode7.W) break;
        if(ZBUF_GET(x, y) > sz) {

            unsigned int c = col;
            if(mode7.fog_enabled)
                c = bm_lerp(c, mode7.fog_color, sz * FOG_RATIO);

            bm_set(bmp, x + mode7.X, y + mode7.Y, c);
            ZBUF_SET(x, y, sz);
            if(mode7.stencil_enabled)
                bm_putpixel(mode7.stencil, x, y);
        }
    }
}

static void rasterize_interp(Bitmap *bmp, const Vector3 v0, const Vector3 v1, const Vector3 v2) {

    int v0y = v0.y, v1y = v1.y, v2y = v2.y;

    int y = v0y;

    if(v2y == v0y) return; /* degenerate case. */

    double mxB = (v2.x - v0.x)/(v2.y - v0.y);
    double mzB = (v2.z - v0.z)/(v2.y - v0.y);

    double xA = v0.x, xB = v0.x;
    double zA = v0.z, zB = v0.z;

    if(v1y > v0y) {
        double mxA = (v1.x - v0.x)/(v1.y - v0.y);
        double mzA = (v1.z - v0.z)/(v1.y - v0.y);
        for(; y < v1y; y++, xA += mxA, zA += mzA, xB += mxB, zB += mzB) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, xA, xB, y, zA, zB, bmp->color);
        }
    }
    xA = v1.x;
    zA = v1.z;

    if(v2y > v1y) {
        double mxC = (v2.x - v1.x)/(v2.y - v1.y);
        double mzC = (v2.z - v1.z)/(v2.y - v1.y);
        for(; y < v2y; y++, xA += mxC, zA += mzC, xB += mxB, zB += mzB) {
            if(y < 0) continue;
            else if(y >= mode7.H) return;
            horiz_line(bmp, xA, xB, y, zA, zB, bmp->color);
        }
    }
}

static void rasterize(Bitmap *bmp, const Vector3 v[3]) {

    if(v[0].y < v[1].y) {
        if(v[0].y < v[2].y) {
            if(v[1].y < v[2].y)
                rasterize_interp(bmp, v[0], v[1], v[2]);
            else
                rasterize_interp(bmp, v[0], v[2], v[1]);
        } else {
            assert(v[2].y <= v[1].y);
            rasterize_interp(bmp, v[2], v[0], v[1]);
        }
    } else {
        if(v[1].y < v[2].y) {
            if(v[0].y < v[2].y)
                rasterize_interp(bmp, v[1], v[0], v[2]);
            else
                rasterize_interp(bmp, v[1], v[2], v[0]);
        } else {
            assert(v[2].y <= v[0].y);
            rasterize_interp(bmp, v[2], v[1], v[0]);
        }
    }
}

void m7_draw_tri(Bitmap *bmp, Vector3 tri[3]) {
    Vector3 proj[3];
    int i;

    for(i = 0; i < 3; i++) {
        Vector3 xp;
        Vector3 r = v3_sub(tri[i], mode7.acw);

        double zc = v3_dot(mode7.C.w,r);
        double l = -zc/mode7.D;
        double x_p = v3_dot(mode7.C.u,r)/l;
        double y_p = v3_dot(mode7.C.v,r)/l;

        if(-zc < mode7.N || -zc > mode7.F) return; /* near/far clipping plane */

        xp = v3(x_p, y_p, -mode7.D);
        proj[i] = v3_sub(xp, mode7.asp);
        proj[i].y = -proj[i].y;
        proj[i].z = l;
    }

    if(!mode7.backface_enabled) {
        /* Backface culling */
        Vector3 n = (v3_cross(v3_sub(proj[2],proj[0]), v3_sub(proj[1],proj[0])));
        Vector3 view_dir = v3(0,0,1);
        if(v3_dot(n, view_dir) < 0)
            return;
    }
    rasterize(bmp, proj);
}

void m7_draw_obj(Bitmap *dst, ObjMesh *obj, Vector3 pos, double yrot, unsigned int color) {
    int i, j;
    if(!obj)
        return;

    Vector3 light_dir = v3_normalize(v3(3.0, -5.0, -1.0));

    double sinAngle = sin(yrot);
    double cosAngle = cos(yrot);

    for(i = 0; i < obj->nfaces; i++) {
        ObjFace *f = obj->faces + i;
        Vector3 tri[3];

        for(j = 0; j < 3; j++) {
            assert(f->verts[j].v < obj->nverts);
            ObjVert *v = obj->verts + f->verts[j].v;
            tri[j].x = -cosAngle * v->x + sinAngle * v->z + pos.x;
            tri[j].y = v->y + pos.y;
            tri[j].z = -sinAngle * v->x - cosAngle * v->z + pos.z;
        }
        Vector3 n = v3_normalize(v3_cross(v3_sub(tri[2],tri[0]), v3_sub(tri[1],tri[0])));

        double intensity = (v3_dot(n, light_dir) + 1) / 2.0;

        intensity = intensity * 0.2 + 0.8;
        unsigned int light = bm_lerp(0, color, intensity);
        bm_set_color(dst, light);

        m7_draw_tri(dst, tri);
    }

#if M7_OBJ_LINES
    bm_set_color(dst, 0);
    for(i = 0; i < obj->nlines; i++) {
        Vector3 p[2];
        int j;
        ObjLine *l = &obj->lines[i];
        ObjVert *v = obj->verts + l->idx[0];
        p[1].x = (-cosAngle * v->x + sinAngle * v->z) + pos.x;
        p[1].y = (v->y) + pos.y;
        p[1].z = (-sinAngle * v->x - cosAngle * v->z) + pos.z;
        for(j = 1; j < l->n; j++) {
            p[0] = p[1];
            v = obj->verts + l->idx[j];
            p[1].x = (-cosAngle * v->x + sinAngle * v->z) + pos.x;
            p[1].y = (v->y) + pos.y;
            p[1].z = (-sinAngle * v->x - cosAngle * v->z) + pos.z;
            m7_line(dst, p[0], p[1]);
        }
    }
#endif
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
        if(x >= 0 && x < mode7.W && y >= 0 && y < mode7.H){
            if(ZBUF_GET(x, y) > z0 - LINE_ZBUF_FACTOR) {
                unsigned int c1 = c;
                if(mode7.fog_enabled) /* [tonc] section 21.5.2 */
                    c1 = bm_lerp(c, mode7.fog_color, z0 * FOG_RATIO);
                ZBUF_SET(x, y, z0);
                bm_set(b, x + mode7.X, y + mode7.Y, c1);
                if(mode7.stencil_enabled)
                    bm_putpixel(mode7.stencil, x, y);
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

int m7_project(Vector3 p, Vector3 *o) {
    Vector3 xp;
    Vector3 r = v3_sub(p, mode7.acw);

    double zc = v3_dot(mode7.C.w,r);
    double l = -zc/mode7.D;
    double x_p = v3_dot(mode7.C.u,r)/l;
    double y_p = v3_dot(mode7.C.v,r)/l;

    if(-zc < mode7.N || -zc > mode7.F) return 0; /* near/far clipping plane */

    xp = v3(x_p, y_p, -mode7.D);
    *o = v3_sub(xp, mode7.asp);
    o->y = - o->y;
    o->z = l;
    return 1;
}

void m7_line(Bitmap *bmp, Vector3 p0, Vector3 p1) {
    Vector3 proj[2];

    if(!m7_project(p0, &proj[0]) ||
        !m7_project(p1, &proj[1]))
        return;

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

/* Enables/disables the stencil buffer */
void m7_stencil_enable(int enable) {
    mode7.stencil_enabled = enable;
}

void m7_set_stencil(unsigned int color) {
    bm_set_color(mode7.stencil, color);
}
void m7_clear_stencil() {
    bm_set_color(mode7.stencil, 0);
    bm_clear(mode7.stencil);
}
Bitmap *m7_get_stencil() {
    return mode7.stencil;
}
unsigned int m7_stencil_enable_at(int x, int y) {
    assert(x >= 0 && x < mode7.stencil->w);
    assert(y >= 0 && y < mode7.stencil->h);
    return bm_get(mode7.stencil, x, y);
}
