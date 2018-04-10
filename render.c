#include <stdlib.h>
#include <math.h>
#include <assert.h>

#if defined(SDL2) || defined(SDL) || defined(__EMSCRIPTEN__)
#  include "pocadv.h"
#else
#  include "gdi.h"
#endif

#include "mode7.h"


static Bitmap *tiles = NULL;
static Bitmap *props = NULL;
static Bitmap *sky = NULL;
static Bitmap *guy = NULL;

/* Some variables for the sprite. fix later */
static double spr_frame = 0;
static double spr_angle = 0;
static double spr_x, spr_z;
static int is_animating = 0;

static double fwd_speed = 100.0;
static double rot_speed = 180.0;

ObjMesh *obj = NULL;
Vector3 obj_pos = {160, 0, 160};
double obj_rot = 0.0;

Vector3 prop_positions[] = {
    {0, 0, 320},
    {0, 0, 0},
    {320, 0, 0},
    {320, 0, 320},
};
int npoints = (sizeof prop_positions)/sizeof(prop_positions[0]);

void draw_props() {
    int i;
    for(i = 0; i < npoints; i++) {
        m7_set_stencil(0xFF0000 + (i << 8));
        m7_draw_sprite(screen, prop_positions[i].x, prop_positions[i].y, prop_positions[i].z, props, (i & 0x01)?40:0, 0, 40, 60);
    }
}

unsigned int the_plotfun(void *data, double pwx, double pwz) {
    int c = 0x8090A0;

    if(pwx >= 0 && pwx < 320 && pwz >= 0 && pwz < 320) {

        int mx = pwx / 20;
        int my = pwz / 20;
        // TODO: Use mx/my to look up the tile on a "map" data structure
        int tile_col = 4;
        int tile_row = 1;
        int tx = tile_col * 20;
        int ty = tile_row * 20;
        (void)mx;(void)my;
        (void)tx;(void)ty;

        int sx = (int)pwx % 20;
        int sy = (int)pwz % 20;

        c = bm_picker(tiles, tx + sx, ty + sy);
    }
    return c;
}

void update_all(double elapsed) {

    double x, y, z;
    double phi, theta;

    m7_get_camera_pos(&x, &y, &z);
    m7_get_camera_ang(&phi, &theta);

    //obj_rot += elapsed * (rot_speed/4.0) * M_PI / 180;

    /* Press up and down to move foraward/backward */
    if(keys[KCODE(UP)]) {
        z -= elapsed * fwd_speed * cos(phi);
        x += elapsed * fwd_speed * sin(phi);
    }
    if(keys[KCODE(DOWN)]) {
        z += elapsed * fwd_speed * cos(phi);
        x -= elapsed * fwd_speed * sin(phi);
    }

    /* Left and right to turn the camera */
    if(keys[KCODE(RIGHT)])
        phi += elapsed * rot_speed * M_PI / 180;
    if(keys[KCODE(LEFT)])
        phi -= elapsed * rot_speed * M_PI / 180;

    /* PgUp and PgDn to look up and down */
    if(keys[KCODE(PAGEUP)])
        theta -= elapsed * rot_speed * M_PI / 180;
    if(keys[KCODE(PAGEDOWN)])
        theta += elapsed * rot_speed * M_PI / 180;

    /* Press 'q' and 'z' to raise and lower the camera */
    if(keys[KCODEA(q,Q)])
        y += elapsed * 100;
    if(keys[KCODEA(z,Z)] && y > 0)
        y -= elapsed * 100;

    is_animating = 0;
    if(keys[KCODEA(h,H)]) {
        is_animating = 1;
        spr_angle = 270 * M_PI / 180;
        spr_x -= elapsed * fwd_speed;
    }
    if(keys[KCODEA(j,J)]) {
        is_animating = 1;
        spr_angle = 0 * M_PI / 180;
        spr_z -= elapsed * fwd_speed;
    }
    if(keys[KCODEA(k,K)]) {
        is_animating = 1;
        spr_angle = 180 * M_PI / 180;
        spr_z += elapsed * fwd_speed;
    }
    if(keys[KCODEA(l,L)]) {
        is_animating = 1;
        spr_angle = 90 * M_PI / 180;
        spr_x += elapsed * fwd_speed;
    }

    /* Set the position and angle of the camera */
    m7_set_camera_pos(x, y, z);
    m7_set_camera_ang(phi, theta);

    /* Press 'l' to look at the sprite - only call `m7_lookat`
    after you've called `m7_set_camera_pos` */
    if(keys[KCODEA(a,A)])
        m7_lookat(spr_x, 16, spr_z);

    /* Press 'b' to dump the stencil buffer to a file */
    if(keys[KCODEA(b,B)]) {
        Bitmap *stencil = m7_get_stencil();
        bm_save(stencil, "stencil.bmp");
    }
}

int render(double elapsed) {
    double x, y, z;
    double phi, theta;

    int row = 1, col = 1;

    double rel_angle = m7_rel_angle(spr_angle);

    update_all(elapsed);

    bm_set_color(screen, bm_atoi("black"));
    bm_clear(screen);

    bm_set_color(screen, bm_atoi("white"));
    bm_rect(screen, 9, 9, SCREEN_WIDTH - 10, SCREEN_HEIGHT - 10);

    m7_clear_zbuf();

    m7_clear_stencil();

    m7_draw_skybox(screen, sky, sky->h * 4, 0xC0D0E0);

    m7_set_stencil(0x00FF00);
    m7_draw_floor(screen, the_plotfun, NULL);

    draw_props();

    m7_set_stencil(0x0000FF);
    m7_draw_obj(screen, obj, obj_pos, obj_rot, bm_atoi("#51687F"));

    Vector3 p[2];
    bm_set_color(screen, 0);

    p[0] = v3(140, 0, 140);
    p[1] = v3(140, 20, 140);
    m7_line(screen, p);
    p[0] = v3(180, 0, 140);
    p[1] = v3(180, 20, 140);
    m7_line(screen, p);
    p[0] = v3(180, 0, 180);
    p[1] = v3(180, 20, 180);
    m7_line(screen, p);
    p[0] = v3(140, 0, 180);
    p[1] = v3(140, 20, 180);
    m7_line(screen, p);

    if(is_animating) {
        spr_frame += elapsed;
        switch(((int)(spr_frame * 10)) % 4) {
        case 0: col = 0; break;
        case 1: col = 1; break;
        case 2: col = 2; break;
        case 3: col = 1; break;
        }
    }

    rel_angle += M_PI / 4;
    if(rel_angle >= 2.0 * M_PI) rel_angle -= 2.0 * M_PI;
    switch((int)(4.0 * rel_angle / (2.0 * M_PI))) {
        case 0: row = 0; break;
        case 1: row = 1; break;
        case 2: row = 2; break;
        case 3: row = 3; break;
    }

    m7_set_stencil(0xFFFF00);
    m7_draw_sprite(screen, spr_x, 0, spr_z, guy, col * 24, row * 32, 24, 32);

    m7_get_camera_pos(&x, &y, &z);
    m7_get_camera_ang(&phi, &theta);

#if 0
    bm_set_color(screen, bm_atoi("white"));
    bm_printf(screen, 5, screen->h - 8, "rel: %.2f (%.2f): %d => %d", rel_angle*180/M_PI, spr_angle*180/M_PI, (int)(4.0 * rel_angle / (2.0 * M_PI)), row);
    bm_printf(screen, 10, 20, "p:%.f t:%.f", phi*180/M_PI, theta*180/M_PI);
    bm_printf(screen, 10, 30, "x:%.f y:%.f z:%.f", x, y, z);
#endif

    return 1;
}

void init_game(int argc, char *argv[]) {

    spr_x = 160; spr_z = 200;
    spr_angle = 180 * M_PI / 180;

    m7_init(10, 10, SCREEN_WIDTH - 20, SCREEN_HEIGHT - 20);

    m7_enable_fog(0x8090A0);

    m7_set_camera_pos(160, 25, 320);
    m7_set_camera_ang(0, 0 * M_PI / 180);

    /*m7_lookat(spr_x, 12, spr_z);*/

    tiles = bm_load("res/HardVacuum.gif");
    if(!tiles) {
        exit_error("Unable to load tiles");
    }

    sky = bm_load("res/sky.gif");
    if(!sky) {
        exit_error("Unable to load sky");
    }

    guy = bm_load("res/char.gif");
    if(!guy) {
        exit_error("Unable to load guy");
    }
    bm_set_color(guy, bm_atoi("#209C00"));

    props = bm_load("res/buildings.gif");
    if(!props) {
        exit_error("Unable to load props");
    }
    bm_set_color(props, bm_atoi("#008A76"));

    obj = obj_load("res/pyramid.obj");
    if(!obj)
        exit_error("Unable to load OBJ");
}

void deinit_game() {
    bm_free(tiles);
    bm_free(props);
    bm_free(guy);
    bm_free(sky);

    obj_free(obj);

    m7_deinit();
}
