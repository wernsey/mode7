#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <errno.h>
#include <setjmp.h>
#include <assert.h>

#include "obj.h"

#if defined(WIN32) && defined(_MSC_VER)
#  define SAFE_C11
#  define strdup _strdup
#endif

typedef struct {
    int sym,
        lastsym;
    char word[64],
        lastword[64];
    char *text;
    char *lex;

    int s; /* smoothing group */

    enum {MODE_OBJ, MODE_MTL} mode;

    jmp_buf buf; /* Exception handling */
    int line;    /* Error reporting */
    char error_buf[128];

} OBJ_Parser;

enum { V = 128, VT, VN, VP, F, L, G, S, MTLLIB, USEMTL, ID, NUM ,
    NEWMTL, KA, KD, KS, TF, TR, ILLUM, NS, MAP_KD
} Types;

struct {
    int t;
    const char *name;
} OBJ_Command_Names[] = {
    {V,"v"},
    {VT,"vt"},
    {VN,"vn"},
    {VP,"vp"},
    {F,"f"},
    {L,"l"},
    {G,"g"},
    {S,"s"},
    {MTLLIB,"mtllib"},
    {USEMTL,"usemtl"},
    {-1, NULL}
},
MTL_Command_Names[] = {
    {NEWMTL,"newmtl"},
    {KA,"Ka"},
    {KD,"Kd"},
    {KS,"Ks"},
    {TF,"Tf"},
    {TR,"Tr"},
    {ILLUM,"illum"},
    {NS,"Ns"},
    {MAP_KD,"map_Kd"},
    {-1, NULL}
};
static const char *typename(int type) {
    int i;
    if(type == ID)
        return "identifier";
    if(type == NUM)
        return "number";
    if(type == '\n')
        return "'\\n'";
    if(type == '\0')
        return "EOF";
    for(i=0; OBJ_Command_Names[i].name; i++) {
        if(OBJ_Command_Names[i].t == type) return OBJ_Command_Names[i].name;
    }
    for(i=0; MTL_Command_Names[i].name; i++) {
        if(MTL_Command_Names[i].t == type) return MTL_Command_Names[i].name;
    }
    return "--";
}

static void nextsym(OBJ_Parser *p);

static void error(OBJ_Parser *p, const char *msg, ...) {
    va_list ap;
    unsigned int len;
    va_start(ap, msg);
    snprintf(p->error_buf, sizeof p->error_buf, "OBJ: error:%d: ", p->line);
    len = strlen(p->error_buf);
    assert(sizeof p->error_buf - len > 0);
    vsnprintf(p->error_buf + len, sizeof p->error_buf - len, msg, ap);
    va_end(ap);
    longjmp(p->buf, 1);
}

static int accept(OBJ_Parser *p, int s) {
    if(p->sym == s) {
        nextsym(p);
        return 1;
    }
    return 0;
}

static int expect(OBJ_Parser *p, int s) {
    if(accept(p, s))
        return 1;
    error(p,"`%s` expected\n", typename(s));
    return 0;
}

static void _init(OBJ_Parser *p, char *text) {
    p->sym = 0;
    p->word[0] = 0;
    p->lastsym = 0;
    p->lastword[0] = 0;
    p->text = text;
    p->lex = text;
    p->line = 1;
    p->s = 0;
}
static void init_OBJ(OBJ_Parser *p, char *text) {
    _init(p, text);
    p->mode = MODE_OBJ;
    nextsym(p);
}
static void init_MTL(OBJ_Parser *p, char *text) {
    _init(p, text);
    p->mode = MODE_MTL;
    nextsym(p);
}

#define SYMBOL(s) do{p->sym=s;return;}while(0)
static void nextsym(OBJ_Parser *p) {
    int l = 0;
    p->lastsym = p->sym;
    memcpy(p->lastword, p->word, sizeof p->lastword);
    p->word[0] = 0;
space:
    while(isspace(p->lex[0])) {
        if(p->lex[0] == '\n') {
            p->line++;
            p->lex++;
            SYMBOL('\n');
        }
        p->lex++;
    }
    if(!p->lex[0])
        SYMBOL('\0');
    else if(p->lex[0] == '#') {
        while(p->lex[0] && p->lex[0] != '\n') p->lex++;
        goto space;
    } else if(isalpha(p->lex[0])) {
        do {
            p->word[l++] = p->lex[0];
            p->lex++;
            if(l >= sizeof p->word - 1)
                error(p, "identifier too long\n");
        } while(isgraph(p->lex[0]));
        p->word[l] = '\0';
        if(p->mode == MODE_OBJ) {
            int i;
            for(i=0; OBJ_Command_Names[i].name; i++) {
                if(!strcmp(p->word, OBJ_Command_Names[i].name))
                    SYMBOL(OBJ_Command_Names[i].t);
            }
        } else if(p->mode == MODE_MTL) {
            int i;
            for(i=0; MTL_Command_Names[i].name; i++) {
                if(!strcmp(p->word, MTL_Command_Names[i].name))
                    SYMBOL(MTL_Command_Names[i].t);
            }
        }
        SYMBOL(ID);
    } else if(isdigit(p->lex[0])) {
        do {
            p->word[l++] = p->lex[0];
            p->lex++;
            if(l >= sizeof p->word - 1)
                error(p, "number too long\n");
        } while(isdigit(p->lex[0]) || p->lex[0] == '.');
        p->word[l] = '\0';
        SYMBOL(NUM);
    }
    p->sym = *p->lex++;
    p->word[0] = p->sym; p->word[1] = '\0';
}

#define VEC_ADD_FUNCTION(NAME, TYPE)                                              \
TYPE *obj_add_ ## NAME(ObjMesh *obj) {                                            \
    if(obj->n ## NAME ## s == obj->a ## NAME ## s) {                              \
        obj->a ## NAME ## s += obj->a ## NAME ## s >> 1;                          \
        obj->NAME ## s = realloc(obj->NAME ## s,                                  \
                            obj->a ## NAME ## s * sizeof *obj->NAME ## s);        \
    }                                                                             \
    return &obj->NAME ## s[obj->n ## NAME ## s++];                                \
}

VEC_ADD_FUNCTION(vert, ObjVert)
VEC_ADD_FUNCTION(norm, ObjNorm)
VEC_ADD_FUNCTION(tex, ObjUVW)
VEC_ADD_FUNCTION(pspace, ObjUVW)
/*VEC_ADD_FUNCTION(face, ObjFace)*/
ObjFace *obj_add_face(ObjMesh *obj) {
    if(obj->nfaces == obj->afaces) {
        obj->afaces += obj->afaces >> 1;
        obj->faces = realloc(obj->faces, obj->afaces * sizeof *obj->faces);
    }
    ObjFace *f = &obj->faces[obj->nfaces++];
    f->g = obj->ngroups - 1;
    f->s = 0;
    return f;
}
VEC_ADD_FUNCTION(group, ObjGroup)
ObjLine *obj_add_line(ObjMesh *obj) {
    if(obj->nlines == obj->alines) {
        obj->alines += obj->alines >> 1;
        obj->lines = realloc(obj->lines, obj->alines * sizeof *obj->lines);
    }
    ObjLine *l = &obj->lines[obj->nlines++];
    l->a = 4;
    l->idx = calloc(l->a, sizeof *l->idx);
    l->n = 0;
    return l;
}
void obj_line_add_vtx(ObjLine *l, int i) {
    if(l->n == l->a) {
        l->a += l->a >> 1;
        l->idx = realloc(l->idx, l->a * sizeof *l->idx);
    }
    l->idx[l->n++] = i;
}

static double read_float(OBJ_Parser *p) {
    double sign = accept(p,'-') ? -1.0 : 1.0;
    expect(p,NUM);
    return sign * atof(p->lastword);
}

static ObjFaceVert *parse_face_vert(OBJ_Parser *p, ObjMesh *obj, ObjFaceVert *fv) {
    fv->vt = -1;
    fv->vn = -1;
    int neg = accept(p, '-');
    expect(p, NUM);
    fv->v = atoi(p->lastword);
    if(neg) {
        fv->v = obj->nverts - fv->v;
        if(fv->v < 0)
            error(p, "negative v index\n");
    } else
        fv->v--;
    assert(fv->v >= 0);/* unsupported at the moment */
    if(accept(p, '/')) {
        neg = accept(p, '-');
        if(accept(p, NUM)) {
            fv->vt = atoi(p->lastword);
            if(neg) {
                fv->vt = obj->ntexs - fv->vt;
                if(fv->vt < 0)
                    error(p, "negative vt index\n");
            } else
                fv->vt--;
        }
        if(accept(p, '/')) {
            neg = accept(p, '-');
            expect(p, NUM);
            fv->vn = atoi(p->lastword);
            if(neg) {
                fv->vn = obj->nnorms - fv->vn;
                if(fv->vn < 0)
                    error(p, "negative vn index\n");
            } else
                fv->vn--;
        }
    }
    return fv;
}

static void parse_element(OBJ_Parser *p, ObjMesh *obj) {
    if(accept(p, V)) {
        ObjVert *v = obj_add_vert(obj);
        v->x = read_float(p);
        v->y = read_float(p);
        v->z = read_float(p);
        if(p->sym == NUM) {
            v->w = read_float(p);
        } else
            v->w = 1.0;
        if(v->x < obj->xmin) obj->xmin = v->x;
        if(v->x > obj->xmax) obj->xmax = v->x;
        if(v->y < obj->ymin) obj->ymin = v->y;
        if(v->y > obj->ymax) obj->ymax = v->y;
        if(v->z < obj->zmin) obj->zmin = v->z;
        if(v->z > obj->zmax) obj->zmax = v->z;
    } else if(accept(p, VN)) {
        ObjNorm *vn = obj_add_norm(obj);
        vn->x = read_float(p);
        vn->y = read_float(p);
        vn->z = read_float(p);
        /* According to [wiki][], the vector is not necessarily normalized */
        double _len = 1.0/sqrt(vn->x*vn->x + vn->y*vn->y + vn->z*vn->z);
        vn->x *= _len;
        vn->y *= _len;
        vn->z *= _len;
    } else if(accept(p, VT)) {
        ObjUVW *vt = obj_add_tex(obj);
        vt->u = read_float(p);
        vt->v = read_float(p);
        if(p->sym == NUM) {
            vt->w = read_float(p);
        } else
            vt->w = 0.0;
    } else if(accept(p, VP)) {
        ObjUVW *vp = obj_add_pspace(obj);
        vp->u = read_float(p);
        vp->v = read_float(p);
        vp->w = read_float(p);
    } else if(accept(p, G)) {
        expect(p, ID);
        ObjGroup *g = obj_add_group(obj);
        g->name = strdup(p->lastword);
    } else if(accept(p, S)) {
        if(accept(p, ID)) {
            if(strcmp(p->lastword, "off"))
                error(p, "number or 'off' expected\n");
            p->s = 0;
        } else {
            expect(p, NUM);
            p->s = atoi(p->lastword);
        }
    } else if(accept(p, MTLLIB)) {
        expect(p, ID);
        /* TODO */
    } else if(accept(p, USEMTL)) {
        expect(p, ID);
        /* TODO */
    } else if(accept(p, L)) {
        ObjLine *l = obj_add_line(obj);
        while(accept(p, '-') || accept(p, NUM)) {
            unsigned int idx;
            if(p->lastsym == '-') {
                expect(p, NUM);
                idx = obj->nverts - atoi(p->lastword);
                if(idx < 0)
                    error(p, "negative index\n");
            } else {
                idx = atoi(p->lastword) - 1;
                if(idx >= obj->nverts)
                    error(p, "invalid index\n");
            }
            obj_line_add_vtx(l, idx);
        }
    } else if(accept(p, F)) {
        /* Wikipedia says that negative indices might be present. See #Relative_and_absolute_indices
        on the wiki page.
        I use negative indices to indicate that the texture/normal index
        is unspecified, so if a negative index is encountered, I fix it here. */
        int i;
        ObjFace *f = obj_add_face(obj);

        ObjFaceVert *first, *last;
        for(i = 0; i < 3; i++) {
            last = parse_face_vert(p, obj, &f->verts[i]);
        }
        first = &f->verts[0];

        /* If more than 3 vertices are present, it is a triangle fan...
        https://stackoverflow.com/a/23724231/115589 */
        while(p->sym != '\0' && p->sym != '\n') {
            f = obj_add_face(obj);
            f->verts[0] = *first;
            f->verts[1] = *last;
            last = parse_face_vert(p, obj, &f->verts[2]);
        }

    } else {
        /* Maybe I can just ignore commands I don't understand, as
        I do for the materials. */
        error(p, "Unexpected %s\n", typename(p->sym));
    }
    if(p->sym == '\0') return;
    expect(p, '\n');
}

static void parse(OBJ_Parser *p, ObjMesh *obj) {
    while(p->sym != '\0') {
        if(accept(p, '\n')) continue;
        parse_element(p, obj);
    }
}

static char *slurp(const char *fname) {
    FILE *f;
    long len;
    char *str;

#ifdef SAFE_C11
    errno_t err = fopen_s(&f, fname, "rb");
    if (err != 0)
        return NULL;
#else
    if(!(f = fopen(fname, "rb")))
        return NULL;
#endif

    fseek(f, 0, SEEK_END);
    len = ftell(f);
    rewind(f);

    if(!(str = malloc(len+2)))
        return NULL;
    if(fread(str, 1, len, f) != len) {
        free(str);
        return NULL;
    }

    fclose(f);
    str[len] = '\0';
    return str;
}

#define INIT(NAME)                                         \
    obj->a ## NAME = 4;                                    \
    obj->NAME = calloc(obj->a ## NAME, sizeof *obj->NAME); \
    obj->n ## NAME = 0;

ObjMesh *obj_create() {
    ObjMesh *obj = malloc(sizeof *obj);
    memset(obj, 0, sizeof *obj);

    INIT(verts)
    INIT(norms)
    INIT(texs)
    INIT(pspaces)
    INIT(faces)
    INIT(groups)
    INIT(lines)

    obj->xmin = DBL_MAX; obj->xmax = DBL_MIN;
    obj->ymin = DBL_MAX; obj->ymax = DBL_MIN;
    obj->zmin = DBL_MAX; obj->zmax = DBL_MIN;
    return obj;
}

static char Error_Buf[128];

const char *obj_last_error() {
    return Error_Buf;
}

ObjMesh *obj_load(const char *fname) {
    OBJ_Parser parser;
    char *text = slurp(fname);
    if(!text) {
        snprintf(Error_Buf, sizeof Error_Buf, "OBJ: couldn't open %s: %s\n", fname, strerror(errno));
        return NULL;
    }
    init_OBJ(&parser, text);

    ObjMesh *obj = obj_create();

    if(!setjmp(parser.buf)) {
        Error_Buf[0] = '\0';
        parse(&parser, obj);
    } else {
        strncpy(Error_Buf, parser.error_buf, sizeof Error_Buf);
        Error_Buf[sizeof Error_Buf - 1] = '\0';
        obj_free(obj);
        obj = NULL;
    }

    free(text);
    return obj;
}

void obj_free(ObjMesh *obj) {
    unsigned int i;
    free(obj->verts);
    free(obj->norms);
    free(obj->texs);
    free(obj->pspaces);
    free(obj->faces);
    for(i = 0; i < obj->ngroups; i++)
        free(obj->groups[i].name);
    for(i = 0; i < obj->nlines; i++)
        free(obj->lines[i].idx);
    free(obj->groups);
    free(obj);
}

static void _obj_out(ObjMesh *obj, FILE *o) {
    unsigned int i;
    int g = -1, s = 0;
    fprintf(o,"# x extents: %g %g  (%g)\n", obj->xmin, obj->xmax, obj->xmax - obj->xmin);
    fprintf(o,"# y extents: %g %g  (%g)\n", obj->ymin, obj->ymax, obj->ymax - obj->ymin);
    fprintf(o,"# z extents: %g %g  (%g)\n", obj->zmin, obj->zmax, obj->zmax - obj->zmin);
    fprintf(o,"# %d vertices\n", obj->nverts);
    for(i = 0; i < obj->nverts; i++) {
        ObjVert *v = obj->verts + i;
        fprintf(o,"v %g %g %g %g\n", v->x, v->y, v->z, v->w);
    }
    if(obj->nnorms) {
        fprintf(o,"# %d normals\n", obj->nnorms);
        for(i = 0; i < obj->nnorms; i++) {
            ObjNorm *vn = obj->norms + i;
            fprintf(o,"vn %g %g %g\n", vn->x, vn->y, vn->z);
        }
    }
    if(obj->ntexs) {
        fprintf(o,"# %d texture coords\n", obj->ntexs);
        for(i = 0; i < obj->ntexs; i++) {
            ObjUVW *vt = obj->texs + i;
            fprintf(o,"vt %g %g %g\n", vt->u, vt->v, vt->w);
        }
    }
    if(obj->npspaces) {
        fprintf(o,"# %d param space verts\n", obj->npspaces);
        for(i = 0; i < obj->npspaces; i++) {
            ObjUVW *vp = obj->pspaces + i;
            fprintf(o,"vp %g %g %g\n", vp->u, vp->v, vp->w);
        }
    }
    fprintf(o,"# %d faces; %d groups\n", obj->nfaces, obj->ngroups);
    for(i = 0; i < obj->nfaces; i++) {
        int j;
        ObjFace *f = obj->faces + i;
        if(f->g != g) {
            g = f->g;
            fprintf(o, "g %s\n", obj->groups[g].name);
        }
        if(f->s != s) {
            s = f->s;
            fprintf(o, "s %d\n", s);
        }
        fputs("f ", o);
        for(j = 0; j < 3; j++) {
            if(f->verts[j].vt >= 0) {
                if(f->verts[j].vn >= 0) {
                    fprintf(o,"%d/%d/%d ", f->verts[j].v + 1, f->verts[j].vt + 1, f->verts[j].vn + 1);
                } else {
                    fprintf(o,"%d/%d ", f->verts[j].v + 1, f->verts[j].vt + 1);
                }
            } else {
                if(f->verts[j].vn >= 0) {
                    fprintf(o,"%d//%d ", f->verts[j].v + 1, f->verts[j].vn + 1);
                } else {
                    fprintf(o,"%d ", f->verts[j].v + 1);
                }
            }
        }
        fputc('\n', o);
    }
    if(obj->nlines) {
        fprintf(o,"# %d lines\n", obj->nlines);
        for(i = 0; i < obj->nlines; i++) {
            unsigned int j;
            ObjLine *l = &obj->lines[i];
            if(!l->n) continue;
            fputs("l ", o);
            for(j = 0; j < l->n; j++) {
                fprintf(o, "%d", l->idx[j] + 1);
                if(j < l->n-1) fputc(' ', o);
            }
            fputc('\n', o);
        }
    }
}

int obj_save(ObjMesh *obj, const char *fname) {
    FILE *f;
#ifdef SAFE_C11
    errno_t err = fopen_s(&f, fname, "w");
    if (err != 0) {
        snprintf(Error_Buf, sizeof Error_Buf, "OBJ: couldn't open %s: %s\n", fname, strerror(errno));
        return 0;
    }
#else
    f = fopen(fname, "w");
    if (!f) {
        snprintf(Error_Buf, sizeof Error_Buf, "OBJ: couldn't open %s: %s\n", fname, strerror(errno));
        return 0;
    }
#endif
    _obj_out(obj, f);
    fclose(f);
    return 1;
}

MtlLibrary *mtl_create() {
    MtlLibrary *lib = malloc(sizeof *lib);
    memset(lib, 0, sizeof *lib);
    lib->n = 0;
    lib->a = 4;
    lib->mtls = calloc(lib->a, sizeof *lib->mtls);
    return lib;
}
void mtl_free(MtlLibrary *lib) {
    unsigned int i;
    for(i = 0; i < lib->n; i++) {
        Material *m = &lib->mtls[i];
        free(m->name);
        if(m->Ka.type == mtl_spectral) free(m->Ka.spec.name);
        if(m->Kd.type == mtl_spectral) free(m->Kd.spec.name);
        if(m->Ks.type == mtl_spectral) free(m->Ks.spec.name);
        if(m->Tf.type == mtl_spectral) free(m->Tf.spec.name);
    }
    free(lib->mtls);
    free(lib);
}
Material *mtl_add(MtlLibrary *lib, const char *name) {
    if(lib->n == lib->a) {
        lib->a += lib->a >> 1;
        lib->mtls = realloc(lib->mtls, lib->a * sizeof *lib->mtls);
    }
    Material *m = &lib->mtls[lib->n++];
    m->name = _strdup(name);
    return m;
}
static MtlColor parse_mtl_color(OBJ_Parser *p) {
    MtlColor col;
    col.type = mtl_rgb;
    col.rgb.r = 1.0;
    col.rgb.g = 1.0;
    col.rgb.b = 1.0;

    if(accept(p, ID)) {
        // FIXME: Expect `spectral` or `xyz`
    } else {
        col.rgb.r = read_float(p);
        col.rgb.g = read_float(p);
        col.rgb.b = read_float(p);
    }
    return col;
}
static void mtl_parse(OBJ_Parser *p, MtlLibrary *lib) {
    Material *mtl = NULL;
    while(p->sym != '\0') {
        if(accept(p, '\n')) continue;
        if(accept(p, NEWMTL)) {
            expect(p, ID);
            mtl = mtl_add(lib, p->lastword);
        } else {
            if(!mtl)
                error(p, "`newmtl` expected, got `%s`\n", typename(p->sym));
            if(accept(p, KA))
                mtl->Ka = parse_mtl_color(p);
            else if(accept(p, KD))
                mtl->Kd = parse_mtl_color(p);
            else if(accept(p, KS))
                mtl->Ks = parse_mtl_color(p);
            else if(accept(p, TF))
                mtl->Tf = parse_mtl_color(p);
            else {
#if 1
                /* Just ignore commands you don't understand */
                while(p->sym != '\n' && p->sym != '\0')
                    nextsym(p);
#else
                error(p, "Unknown command `%s`\n", p->word);
#endif
            }
        }
    }
}
MtlLibrary *mtl_load(const char *fname) {
    OBJ_Parser parser;
    char *text = slurp(fname);
    if(!text) {
        snprintf(Error_Buf, sizeof Error_Buf, "MTL: couldn't open %s: %s\n", fname, strerror(errno));
        return NULL;
    }
    init_MTL(&parser, text);

    MtlLibrary *lib = mtl_create();

    if(!setjmp(parser.buf)) {
        Error_Buf[0] = '\0';
        mtl_parse(&parser, lib);
    } else {
        strncpy(Error_Buf, parser.error_buf, sizeof Error_Buf);
        Error_Buf[sizeof Error_Buf - 1] = '\0';
        mtl_free(lib);
        lib = NULL;
    }

    free(text);
    return lib;
}
static void col_out(const char *name, MtlColor col, FILE *out) {
    if(col.type == mtl_rgb) {
        fprintf(out, "%s %g %g %g\n", name, col.rgb.r, col.rgb.g, col.rgb.b);
    } else if(col.type == mtl_spectral) {
        fprintf(out, "%s spectral %s %g\n", name, col.spec.name, col.spec.factor);
    } else if(col.type == mtl_xyz) {
        fprintf(out, "%s xyz %g %g %g\n", name, col.xyz.x, col.xyz.y, col.xyz.z);
    }
}
static void _mtl_out(MtlLibrary *lib, FILE *out) {
    unsigned int i;
    for(i = 0; i < lib->n; i++) {
        Material *mtl = &lib->mtls[i];
        fprintf(out, "newmtl %s\n", mtl->name);
        col_out("Ka", mtl->Ka, out);
        col_out("Kd", mtl->Kd, out);
        col_out("Ks", mtl->Ks, out);
        col_out("Tf", mtl->Tf, out);
    }
}
int mtl_save(MtlLibrary *lib, const char *fname) {
    FILE *f = fopen(fname, "w");
    if(!f) {
        snprintf(Error_Buf, sizeof Error_Buf, "MTL: couldn't open %s: %s\n", fname, strerror(errno));
        return 0;
    }
    _mtl_out(lib, f);
    fclose(f);
    return 1;
}
#if defined(TEST)
/* gcc -g -Wall -DTEST obj.c */
int main(int argc, char *argv[]) {
    if(argc < 2) return 0;

    ObjMesh *obj = obj_load(argv[1]);
    if(!obj) {
        fputs(obj_last_error(), stderr);
        fprintf(stderr, "Unable to load OBJ %s", argv[1]);
        return 1;
    }
    _obj_out(obj, stdout);
    obj_free(obj);
    return 0;
}
#elif defined(MTL_TEST)
/* gcc -g -Wall -DMTL_TEST obj.c */
int main(int argc, char *argv[]) {
    if(argc < 2) return 0;

    MtlLibrary *lib = mtl_load(argv[1]);
    if(!lib) {
        fputs(obj_last_error(), stderr);
        fprintf(stderr, "Unable to load MTL %s", argv[1]);
        return 1;
    }
    _mtl_out(lib, stdout);
    mtl_free(lib);
    return 0;
}
#endif
