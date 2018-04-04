#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <assert.h>

#include "obj.h"

typedef struct {
    int sym,
        lastsym;
    char word[64],
        lastword[64];
    char *text;
    char *lex;
    int line;
    int s; /* smoothing group */
} OBJ_Parser;

enum { V = 128, VT, VN, VP, F, L, G, S, MTLLIB, USEMTL, ID, NUM } Types;

static void nextsym(OBJ_Parser *p);

static void error(OBJ_Parser *p, const char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    fprintf(stderr, "\nerror:%d: ", p->line);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
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
    error(p,"%d expected", s);
    return 0;
}

static void init(OBJ_Parser *p, char *text) {
    p->sym = 0;
    p->word[0] = 0;
    p->lastsym = 0;
    p->lastword[0] = 0;
    p->text = text;
    p->lex = text;
    p->line = 1;
    p->s = 0;
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
                error(p, "identifier too long");
        } while(isgraph(p->lex[0]));
        p->word[l] = '\0';
        if(!strcmp(p->word, "v"))  SYMBOL(V);
        if(!strcmp(p->word, "vt")) SYMBOL(VT);
        if(!strcmp(p->word, "vn")) SYMBOL(VN);
        if(!strcmp(p->word, "f"))  SYMBOL(F);
        if(!strcmp(p->word, "l"))  SYMBOL(L);
        if(!strcmp(p->word, "vp"))  SYMBOL(VP);
        if(!strcmp(p->word, "g"))  SYMBOL(G);
        if(!strcmp(p->word, "s"))  SYMBOL(S);
        if(!strcmp(p->word, "usemtl"))  SYMBOL(USEMTL);
        if(!strcmp(p->word, "mtllib"))  SYMBOL(MTLLIB);
        SYMBOL(ID);
    } else if(isdigit(p->lex[0])) {
        do {
            p->word[l++] = p->lex[0];
            p->lex++;
            if(l >= sizeof p->word - 1)
                error(p, "number too long");
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
            error(p, "negative v index");
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
                    error(p, "negative vt index");
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
                    error(p, "negative vn index");
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
                error(p, "number or 'off' expected");
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
            int idx;
            if(p->lastsym == '-') {
                expect(p, NUM);
                idx = obj->nverts - atoi(p->lastword);
                if(idx < 0)
                    error(p, "negative index");
            } else {
                idx = atoi(p->lastword) - 1;
                if(idx >= obj->nverts)
                    error(p, "invalid index");
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
        error(p, "Unexpected %d", p->sym);
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

static char *_obj_readf(const char *fname) {
    FILE *f;
    long len;
    char *str;

    if(!(f = fopen(fname, "rb")))
        return NULL;

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

ObjMesh *obj_load(const char *fname) {
    OBJ_Parser parser;
    char *text = _obj_readf(fname);
    if(!text)
        return NULL;
    init(&parser, text);

    ObjMesh *obj = obj_create();

    parse(&parser, obj);

    free(text);
    return obj;
}

void obj_free(ObjMesh *obj) {
    int i;
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
    int i, g = -1, s = 0;
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
            int j;
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

void obj_out(ObjMesh *obj, const char *fname) {
    FILE *f = fopen(fname, "w");
    if(f) {
        _obj_out(obj, f);
        fclose(f);
    }
}

#ifdef TEST
/* gcc -g -Wall -DTEST obj.c */
int main(int argc, char *argv[]) {
    if(argc < 2) return 0;

    ObjMesh *obj = obj_load(argv[1]);
    if(!obj) {
        fprintf(stderr, "Unable to load OBJ %s", argv[1]);
        return 1;
    }
    _obj_out(obj, stdout);
    obj_free(obj);
    return 0;
}
#endif
