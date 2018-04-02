#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <assert.h>

#include "obj.h"

static int Sym = 0, LastSym = 0;
static char Word[64], LastWord[64];
static char *Text = NULL;
static char *Lex = NULL;
int Line = 1;

enum { V = 'A', VT, VN, VP, F, L } Types;

static void nextsym();
void error(const char *msg, ...) {
    va_list ap;
    va_start(ap, msg);
    fprintf(stderr, "\nerror:%d: ", Line);
    vfprintf(stderr, msg, ap);
    va_end(ap);
    exit(EXIT_FAILURE);
}

static int accept(int s) {
    if(Sym == s) {
        nextsym();
        return 1;
    }
    return 0;
}

static int expect(int s) {
    if(accept(s))
        return 1;
    error("%d expected", s);
    return 0;
}

static void init(char *text) {
    Text = text;
    Lex = text;
    Line = 1;
    nextsym();
}

#define SYMBOL(s) do{Sym=s;return;}while(0)
static void nextsym() {
    int l = 0;
    LastSym = Sym;
    memcpy(LastWord, Word, sizeof LastWord);
    Word[0] = 0;
space:
    while(isspace(Lex[0])) {
        if(Lex[0] == '\n') {
            Line++;
            Lex++;
            SYMBOL('\n');
        }
        Lex++;
    }
    if(!Lex[0])
        SYMBOL('\0');
    else if(Lex[0] == '#') {
        while(Lex[0] && Lex[0] != '\n') Lex++;
        goto space;
    } else if(isalpha(Lex[0])) {
        do {
            Word[l++] = Lex[0];
            Lex++;
            if(l >= sizeof Word - 1)
                error("identifier too long");
        } while(isalnum(Lex[0]));
        Word[l] = '\0';
        if(!strcmp(Word, "v"))  SYMBOL(V);
        if(!strcmp(Word, "vt")) SYMBOL(VT);
        if(!strcmp(Word, "vn")) SYMBOL(VN);
        if(!strcmp(Word, "f"))  SYMBOL(F);
        if(!strcmp(Word, "l"))  SYMBOL(L);
        if(!strcmp(Word, "vp"))  SYMBOL(VP);
        error("unsupported type '%s'", Word);
    } else if(isdigit(Lex[0])) {
        do {
            Word[l++] = Lex[0];
            Lex++;
            if(l >= sizeof Word - 1)
                error("number too long");
        } while(isdigit(Lex[0]) || Lex[0] == '.');
        Word[l] = '\0';
        SYMBOL('0');
    }
    Sym = *Lex++;
    Word[0] = Sym; Word[1] = '\0';
}

static char *readfile(const char *fname) {
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

#define VEC_ADD_FUNCTION(NAME, TYPE)                                          \
static TYPE *add_ ## NAME(ObjMesh *obj) {                                     \
    if(obj->n ## NAME == obj->a ## NAME) {                                    \
        obj->a ## NAME += obj->a ## NAME >> 1;                                \
        obj->NAME = realloc(obj->NAME, obj->a ## NAME * sizeof *obj->NAME);   \
    }                                                                         \
    return &obj->NAME[obj->n ## NAME++];                                      \
}

VEC_ADD_FUNCTION(verts, ObjVert)
VEC_ADD_FUNCTION(norms, ObjNorm)
VEC_ADD_FUNCTION(texs, ObjUVW)
VEC_ADD_FUNCTION(pspaces, ObjUVW)
VEC_ADD_FUNCTION(faces, ObjFace)

static double read_float() {
    double sign = accept('-') ? -1.0 : 1.0;
    expect('0');
    return sign * atof(LastWord);
}

static void parse_element(ObjMesh *obj) {
    if(accept(V)) {
        ObjVert *v = add_verts(obj);
        v->x = read_float();
        v->y = read_float();
        v->z = read_float();
        if(Sym == '0') {
            v->w = read_float();
        } else
            v->w = 1.0;
        if(v->x < obj->xmin) obj->xmin = v->x;
        if(v->x > obj->xmax) obj->xmax = v->x;
        if(v->y < obj->ymin) obj->ymin = v->y;
        if(v->y > obj->ymax) obj->ymax = v->y;
        if(v->z < obj->zmin) obj->zmin = v->z;
        if(v->z > obj->zmax) obj->zmax = v->z;
    } else if(accept(VN)) {
        ObjNorm *vn = add_norms(obj);
        vn->x = read_float();
        vn->y = read_float();
        vn->z = read_float();
        // TODO: Wikipedia says the normal vector is not necessarily a unit vector
    } else if(accept(VT)) {
        ObjUVW *vt = add_texs(obj);
        vt->u = read_float();
        vt->v = read_float();
        if(Sym == '0') {
            vt->w = read_float();
        } else
            vt->w = 0.0;
    } else if(accept(VP)) {
        ObjUVW *vp = add_pspaces(obj);
        vp->u = read_float();
        vp->v = read_float();
        vp->w = read_float();
    } else if(accept(L)) {
        while(accept('0')) {
            // TODO: something something atoi(LastWord)
        }
    } else if(accept(F)) {
        /* TODO: Wikipedia says that negative indices might be present.
        https://en.wikipedia.org/wiki/Wavefront_.obj_file#Relative_and_absolute_indices
        I use negative indices to indicate that the texture/normal index is unspecified,
        which may lead to incompatibilities.
        */
        int i;
        ObjFace *f = add_faces(obj);
        for(i = 0; i < 3; i++) {
            f->verts[i].vt = -1;
            f->verts[i].vn = -1;
            expect('0');
            f->verts[i].v = atoi(LastWord);
            assert(f->verts[i].v >= 0);/* unsupported at the moment */
            if(accept('/')) {
                if(accept('0'))
                    f->verts[i].vt = atoi(LastWord);
                if(accept('/')) {
                    expect('0');
                    f->verts[i].vn = atoi(LastWord);
                }
            }
        }
    } else {
        error("Unexpected %d", Sym);
    }
    if(Sym == '\0') return;
    expect('\n');
}

static void parse(ObjMesh *obj) {
    while(Sym != '\0') {
        if(accept('\n')) continue;
        parse_element(obj);
    }
}

#define INIT(NAME)                                         \
    obj->a ## NAME = 2;                                    \
    obj->NAME = calloc(obj->a ## NAME, sizeof *obj->NAME); \
    obj->n ## NAME = 0;

ObjMesh *obj_load(const char *fname) {
    char *text = readfile(fname);
    if(!text)
        return NULL;
    init(text);

    ObjMesh *obj = malloc(sizeof *obj);
    memset(obj, 0, sizeof *obj);

    INIT(verts)
    INIT(norms)
    INIT(texs)
    INIT(pspaces)
    INIT(faces)

    obj->xmin = DBL_MAX; obj->xmax = DBL_MIN;
    obj->ymin = DBL_MAX; obj->ymax = DBL_MIN;
    obj->zmin = DBL_MAX; obj->zmax = DBL_MIN;

    parse(obj);

    free(text);
    return obj;
}

void obj_free(ObjMesh *obj) {
    free(obj->verts);
    free(obj->norms);
    free(obj->texs);
    free(obj->pspaces);
    free(obj->faces);
    free(obj);
}

static void _obj_out(ObjMesh *obj, FILE *o) {
    int i;
    fprintf(o,"# x extents: %g %g  (%g)\n", obj->xmin, obj->xmax, obj->xmax - obj->xmin);
    fprintf(o,"# y extents: %g %g  (%g)\n", obj->ymin, obj->ymax, obj->ymax - obj->ymin);
    fprintf(o,"# z extents: %g %g  (%g)\n", obj->zmin, obj->zmax, obj->zmax - obj->zmin);
    fprintf(o,"# %d vertices\n", obj->nverts);
    for(i = 0; i < obj->nverts; i++) {
        ObjVert *v = obj->verts + i;
        fprintf(o,"v %g %g %g %g\n", v->x, v->y, v->z, v->w);
    }
    fprintf(o,"# %d normals\n", obj->nnorms);
    for(i = 0; i < obj->nnorms; i++) {
        ObjNorm *vn = obj->norms + i;
        fprintf(o,"vn %g %g %g\n", vn->x, vn->y, vn->z);
    }
    fprintf(o,"# %d texture coords\n", obj->ntexs);
    for(i = 0; i < obj->ntexs; i++) {
        ObjUVW *vt = obj->texs + i;
        fprintf(o,"vt %g %g %g\n", vt->u, vt->v, vt->w);
    }
    fprintf(o,"# %d param space verts\n", obj->npspaces);
    for(i = 0; i < obj->npspaces; i++) {
        ObjUVW *vp = obj->pspaces + i;
        fprintf(o,"vp %g %g %g\n", vp->u, vp->v, vp->w);
    }
    fprintf(o,"# %d faces\n", obj->nfaces);
    for(i = 0; i < obj->nfaces; i++) {
        int j;
        ObjFace *f = obj->faces + i;
        fprintf(o,"f ");
        for(j = 0; j < 3; j++) {
            if(f->verts[j].vt >= 0) {
                if(f->verts[j].vn >= 0) {
                    fprintf(o,"%d/%d/%d ", f->verts[j].v, f->verts[j].vt, f->verts[j].vn);
                } else {
                    fprintf(o,"%d/%d ", f->verts[j].v, f->verts[j].vt);
                }
            } else {
                if(f->verts[j].vn >= 0) {
                    fprintf(o,"%d//%d ", f->verts[j].v, f->verts[j].vn);
                } else {
                    fprintf(o,"%d ", f->verts[j].v);
                }
            }
        }
        fprintf(o,"\n");
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
int main(int argc, char *argv[]) {
    if(argc < 2) return 0;

    ObjMesh *obj = obj_load(argv[1]);
    if(!obj) {
        error("Unable to load OBJ %s", argv[1]);
        return 1;
    }
    _obj_out(obj, stdout);
    obj_free(obj);
    return 0;
}
#endif
