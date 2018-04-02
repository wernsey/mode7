#ifndef OBJ_H
#define OBJ_H

/*
 * OBJ
 * ===
 *
 * Reads [Wavefront .OBJ][wiki] 3D models.
 *
 * ## Links:
 *
 *  * [Wavefront .obj file][wiki]
 *  * [Object Files (.obj)][bourke]
 *  * <https://www.cs.cmu.edu/~mbz/personal/graphics/obj.html>
 *
 * ## TODO:
 *
 *  * Groups with the `g` command
 *  * Materials: `mtllib`, `usemtl` etc.
 *
 * [wiki]: https://en.wikipedia.org/wiki/Wavefront_.obj_file
 * [bourke]: http://paulbourke.net/dataformats/obj/
 */

typedef struct {
    double x, y, z, w;
} ObjVert;

typedef struct {
    double x, y, z;
} ObjNorm;

typedef struct {
    double u, v, w;
} ObjUVW;

typedef struct {
    /* Technically, OBJ supports more
    than 3 vertices per face, but I don't. */
    struct {
        int v, vt, vn;
    } verts[3];
} ObjFace;

typedef struct ObjMesh {

    unsigned int nverts, averts;
    ObjVert *verts; /* Vertices */

    unsigned int nnorms, anorms;
    ObjNorm *norms; /* Normals */

    unsigned int ntexs, atexs;
    ObjUVW *texs;  /* Texture coordinates */

    unsigned int npspaces, apspaces;
    ObjUVW *pspaces;  /* Parameter space vertices */

    unsigned int nfaces, afaces;
    ObjFace *faces; /* Faces */

    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;

} ObjMesh;

ObjMesh *obj_load(const char *fname);
void obj_free(ObjMesh *obj);
void obj_out(ObjMesh *obj, const char *fname);

#endif
