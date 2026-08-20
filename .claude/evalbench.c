/* evalbench.c -- how long does one (value, gradient) evaluation of a PLQ conjugate cost?
 *
 * The question is SCIP's: it calls a nonlinear constraint's eval/grad millions of times, so the
 * per-call cost decides whether a CCA2 conjugate can be a constraint at all. Three strategies,
 * all in double precision, all returning value AND gradient:
 *
 *   A  mesh + linear scan        locate the cell by testing cells in order (early exit per cell)
 *   B  mesh + O(1) bucket        a uniform grid over the bounding box names 1..3 candidate cells
 *   C  no mesh, per-piece sup    f*(s) = max_k sup_{x in T_k} <s,x> - q_k(x), closed form
 *   D  baseline                  a small expression tree walked like SCIPexprEval does
 *
 * Every cell predicate is evaluated with the FULL conic form (a x^2 + b xy + c y^2 + d x + e y + f)
 * even when the constraint is affine, so A and B are upper bounds on the real per-constraint cost.
 *
 * Build:  gcc -O2 -march=native -o evalbench evalbench.c -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define REPS 3

static double now_s(void) {
    return (double)clock() / (double)CLOCKS_PER_SEC;
}

/* ---------------------------------------------------------------- mesh strategies (A, B) */
typedef struct { double a, b, c, d, e, f; } conic;      /* a x^2 + b xy + c y^2 + d x + e y + f */
typedef struct { conic g[4]; double q[6]; } cell;       /* 4 constraints <= 0, one quadratic face */

static inline double conic_at(const conic *g, double x, double y) {
    return ((g->a * x + g->b * y + g->d) * x) + ((g->c * y + g->e) * y) + g->f;
}

/* value and gradient of the face quadratic q = [A B C D E F] : A x^2 + B xy + C y^2 + D x + E y + F */
static inline double face_at(const double *q, double x, double y, double *gx, double *gy) {
    *gx = 2.0 * q[0] * x + q[1] * y + q[3];
    *gy = q[1] * x + 2.0 * q[2] * y + q[4];
    return ((q[0] * x + q[1] * y + q[3]) * x) + ((q[2] * y + q[4]) * y) + q[5];
}

static cell *mesh;
static int   NCELL, GRID;                 /* NCELL = GRID*GRID boxes over [0,1]^2 */

static void build_mesh(int g) {
    GRID = g; NCELL = g * g;
    mesh = malloc(sizeof(cell) * NCELL);
    double h = 1.0 / g;
    for (int i = 0; i < g; i++) for (int j = 0; j < g; j++) {
        cell *C = &mesh[i * g + j];
        double x0 = i * h, x1 = (i + 1) * h, y0 = j * h, y1 = (j + 1) * h;
        /* four affine constraints, stored (and evaluated) in the conic form */
        C->g[0] = (conic){0,0,0, -1, 0,  x0};   /* x0 - x <= 0 */
        C->g[1] = (conic){0,0,0,  1, 0, -x1};   /* x - x1 <= 0 */
        C->g[2] = (conic){0,0,0,  0,-1,  y0};
        C->g[3] = (conic){0,0,0,  0, 1, -y1};
        for (int k = 0; k < 6; k++) C->q[k] = 0.25 + 0.1 * ((i * 7 + j * 13 + k) % 5);
    }
}

static double eval_scan(double x, double y, double *gx, double *gy) {
    for (int c = 0; c < NCELL; c++) {
        const cell *C = &mesh[c];
        if (conic_at(&C->g[0], x, y) > 0) continue;
        if (conic_at(&C->g[1], x, y) > 0) continue;
        if (conic_at(&C->g[2], x, y) > 0) continue;
        if (conic_at(&C->g[3], x, y) > 0) continue;
        return face_at(C->q, x, y, gx, gy);
    }
    *gx = *gy = 0; return 0;
}

static double eval_bucket(double x, double y, double *gx, double *gy, int cand) {
    int i = (int)(x * GRID); if (i < 0) i = 0; if (i >= GRID) i = GRID - 1;
    int j = (int)(y * GRID); if (j < 0) j = 0; if (j >= GRID) j = GRID - 1;
    int base = i * GRID + j;
    for (int t = 0; t < cand; t++) {                 /* cand-1 near misses, then the hit */
        int c = (base + NCELL - t) % NCELL;
        const cell *C = &mesh[c];
        if (conic_at(&C->g[0], x, y) > 0) continue;
        if (conic_at(&C->g[1], x, y) > 0) continue;
        if (conic_at(&C->g[2], x, y) > 0) continue;
        if (conic_at(&C->g[3], x, y) > 0) continue;
        return face_at(C->q, x, y, gx, gy);
    }
    const cell *C = &mesh[base];
    return face_at(C->q, x, y, gx, gy);
}

/* ------------------------------------------------------- strategy C: no mesh, per-piece sup */
typedef struct { double V[3][2]; double Q[3]; double be[2]; double ga; } piece;  /* Q = [q11 q12 q22] */
static piece *pieces;
static int NPIECE;

static void build_pieces(int n) {
    NPIECE = n;
    pieces = malloc(sizeof(piece) * n);
    for (int k = 0; k < n; k++) {
        double t = 2.0 * M_PI * k / n;
        pieces[k].V[0][0] = 0;              pieces[k].V[0][1] = 0;
        pieces[k].V[1][0] = cos(t);         pieces[k].V[1][1] = sin(t);
        pieces[k].V[2][0] = cos(t + 1.0);   pieces[k].V[2][1] = sin(t + 1.0);
        pieces[k].Q[0] = 1.0 + 0.1 * k; pieces[k].Q[1] = 0.3; pieces[k].Q[2] = 2.0 - 0.05 * k;
        pieces[k].be[0] = 0.2 * k; pieces[k].be[1] = -0.1 * k;
        pieces[k].ga = 0.5;
    }
}

/* sup over the triangle of <s,x> - q(x), with q(x) = 1/2 x'Qx + be'x + ga.  The objective is
 * concave here (Q positive definite), so: the unconstrained critical point when it is inside,
 * otherwise the max over the three edges, each a 1-D quadratic on a segment -- closed form, and
 * the argmax is the gradient (Danskin). */
static double piece_sup(const piece *P, double sx, double sy, double *ax, double *ay) {
    double q11 = P->Q[0], q12 = P->Q[1], q22 = P->Q[2];
    double r0 = sx - P->be[0], r1 = sy - P->be[1];
    double det = q11 * q22 - q12 * q12;
    double best = -INFINITY, bx = 0, by = 0;
    /* unconstrained critical point x* = Q^{-1} r */
    double cx = ( q22 * r0 - q12 * r1) / det;
    double cy = (-q12 * r0 + q11 * r1) / det;
    /* inside test: barycentric against the three edges */
    const double (*V)[2] = P->V;
    int inside = 1;
    for (int e = 0; e < 3; e++) {
        double x0 = V[e][0], y0 = V[e][1], x1 = V[(e+1)%3][0], y1 = V[(e+1)%3][1];
        double x2 = V[(e+2)%3][0], y2 = V[(e+2)%3][1];
        double s1 = (x1-x0)*(cy-y0) - (y1-y0)*(cx-x0);
        double s2 = (x1-x0)*(y2-y0) - (y1-y0)*(x2-x0);
        if (s1 * s2 < 0) { inside = 0; break; }
    }
    if (inside) {
        double v = sx*cx + sy*cy - (0.5*(q11*cx*cx + 2*q12*cx*cy + q22*cy*cy)
                                    + P->be[0]*cx + P->be[1]*cy + P->ga);
        best = v; bx = cx; by = cy;
    } else {
        for (int e = 0; e < 3; e++) {
            double x0 = V[e][0], y0 = V[e][1];
            double dx = V[(e+1)%3][0] - x0, dy = V[(e+1)%3][1] - y0;
            /* phi(t) = <s, x0+t d> - q(x0+t d) : concave 1-D quadratic, maximise on [0,1] */
            double A = -0.5 * (q11*dx*dx + 2*q12*dx*dy + q22*dy*dy);
            double B = (sx*dx + sy*dy) - (q11*x0*dx + q12*(x0*dy + y0*dx) + q22*y0*dy)
                       - (P->be[0]*dx + P->be[1]*dy);
            double C = sx*x0 + sy*y0 - (0.5*(q11*x0*x0 + 2*q12*x0*y0 + q22*y0*y0)
                                        + P->be[0]*x0 + P->be[1]*y0 + P->ga);
            double t = (A < 0) ? -B / (2*A) : (B > 0 ? 1.0 : 0.0);
            if (t < 0) t = 0; if (t > 1) t = 1;
            double v = (A*t + B)*t + C;
            if (v > best) { best = v; bx = x0 + t*dx; by = y0 + t*dy; }
        }
    }
    *ax = bx; *ay = by;
    return best;
}

static double eval_pieces(double sx, double sy, double *gx, double *gy) {
    double best = -INFINITY, bx = 0, by = 0;
    for (int k = 0; k < NPIECE; k++) {
        double ax, ay;
        double v = piece_sup(&pieces[k], sx, sy, &ax, &ay);
        if (v > best) { best = v; bx = ax; by = ay; }
    }
    *gx = bx; *gy = by;                   /* Danskin: the argmax IS the gradient */
    return best;
}

/* ----------------------------------------------------- baseline D: a small expression tree */
enum { OP_CONST, OP_VARX, OP_VARY, OP_ADD, OP_MUL, OP_SQR };
typedef struct { int op, l, r; double val; } node;
static node tree[24];
static int NNODE;

static void build_tree(void) {
    /* ~20 nodes: a quadratic-ish expression, the size SCIP sees for a small nonlinear constraint */
    int n = 0;
    tree[n++] = (node){OP_VARX,0,0,0};        /* 0 */
    tree[n++] = (node){OP_VARY,0,0,0};        /* 1 */
    tree[n++] = (node){OP_CONST,0,0,2.5};     /* 2 */
    tree[n++] = (node){OP_MUL,0,1,0};         /* 3  x*y */
    tree[n++] = (node){OP_SQR,0,0,0};         /* 4  x^2 */
    tree[n++] = (node){OP_SQR,1,0,0};         /* 5  y^2 */
    tree[n++] = (node){OP_MUL,2,4,0};         /* 6 */
    tree[n++] = (node){OP_ADD,6,3,0};         /* 7 */
    tree[n++] = (node){OP_ADD,7,5,0};         /* 8 */
    tree[n++] = (node){OP_MUL,8,0,0};         /* 9 */
    tree[n++] = (node){OP_ADD,9,8,0};         /* 10 */
    tree[n++] = (node){OP_MUL,10,1,0};        /* 11 */
    tree[n++] = (node){OP_ADD,11,3,0};        /* 12 */
    tree[n++] = (node){OP_ADD,12,4,0};        /* 13 */
    tree[n++] = (node){OP_MUL,13,2,0};        /* 14 */
    tree[n++] = (node){OP_ADD,14,5,0};        /* 15 */
    tree[n++] = (node){OP_ADD,15,0,0};        /* 16 */
    tree[n++] = (node){OP_MUL,16,16,0};       /* 17 */
    tree[n++] = (node){OP_ADD,17,12,0};       /* 18 */
    tree[n++] = (node){OP_ADD,18,2,0};        /* 19 */
    NNODE = n;
}

static double eval_tree(double x, double y, double *gx, double *gy) {
    double v[24], dx[24], dy[24];              /* value + forward-mode gradient, as SCIP does */
    for (int i = 0; i < NNODE; i++) {
        const node *t = &tree[i];
        switch (t->op) {
            case OP_CONST: v[i] = t->val; dx[i] = 0; dy[i] = 0; break;
            case OP_VARX:  v[i] = x; dx[i] = 1; dy[i] = 0; break;
            case OP_VARY:  v[i] = y; dx[i] = 0; dy[i] = 1; break;
            case OP_ADD:   v[i] = v[t->l] + v[t->r]; dx[i] = dx[t->l] + dx[t->r];
                           dy[i] = dy[t->l] + dy[t->r]; break;
            case OP_MUL:   v[i] = v[t->l] * v[t->r];
                           dx[i] = dx[t->l]*v[t->r] + v[t->l]*dx[t->r];
                           dy[i] = dy[t->l]*v[t->r] + v[t->l]*dy[t->r]; break;
            case OP_SQR:   v[i] = v[t->l]*v[t->l]; dx[i] = 2*v[t->l]*dx[t->l];
                           dy[i] = 2*v[t->l]*dy[t->l]; break;
        }
    }
    *gx = dx[NNODE-1]; *gy = dy[NNODE-1];
    return v[NNODE-1];
}

/* ------------------------------------------------------------------------------- driver */
static unsigned long long rs = 88172645463325252ULL;
static inline double rnd(void) { rs ^= rs<<13; rs ^= rs>>7; rs ^= rs<<17; return (rs>>11)*(1.0/9007199254740992.0); }

#define TIMEIT(NAME, CALL, NITER)                                        \
    do {                                                                 \
        double bestt = 1e30;                                             \
        for (int r = 0; r < REPS; r++) {                                 \
            rs = 88172645463325252ULL;                                    \
            double acc = 0, t0 = now_s();                                \
            for (long it = 0; it < (NITER); it++) {                      \
                double X = rnd(), Y = rnd(), gx, gy;                     \
                acc += CALL; acc += gx + gy;                             \
            }                                                            \
            double el = now_s() - t0;                                    \
            if (el < bestt) bestt = el;                                  \
            if (acc == 12345.6789) printf("");                           \
        }                                                                \
        printf("%-34s %8.1f ns/eval\n", NAME, 1e9 * bestt / (NITER));    \
    } while (0)

int main(void) {
    char buf[64];
    build_tree();
    printf("--- baseline: a ~20-node expression tree, value + forward-mode gradient\n");
    TIMEIT("D  expression tree (20 nodes)", eval_tree(X, Y, &gx, &gy), 5000000L);

    printf("\n--- mesh strategies: cell location, then one quadratic face\n");
    int gs[] = {3, 6, 9, 16, 32};             /* 9, 36, 81, 256, 1024 cells */
    for (int i = 0; i < 5; i++) {
        build_mesh(gs[i]);
        snprintf(buf, sizeof buf, "A  linear scan, %4d cells", NCELL);
        TIMEIT(buf, eval_scan(X, Y, &gx, &gy), 2000000L);
        free(mesh);
    }
    for (int i = 0; i < 5; i++) {
        build_mesh(gs[i]);
        snprintf(buf, sizeof buf, "B  bucket + 1 candidate, %4d cells", NCELL);
        TIMEIT(buf, eval_bucket(X, Y, &gx, &gy, 1), 5000000L);
        snprintf(buf, sizeof buf, "B  bucket + 3 candidates, %4d cells", NCELL);
        TIMEIT(buf, eval_bucket(X, Y, &gx, &gy, 3), 5000000L);
        free(mesh);
    }

    printf("\n--- no mesh: max over pieces of a closed-form 2-D sup (gradient = argmax)\n");
    int ps[] = {3, 6, 12, 24};
    for (int i = 0; i < 4; i++) {
        build_pieces(ps[i]);
        snprintf(buf, sizeof buf, "C  per-piece sup, %3d pieces", NPIECE);
        TIMEIT(buf, eval_pieces(X, Y, &gx, &gy), 2000000L);
        free(pieces);
    }
    return 0;
}
