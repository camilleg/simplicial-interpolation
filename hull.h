/*
 * Ken Clarkson wrote this.  Copyright (c) 1995 by AT&T..
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#pragma once
#include <array>
#include <cstdio>
#include "points.h"
#include "stormacs.h"

#define MAXDIM 8
#define BLOCKSIZE 10000
#define MAXBLOCKS 3
#define DEBUG 0
#define CHECK_OVERSHOOT 1

extern char tmpfilenam[];

extern short check_overshoot_f;

FILE* efopen(const char *, const char *);

extern FILE *DFILE;

#define DEBS(qq)  {if (DEBUG>qq) {
#define EDEBS }}
#define DEBOUT DFILE
#define DEB(ll,mes)  DEBS(ll) fprintf(DEBOUT,#mes "\n");fflush(DEBOUT); EDEBS
#define DEBEXP(ll,exp) DEBS(ll) fprintf(DEBOUT,#exp "=%G\n", (double) exp); fflush(DEBOUT); EDEBS
#define DEBTR(ll) DEBS(ll) fprintf(DEBOUT, __FILE__ " line %d \n" ,__LINE__);fflush(DEBOUT); EDEBS
#define warning(lev, x) 					\
	{static int messcount;					\
		if (++messcount<=10) {DEB(lev,x) DEBTR(lev)}	\
		if (messcount==10) DEB(lev, consider yourself warned) \
	}							\

#define SBCHECK(s) /*								\
{double Sb_check=0;								\
int i;										\
	for (i=1;i<cdim;i++) if (s->neigh[i].basis)				\
				Sb_check+=s->neigh[i].basis->sqb;		\
	if ((float)(Sb_check - s->Sb) !=0.0)					\
	{DEBTR DEB(bad Sb); DEBEXP(s->Sb) DEBEXP(Sb_check);print_simplex(s); exit(1);}}*/

typedef point site;

extern site hull_p; 			/* the current site */

extern Coord hull_infinity[10];	/* point at infinity for Delaunay triang */

extern int
	rdim,	/* region dimension: (max) number of sites specifying region */
	cdim,	/* number of sites currently specifying region */
	site_size, /* size of malloc needed for a site */
	point_size;  /* size of malloc needed for a point */

typedef struct basis_s {
	struct basis_s *next; /* free list */
	int ref_count;	/* storage management */
	int lscale;    /* the log base 2 of total scaling of vector */
	Coord sqa, sqb; /* sums of squared norms of a part and b part */
	Coord vecs[1]; /* the actual vectors, extended by malloc'ing bigger */
} basis_s;
STORAGE_GLOBALS(basis_s)

typedef struct neighbor {
	site vert; /* vertex of simplex */
	struct simplex *simp; /* neighbor sharing all vertices but vert */
	basis_s *basis; /* derived vectors */
} neighbor;

typedef struct simplex {
	struct simplex *next;	/* free list */
	long visit;		/* number of last site visiting this simplex */
/*	float Sb; */
	short mark;
	basis_s* normal;	/* normal vector pointing inward */
	neighbor peak;		/* if null, remaining vertices give facet */
	neighbor neigh[1];	/* neighbors of simplex */
} simplex;
STORAGE_GLOBALS(simplex)

typedef struct fg_node fg;
typedef struct tree_node Tree;
struct tree_node {
    Tree *left, *right;
    site key;
    int size;   /* maintained to be the number of nodes rooted here */
    fg *fgs;
    Tree *next; /* freelist */
};
STORAGE_GLOBALS(Tree)

typedef struct fg_node {
	Tree *facets;
	double dist, vol;	/* of Voronoi face dual to this */
	fg *next;  		/* freelist */
	short mark;
	int ref_count;
} fg_node;
STORAGE_GLOBALS(fg)

typedef simplex* visit_func(simplex *, void *);
typedef int test_func(simplex *, int, void *);
typedef void out_func(point *, int, FILE*, int);

/* from driver, e.g., hullmain.c */

typedef site gsitef(void);
extern gsitef *get_site;	

typedef long site_n(site);
extern site_n *site_num;

extern double mult_up;

extern Coord mins[MAXDIM], maxs[MAXDIM];

typedef short zerovolf(simplex *);

extern double Huge;

/* from segt.c or ch.c */

simplex *build_convex_hull(gsitef*, site_n*, short, short);
void free_hull_storage(void);
int sees(site, simplex *);
void get_normal(simplex *s);
int out_of_flat(simplex*, site);
void set_ch_root(simplex*);
void print_site(site, FILE*);
void print_normal(simplex*);
visit_func check_marks;

double find_alpha(simplex*);
test_func alph_test;
simplex* visit_outside_ashape(simplex*, visit_func);

void get_basis_sede(simplex *);

	/* for debugging */
int check_perps(simplex*);

void find_volumes(fg*, FILE*);

#define MAXPOINTS 100
extern short mi[MAXPOINTS], mo[MAXPOINTS];

/* from hull.c */

simplex *visit_triang_gen(simplex *, visit_func, test_func);
simplex *visit_triang(simplex *, visit_func);
simplex* visit_hull(simplex *, visit_func);

neighbor *op_simp(simplex *a, simplex *b);

neighbor *op_vert(simplex *a, site b);

simplex *new_simp(void);

void buildhull(simplex *);

/* from io.c */

void panic(const char *fmt, ...);

typedef void print_neighbor_f(FILE*, neighbor*);
extern print_neighbor_f print_neighbor_full, print_neighbor_snum;

void check_triang(simplex*);

void check_new_triangs(simplex *);

void print_extra_facets(void);

simplex *print_facet(FILE*, simplex*, print_neighbor_f*);

void print_basis(FILE*, basis_s*);

simplex *print_simplex_f(simplex*, FILE*, print_neighbor_f*);

simplex *print_simplex(simplex*, void*);

void print_triang(simplex*, FILE*, print_neighbor_f*);

out_func CG_vlist_out, vlist_out, ps_out, cpr_out, mp_out, off_out;
	/* functions for different formats */

visit_func facets_print, afacets_print, ridges_print;
	/* to print facets, alpha facets, ridges */

void print_edge_dat(fg *, FILE *);

/* from pointops.c */

void print_point(FILE*, int, point);
void print_point_int(FILE*, int, point);
Coord maxdist(int,point p1, point p2);

/* from rand.c */

double double_rand(void);
void init_rand(long seed);


/* from fg.c, for face graphs */

fg *build_fg(simplex*);
void print_fg(fg*, FILE *);
void print_fg_alt(fg*, FILE *, int);
void print_hist_fg(simplex *, fg*, FILE *);

/*  void arena_check(void); */	/* from hobby's debugging malloc  */

// Derived from vss's map/mc.h.
using TT = std::array<int, MAXDIM+2>; // d+2 would suffice, but hull.exe doesn't know d.
using HH = std::array<int, MAXDIM+1>; // d+1 would suffice.
