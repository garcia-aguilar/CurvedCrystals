/// Declarations of globals and common objects

#ifndef COMMON_H_INCLUDED
#define COMMON_H_INCLUDED

///// geodesic_algorithm  ---       // include only in init_geometry -> problems with definitions
/*#include <iostream>                  // inside .h geodesic_algorithm files
#include <fstream>
#include <string>
*//*
#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"*/
///// ---- geodesic_algorithm

///// interface.h  ---
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ** random function ran2() --
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-EPS)
//**/ --- random ran2()

#define PI 3.14159265359
#define EPS 1.2e-7
//#define DT_tolerance 0.45           // Fix tolerance for time-step estimation
#define LINESIZE 130                // max amount of characters read "on a line"
#define MAX_SIZE 180000             // max number of vertices, triangles on the mesh
#define MAX_NEIGHBORS 500
///// ---- interface.h

#define OBTUSE_tolerance 30.        // allows only this % of obtuse angles in mesh
#define TIME_DATA_POINTS 1000       //  hi.dat save this amount of time info points
#define FIELD_CONF_FILES 10         // save this amount of s-t*dat configuration files
#define DIVERGENCE_STOP 1e60        // consider the integration algorithm to diverge if energy reaches this value

//#define VORTEX_EXP 2.               // in dislocation density approximation
//#define VORTEX_COS 6.               // in dislocation density approximation

#define  sign(x) ((x < 0) ? -1 : 1)
#define  min(x,y) ((x < y) ? x : y)
#define  max(x,y) ((x > y) ? x : y)
#define  mod(x,y) ((x % y + y) % y)

#define INPUT_FILE input_parameters.txt

/*******************************************************/

// *** Structure declarations  **

typedef struct {

	int label;

	long num_of_neighbors;
	long neighbor[MAX_NEIGHBORS];  //list of the vertex number of the neighboring vertices

  	double x;
	double y;
	double z;
	double r;

	double hx;
	double hy;
	double hz;
	double h2;
	double kg;

	double nx;
	double ny;
	double nz;

	double phi;

	double area;
	double weight[MAX_NEIGHBORS];

	double S;       // screened disclinicity from defect deformation
	double vorticity;   // dislocation density **+++ =0.0 if only disloc-screened charge
	double distToDefect[50];

} Vertex;

typedef struct {

	long v1;
	long v2;
	long v3;

} Triangle;

typedef struct {

	double x;
	double y;
	double z;

	double q;
	float qOff;

	double core_r;
	long closest_v;
	//long num_of_connections;
    //long connection[MAX_CONNECTIONS];

} Defect;

typedef struct {

  	long d;
  	long h;
  	long m;
  	long s;

} CPU_Time;

/*******************************************************/

// *** Mesh variables used throughout the code **
//extern long num_of_meshpoint;
extern Vertex vertex[MAX_SIZE];
extern double total_area;
extern double total_volume;
extern double dl_min;

/*******************************************************/

// *** File variables used throughout the code **
extern FILE *f_ou;


/*******************************************************/

// *** Integration parameters read from input **
extern double init_stress;          // initial total pre-stress
extern double phi_range;            // set the initial stress at each point randomly within this range
extern long seed;                   // for random fx to set initial stress value
//extern long num_of_iteration;
//extern double run_time;             // integration time step set accordingly  // change to read DT**
extern double time_step;            // integration time step
extern int integration_method;      // 1. Euler   2. Runge-Kutta (2 steps)  3. Runge-Kutta (4 steps) -> [bug]
extern double convergence_stop;     // stop integration when num_iterations completed or reaching this error

/*******************************************************/

// *** integration variables used throughout the code **
extern double convergence;

/*******************************************************/
// ***** main functions **
void init_geometry();                       // 3     in: init_geometry.cpp
//void run_field_integration();               // 5     in: run_field_integration.cpp
//void end(time_t, time_t);                   // 6     in: common.cpp
void get_time(time_t, CPU_Time *);          // 6.1   in: common.cpp

// **** Common functions   **
void export_conf(long, long);                               // 5.2.2
void progress_bar(long, long);                              // 3.6.3
int progress_bar_on_ou(long, long, int);                    // 3.6.3b



#endif // COMMON_H_INCLUDED
