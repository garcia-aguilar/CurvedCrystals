/******************************************************************
**                                                               **
**                           CALCULATE_DISTANCES_I               **
**                                                               **
** Program: calculate_distances_i.cpp                            **
** Version: 2.1                                                  **
**     - > from .msh and defs.tmp                                **
**        get distance from input vertex i to all defect         **
**        with a chosen algorithm                                **
**        and store info in i/d/dist file                        **
** (with common.h)                                               **
** Last Change: Sep 11 2018                                      **
**                                                               **
*******************************************************************/
#include <iostream>
#include <fstream>
#include <string>

#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"

#include "common.h"

time_t t1, t2;

Triangle triangle[MAX_SIZE];
Vertex vertex[MAX_SIZE];
Defect defect[MAX_SIZE];    //+++ defs

long num_of_meshpoint;
long num_of_triangles;

int distance_method;
int mesh_number;

long num_of_defects;     //+++ defs
long num_of_obtuse;
//double dl_min;        // in COMMON.CPP
double dl_max;
double max_DT;

double R0_s=1.0;
double calculated_radius;
//double total_area;        // in COMMON.CPP
//double total_volume;      // in COMMON.CPP
double asphericity = 0;

double max_dZ = 0.0;
int displace = 1;

char f_name[32];
char fD_name[32];   //+++ defs
char f_dist_name[32];
int p_crystal;

void init(int, char **);
void import_mesh(char *);
long vertex_i;

void get_geometry();
void find_neighbors();


void read_defects(char *);
void calculate_defect_distance();
template<class Points, class Faces>
void read_mesh_from_current_data(Points&, Faces&);
double great_circle_distance(long, long);
void write_triMshd(FILE *);
void write_triMshd_head(FILE *);
void write_triMshd_tail(FILE *);

//void progress_bar(long, long);            // in COMMON.CPP
//void get_time(time_t, CPU_Time *);        // in COMMON.CPP

void end_calculate_distances();

/*******************************************************************/

main(int argc, char** argv)
{

    time(&t1);

   // START
	init(argc,argv);

    get_geometry();
    read_defects(fD_name);

    calculate_defect_distance();

	end_calculate_distances();

	fclose(f_ou);

  	return EXIT_SUCCESS;
}

/*******************************************************************/
void init(int argCount, char **argList)         // change to C++ program arguments
{
    if(argCount != 7)
    {
        printf("\n\n Input Params can be passed as program command-line arguments (6 values).\n\n Here input by hand:\n");
      // choose algorithm for surface distance calculations
        printf("  Choose distance calculation method \n\t1.Circle\t2.Dijkstra\t3.Exact\t\t5.Euclidean\n");
        printf("  Input (int):  ");
        scanf("%d", &distance_method);

      // .msh file name
        printf("  Input mesh file:   ");
        scanf("%s",f_name);
        printf("  Mesh Number:   ");
        scanf("%d", &mesh_number);

      // Add defect
        printf("  Input defect file:  ");
        scanf("%s",fD_name);

      // Which meshpoint?
        printf("  Input vertex number:  ");
        scanf("%ld",&vertex_i);

      // File Name for distances
        printf("  Input distances file:  ");
        scanf("%s",f_dist_name);
    }

    if(argCount == 7)
    {
        char *dumm;
        printf("\n\n Reading input params as: 1.DistanceMethod, 2.MeshFileName, 3.MeshID, 4.DefectFileName 5.Vertex_indx 6.DistFileName\n\n");
        distance_method =  strtol (argList[1],&dumm,10);//(int) argList[1]-'0';
        strcpy(f_name, argList[2]);
        mesh_number = strtol (argList[3],&dumm,10);
        strcpy(fD_name,argList[4]);
        vertex_i = strtol (argList[5],&dumm,0);
        strcpy(f_dist_name,argList[6]);
    }

   // Write ou file with params and progress
    std::string f_ou_name;                  // f_ou in common.cpp
    f_ou_name = "ou_ico"+std::to_string(mesh_number)+"_"+std::to_string(distance_method)+"_i"+std::to_string(vertex_i)+".dat";
    //printf("%s",f_ou_name.c_str());
    f_ou = fopen(f_ou_name.c_str(), "w");

    fprintf(f_ou, "<<< Input Parameters >>>\n");
    fprintf(f_ou,"MeshFile  %s\n", f_name);
    fprintf(f_ou,"MeshNumber  %d\n", mesh_number);
    fprintf(f_ou,"DistanceMethod  %d\n", distance_method);
    fprintf(f_ou,"\t1.Circle\t2.Dijkstra\t3.Exact\t\t5.Euclidean\n");
    fprintf(f_ou,"DefectFile  %s\n", fD_name);
    fprintf(f_ou,"MeshpointIndx  %ld\n", vertex_i);
    fflush(f_ou);

    import_mesh(f_name);
}


/*******************************************************************/

void import_mesh(char *f_name)
{
	long i, v1, v2, v3, chi, dummy, num_of_edges=0;


	FILE *f_in = fopen(f_name,"r");


	char line[LINESIZE];

	// Read line by line on the file (or up to LINESIZE characts)
	// Go pass the first 8
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);

	fscanf(f_in,"%ld",&num_of_meshpoint);

	if (num_of_meshpoint>MAX_SIZE){
		printf("Error: the number of vertices (%ld) exceeds MAX_SIZE (%d)\n",num_of_meshpoint,MAX_SIZE);
		exit(0);
	}

	for (i=0; i<num_of_meshpoint; i++){
		fscanf(f_in,"%ld%lg%lg%lg",
		&dummy,
		&vertex[i].x,
		&vertex[i].y,
		&vertex[i].z);
	}

	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);
	fgets(line,LINESIZE-1,f_in);

	fscanf(f_in,"%ld",&num_of_triangles);

	if (num_of_triangles>MAX_SIZE){
		printf("Error: the number of triangles (%ld) exceeds MAX_SIZE (%d)\n",num_of_triangles,MAX_SIZE);
		exit(0);
	}

	for (i=0; i<num_of_triangles; i++){
		fscanf(f_in,"%ld%ld%ld%ld%ld%ld%ld%ld",
		&dummy,
		&dummy,
		&dummy,
		&dummy,
		&dummy,
		&v1,
		&v2,
		&v3);
		triangle[i].v1 = v1-1;
		triangle[i].v2 = v2-1;
		triangle[i].v3 = v3-1;
	}
	fclose(f_in);

	find_neighbors();

	for (i=0; i<num_of_meshpoint; i++){
		num_of_edges += vertex[i].num_of_neighbors;
	}

	if (!num_of_edges%2){
		printf("Error: bad triangulation, 2E = %ld\n",num_of_edges);
		exit(0);
	}
    // Eurler characteristic
	chi = num_of_meshpoint-num_of_edges/2+num_of_triangles;
	printf("Euler characteristic = %ld\n",chi);         //+++ verf
	//printf("meshpoint = %ld\n",num_of_meshpoint);         //+++ verf
	//printf("triangles = %ld\n",num_of_triangles);         //+++ verf
	//printf("edge= %ld\n",num_of_edges);         //+++ verf

    // change for topologies with more handles
	if (chi!=2){
		printf("Error: bad triangulation, chi = %ld\n",chi);
		exit(0);
	}
	printf("\tVertices %ld\n",num_of_meshpoint);
	printf("\tTriangles %ld\n",num_of_triangles);
}

/*******************************************************************/

void find_neighbors()
{
	int is_neighbor(long,long);
	long i, max_neighbors=0;

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].num_of_neighbors = 0;
	}

	for (i=0; i<num_of_triangles; i++){
		if (!is_neighbor(triangle[i].v1,triangle[i].v2)){
			vertex[triangle[i].v1].neighbor[vertex[triangle[i].v1].num_of_neighbors] = triangle[i].v2;
			vertex[triangle[i].v1].num_of_neighbors++;
			if (vertex[triangle[i].v1].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v1;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v1,triangle[i].v3)){
			vertex[triangle[i].v1].neighbor[vertex[triangle[i].v1].num_of_neighbors] = triangle[i].v3;
			vertex[triangle[i].v1].num_of_neighbors++;
			if (vertex[triangle[i].v1].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v1;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v2,triangle[i].v1)){
			vertex[triangle[i].v2].neighbor[vertex[triangle[i].v2].num_of_neighbors] = triangle[i].v1;
			vertex[triangle[i].v2].num_of_neighbors++;
			if (vertex[triangle[i].v2].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v2;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v2,triangle[i].v3)){
			vertex[triangle[i].v2].neighbor[vertex[triangle[i].v2].num_of_neighbors] = triangle[i].v3;
			vertex[triangle[i].v2].num_of_neighbors++;
			if (vertex[triangle[i].v2].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v2;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v3,triangle[i].v1)){
			vertex[triangle[i].v3].neighbor[vertex[triangle[i].v3].num_of_neighbors] = triangle[i].v1;
			vertex[triangle[i].v3].num_of_neighbors++;
			if (vertex[triangle[i].v3].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v3;
				break;
			}
		}
		if (!is_neighbor(triangle[i].v3,triangle[i].v2)){
			vertex[triangle[i].v3].neighbor[vertex[triangle[i].v3].num_of_neighbors] = triangle[i].v2;
			vertex[triangle[i].v3].num_of_neighbors++;
			if (vertex[triangle[i].v3].num_of_neighbors==MAX_NEIGHBORS){
				max_neighbors = triangle[i].v3;
				break;
			}
		}
	}

	// If any vertex gets to the max number of neighbors
	if (max_neighbors){
		printf("Error: vertex %ld reached the maximum number of neighbors\n",max_neighbors);
		exit(0);
	}
}

/*******************************************************************/

int is_neighbor(long i, long j)
{
	long k;

	for (k=0; k<vertex[i].num_of_neighbors; k++){
		if (j==vertex[i].neighbor[k]) return 1;
	}

	return 0;
}

/*******************************************************************/
void get_geometry()
{
	long t, i, j, k, n;

	double theta_i,theta_j,theta_k,cot_i,cot_j,cot_k;
	double x[3],y[3],z[3], dx, dy, base, height,side_i2,side_j2;
	double norm, area_t,mean_radius, radial_deviation;

	long obtuse_t[MAX_SIZE];
	float obtuse_tolerance = 30.0;    //+++ admit only this % of obtuse angles in msh
	float obtuse_percent;

	//// Control
	num_of_obtuse = 0;
	total_area = 0;
	total_volume =0;
	dl_min = 1E5;
	dl_max = 0.0;
	mean_radius = 0;
	radial_deviation = 0;

	// Initialize vector quantities

	for (i=0; i<num_of_meshpoint; i++){
        vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].h2 = 0;
        vertex[i].nx = 0;
		vertex[i].ny = 0;
		vertex[i].nz = 0;
		vertex[i].area = 0;
		vertex[i].S = 0;
		vertex[i].r = sqrt(vertex[i].x*vertex[i].x+vertex[i].y*vertex[i].y+vertex[i].z*vertex[i].z);
		mean_radius += vertex[i].r;    // Asphericity calculations
	}
	mean_radius /= num_of_meshpoint;

    // Loop through triangles to get geometry info for each vertex

	for (t=0; t<num_of_triangles; t++){

        i = triangle[t].v1;
        j = triangle[t].v2;
        k = triangle[t].v3;

        // Construct an approximation for the tangent plane at each vertex

		x[0] = vertex[j].x-vertex[i].x;
		x[1] = vertex[j].y-vertex[i].y;
		x[2] = vertex[j].z-vertex[i].z;

		base = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

		x[0] /= base;
		x[1] /= base;
		x[2] /= base;

		y[0] = vertex[k].x-vertex[i].x;
		y[1] = vertex[k].y-vertex[i].y;
		y[2] = vertex[k].z-vertex[i].z;

		norm = x[0]*y[0]+x[1]*y[1]+x[2]*y[2];

		y[0] -= norm*x[0];
		y[1] -= norm*x[1];
		y[2] -= norm*x[2];

		height = sqrt(y[0]*y[0]+y[1]*y[1]+y[2]*y[2]);

		y[0] /= height;
		y[1] /= height;
		y[2] /= height;

        dx = (vertex[k].x-vertex[i].x)*x[0]
            +(vertex[k].y-vertex[i].y)*x[1]
            +(vertex[k].z-vertex[i].z)*x[2];

        dy = (vertex[k].x-vertex[i].x)*y[0]
            +(vertex[k].y-vertex[i].y)*y[1]
            +(vertex[k].z-vertex[i].z)*y[2];


    // Initially, calculate surface normal vector from averaging over normals of 1-ring faces
        //un normalized cross product
        z[0] = x[1]*y[2]-x[2]*y[1];
        z[1] = x[2]*y[0]-x[0]*y[2];
        z[2] = x[0]*y[1]-x[1]*y[0];
        //check R(dot)N is positive (outwards for convex, with origin 'inside')
        norm = z[0]*vertex[i].x+z[1]*vertex[i].y+z[2]*vertex[i].z;
        if (norm < 0){
            //printf("an inwards N vector was outwards at i = %ld\n",i);    //+++++ CHK
            z[0] *= -1;
            z[1] *= -1;
            z[2] *= -1;
        }
        //vertex i
        vertex[i].nx += z[0];
        vertex[i].ny += z[1];
        vertex[i].nz += z[2];
        //vertex j
        vertex[j].nx += z[0];
        vertex[j].ny += z[1];
        vertex[j].nz += z[2];
        //vertex k
        vertex[k].nx += z[0];
        vertex[k].ny += z[1];
        vertex[k].nz += z[2];

        // Area triangle
        area_t = base*height/2;
        total_area += area_t;       //+++ CHK

        // Calculate angles

        theta_i = atan2(dy,dx);
        cot_i = dx/dy;

        dx = base-dx;       // if theta_i/j obtuse, base < dx, but sign differs
        theta_j = atan2(dy,dx);
        cot_j = dx/dy;
        theta_k = PI-theta_i-theta_j;
        cot_k = 1/tan(theta_k);

        // Besides base (side k), get length of other two sides
        side_j2 = (vertex[k].x-vertex[i].x)*(vertex[k].x-vertex[i].x)
                +(vertex[k].y-vertex[i].y)*(vertex[k].y-vertex[i].y)
                +(vertex[k].z-vertex[i].z)*(vertex[k].z-vertex[i].z);
        side_i2 = (vertex[k].x-vertex[j].x)*(vertex[k].x-vertex[j].x)
                +(vertex[k].y-vertex[j].y)*(vertex[k].y-vertex[j].y)
                +(vertex[k].z-vertex[j].z)*(vertex[k].z-vertex[j].z);

        // Find min and max vertex distance
        if(dl_min > base){ dl_min = base;}
        if(dl_max < base){ dl_max = base;}
        if(dl_min > sqrt(side_i2)){ dl_min = sqrt(side_i2);}
        if(dl_max < sqrt(side_i2)){ dl_max = sqrt(side_i2);}
        if(dl_min > sqrt(side_j2)){ dl_min = sqrt(side_j2);}
        if(dl_max < sqrt(side_j2)){ dl_max = sqrt(side_j2);}


        // Sum area contribution to vertices
        if (theta_i >= PI/2 ){
            vertex[i].area += area_t/2;
            vertex[j].area += area_t/4;
            vertex[k].area += area_t/4;
            obtuse_t[num_of_obtuse] = t;
            num_of_obtuse++;
        } else if (theta_j >= PI/2){
            vertex[i].area += area_t/4;
            vertex[j].area += area_t/2;
            vertex[k].area += area_t/4;
            obtuse_t[num_of_obtuse] = t;
            num_of_obtuse++;
        } else if (theta_k >= PI/2){
            vertex[i].area += area_t/4;
            vertex[j].area += area_t/4;
            vertex[k].area += area_t/2;
            obtuse_t[num_of_obtuse] = t;
            num_of_obtuse++;
        } else {
            vertex[i].area += (cot_j*side_j2+cot_k*base*base)/8;
            vertex[j].area += (cot_i*side_i2+cot_k*base*base)/8;
            vertex[k].area += (cot_i*side_i2+cot_j*side_j2)/8;
        }


        // Sum weighted vectors to compute mean curvature

        // vertex i
        vertex[i].hx += (cot_k*(vertex[i].x-vertex[j].x)+cot_j*(vertex[i].x-vertex[k].x))/4;
        vertex[i].hy += (cot_k*(vertex[i].y-vertex[j].y)+cot_j*(vertex[i].y-vertex[k].y))/4;
        vertex[i].hz += (cot_k*(vertex[i].z-vertex[j].z)+cot_j*(vertex[i].z-vertex[k].z))/4;
        // vertex j
        vertex[j].hx += (cot_k*(vertex[j].x-vertex[i].x)+cot_i*(vertex[j].x-vertex[k].x))/4;
        vertex[j].hy += (cot_k*(vertex[j].y-vertex[i].y)+cot_i*(vertex[j].y-vertex[k].y))/4;
        vertex[j].hz += (cot_k*(vertex[j].z-vertex[i].z)+cot_i*(vertex[j].z-vertex[k].z))/4;
        // vertex k
        vertex[k].hx += (cot_i*(vertex[k].x-vertex[j].x)+cot_j*(vertex[k].x-vertex[i].x))/4;
        vertex[k].hy += (cot_i*(vertex[k].y-vertex[j].y)+cot_j*(vertex[k].y-vertex[i].y))/4;
        vertex[k].hz += (cot_i*(vertex[k].z-vertex[j].z)+cot_j*(vertex[k].z-vertex[i].z))/4;

	}

	printf("\tTotal surface area from triangles area = %g\n",total_area); //++++ chk
    total_area = 0;

    // Loop through vertices to normalize over the area, and get integral quantities
    double max_z = -1.0;
	double min_z = 1.0;

    for (i=0;i<num_of_meshpoint;i++){

        //Geometry height
        if (vertex[i].z<min_z){
            min_z = vertex[i].z;
        }
        if (vertex[i].z > max_z){
            max_z = vertex[i].z;
        }

        // Area
        total_area += vertex[i].area;

        // Mean Curvature
        vertex[i].hx /= vertex[i].area;
        vertex[i].hy /= vertex[i].area;
        vertex[i].hz /= vertex[i].area;
        vertex[i].h2 = vertex[i].hx*vertex[i].hx+vertex[i].hy*vertex[i].hy+vertex[i].hz*vertex[i].hz;

        // Normal vectors (using H2 if 'high')
        if (vertex[i].h2 > 1E-20){
            norm = vertex[i].x*vertex[i].hx+vertex[i].y*vertex[i].hy+vertex[i].z*vertex[i].hz;
            //*** ONLY for ~ convex geometries, with origin inside the closed surface:
            if(norm > 0){
                vertex[i].nx = vertex[i].hx;
                vertex[i].ny = vertex[i].hy;
                vertex[i].nz = vertex[i].hz;
            }
            else if(norm < 0){
                //printf("an inwards H vector was turned to an outwards N vector at i = %ld\n",i);
                vertex[i].nx = -vertex[i].hx;
                vertex[i].ny = -vertex[i].hy;
                vertex[i].nz = -vertex[i].hz;
            } else {printf("R dot N = 0 at i = %ld\n",i);}
		}

		// normalize outwards surface vectors (then nec for both H and cross products)
		norm = sqrt(vertex[i].nx*vertex[i].nx+vertex[i].ny*vertex[i].ny+vertex[i].nz*vertex[i].nz);
        vertex[i].nx /= norm;
        vertex[i].ny /= norm;
        vertex[i].nz /= norm;

		// Volume
		total_volume += (vertex[i].x*vertex[i].nx+vertex[i].y*vertex[i].ny+vertex[i].z*vertex[i].nz)*vertex[i].area/3.;


        // Asphericity calculations
		double diff = vertex[i].r - mean_radius;
		radial_deviation += diff*diff;
    }

    // Asphericity
    asphericity = radial_deviation/(num_of_meshpoint*mean_radius*mean_radius);



    // Geometry Z height
    max_dZ  = max_z-min_z;
    calculated_radius = pow(total_volume*3./(4*PI),1./3);
    double vol_sphere = 4*PI*pow(calculated_radius,3)/3;
    double area_sphere = 4*PI*pow(calculated_radius,2);

	//printf("\tTotal volume with (x,y,z)/3 = %.10g   (sphere should be ~ %.6g)\n",total_volume, vol_sphere);    //+++ CHK


    obtuse_percent = 100.0*num_of_obtuse/num_of_triangles;

	//printf("\n");
	printf("\tVertices %ld\n",num_of_meshpoint);
	printf("\tTriangles %ld\n",num_of_triangles);
	printf("\tObtuse triangles %ld  -> %.3f%%\n",num_of_obtuse,obtuse_percent);
	printf("\tTotal surface area = %g  (sphere should be ~ %.5g)\n",total_area, area_sphere);
	printf("\tTotal volume = %g   (sphere should be ~ %.5g)\n",total_volume, vol_sphere);
	printf("\tAsphericity by radial variance is = %.8G  (sphere should be ~ 0)\n",asphericity);
	printf("\tMax dl = %.4G, Min dl = %.4G\n", dl_max, dl_min);   //++++
	//printf("\tDt < %.3E  (if cut .45)\n", DT_tolerance*dl_min*dl_min); //+++
	printf("\n");
	printf("\n\t Calculated Radius %.7lf",calculated_radius);      //+++++CHK
    printf("\n\t size difference(set-calculated) = %.4lG%%  (should be ~0)\n",(R0_s-calculated_radius)*100./R0_s);


	// Check is percentage of obstuse angles in msh less than tolerance
	if (obtuse_percent > obtuse_tolerance)
	{
        printf("Too many obtuse triangles\nPROBLEM\n");
        exit(1);
	}

}


/*******************************************************************/

void read_defects(char *f_name)
{
    long d;

    FILE *f_in = fopen(f_name,"r");
    char line[LINESIZE];

    // Read crystal ideal coordination number
    fgets(line,LINESIZE-1,f_in);    // $Crystal coordination
    fscanf(f_in,"%d\n",&p_crystal);

    // Read number of defects
    fgets(line,LINESIZE-1,f_in);    // $ Number of defects
    fscanf(f_in,"%ld\n",&num_of_defects);

    if (num_of_defects>MAX_SIZE){
		printf("Error: the number of defects (%ld) exceeds MAX_SIZE (%d)\n",num_of_defects,MAX_SIZE);
		exit(0);
	}

     // Read defect information
    fgets(line,LINESIZE-1,f_in);    // $ Defects
    fgets(line,LINESIZE-1,f_in);    //  info

    for (d =0; d<num_of_defects; d++){
        fscanf(f_in,"%lg %lg %lg %lg %g %lg",
        &defect[d].x,
        &defect[d].y,
        &defect[d].z,
        &defect[d].q,
        &defect[d].qOff,
        &defect[d].core_r);
    }

    fclose(f_in);

    printf("\n\n\tN_defects  %ld\n", num_of_defects);
}


/*******************************************************************/
void calculate_defect_distance()
{
    time(&t1);

    long i, d, i_min;
    double r[3],distance_to_defect, r2, min_dist;


    // Locate each defect (x,y,z) at the mesh vertex closest to the initial coordinates.
    // False if defect centers at exact tips of inscribed icosahedron
    //int displace = 1;
    //if(distance_method==1) displace =0;           //++ Always displace for evolver msh

    /*   printf("\n BEFOREDISPLACED \n");
    for (d=0;d<num_of_defects;d++)
    {
        printf("\n\t d=%ld  %.5f\t%.5f\t%.5f  --  %ld",d,defect[d].x,defect[d].y,defect[d].z,defect[d].closest_v);
    }
    printf("\n");*/    //CHECK pos reading ++RM

    for(d=0;d<num_of_defects;d++){
        min_dist = 1E6;
        i_min = 0;
        for(i=0; i<num_of_meshpoint; i++){
            r[0] = vertex[i].x-defect[d].x;
            r[1] = vertex[i].y-defect[d].y;
            r[2] = vertex[i].z-defect[d].z;
            r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            if(r2<min_dist){
                min_dist = r2;
                i_min = i;
            }
        }
        if(displace){
            defect[d].x = vertex[i_min].x;
            defect[d].y = vertex[i_min].y;
            defect[d].z = vertex[i_min].z;
        }
        defect[d].closest_v = i_min;
    }
    if(displace) {
        printf("\tDefects displaced\n");
    }else{printf("\tDefects NOT displaced\n");}

    /*
    printf("\n\n");
    for (d=0;d<num_of_defects;d++)
    {
        printf("\n\t d=%ld  %.5f\t%.5f\t%.5f  --  %ld",d,defect[d].x,defect[d].y,defect[d].z,defect[d].closest_v);
    }
    printf("\n");*/   // CHECK reposition ++RM


  // Write .mshd file with distances and triangle information
    std::string f_mshd_name;
    f_mshd_name = "tri_ico"+std::to_string(mesh_number)+"_"+std::to_string(distance_method)+".mshd";
    //FILE *f_mshd = fopen(f_mshd_name.c_str(), "w");           // per vertex, just append
    FILE *f_mshd = fopen(f_mshd_name.c_str(), "a");

  // Write .dat file with vertex / defect / distances info
    std::string f_dist_saveName;
    //f_mshd_name = "tri_"+std::string (f_name)+"d";
    f_dist_saveName = std::string(f_dist_name)+"_ico"+std::to_string(mesh_number)+"_"+std::to_string(distance_method)+".dat";
    FILE *f_dist = fopen(f_dist_saveName.c_str(), "a");

    // Variables needed for geodesic_algorithm
    std::vector<double> points;     // keep types because it is read by mesh.initialize_mesh_data
    std::vector<unsigned> faces;
    std::vector<geodesic::SurfacePoint> path;
    geodesic::Mesh mesh;
    unsigned source_vertex_index, target_vertex_index;

  //  Include progress bar on ou.dat to check remotely
    int progress_bar_count = 0;
    fprintf(f_ou, "\n ---- Distance Progress: [ 10 .");
    fflush(f_ou);

    // Choose how to calculate distance on surface
        // 1: great-circle distance on a sphere, 2:dijkstra or 3: exact geodesic (mesh) algortihm.
        // 4: read distances from file, 5: euclidean
    switch(distance_method){
    case 1: // circle
        printf("\n** Using great-circle distance**\n");

        //for (i=0; i<num_of_meshpoint; i++){
        fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", vertex_i, vertex[vertex_i].x,vertex[vertex_i].y,vertex[vertex_i].z);
        for(d=0;d<num_of_defects;d++){
            distance_to_defect = great_circle_distance(vertex_i,d);        //*****
            r2 = distance_to_defect*distance_to_defect;
            vertex[vertex_i].distToDefect[d] = distance_to_defect;   //LLLL   //RM
            fprintf(f_mshd,"\t%.11f", distance_to_defect);
            fprintf(f_dist,"%ld\t%ld\t%.11f\n", vertex_i, d, distance_to_defect);
            progress_bar(d, num_of_defects);
            progress_bar_count = progress_bar_on_ou(d, num_of_defects, progress_bar_count);
        }

        //}
        //write_triMshd(f_mshd);
        fprintf(f_mshd,"\n");
        break;

    case 2:  // dijkstra
    {
        printf("\n **Using Dijkstra algorithm to calculate geodesic distances**\n");

        // SAVE tri*mshd file Line-by-line
        //write_triMshd_head(f_mshd);  // per vertex, do by hand.

        read_mesh_from_current_data(points, faces);
        // + create internal mesh data
        mesh.initialize_mesh_data(points, faces);
        geodesic::GeodesicAlgorithmDijkstra dijkstra_algorithm(&mesh);

        //for (i=0; i<num_of_meshpoint; i++){
          // write to tri*mshd
        fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", vertex_i, vertex[vertex_i].x,vertex[vertex_i].y,vertex[vertex_i].z);
        for(d=0;d<num_of_defects;d++){
            source_vertex_index = defect[d].closest_v;
            target_vertex_index = vertex_i;
            geodesic::SurfacePoint source (&mesh.vertices()[source_vertex_index]);
            geodesic::SurfacePoint target (&mesh.vertices()[target_vertex_index]);
            dijkstra_algorithm.geodesic(source, target, path);
            distance_to_defect = geodesic::length(path);
            //std::cout <<"D" << d <<" i"<< i <<" "<< distance_to_defect << std::endl;

            r2 = distance_to_defect*distance_to_defect;
            vertex[vertex_i].distToDefect[d] = distance_to_defect;   //LLLL
            fprintf(f_mshd,"\t%.11f", distance_to_defect);
            fprintf(f_dist,"%ld\t%ld\t%.11f\n", vertex_i, d, distance_to_defect);
            progress_bar(d, num_of_defects);
            progress_bar_count = progress_bar_on_ou(d, num_of_defects, progress_bar_count);
        }


        fprintf(f_mshd,"\n");
        //}
       // write_triMshd_tail(f_mshd);       //per vertex, do by hand
        break;
    }

    case 3:  // exact geodesic_algortihm
    {
        printf("\n **Using Exact algorithm to calculate geodesic distances**\n");

        // SAVE tri*mshd file Line-by-line
        //write_triMshd_head(f_mshd);       //per vertex, do by hand

        read_mesh_from_current_data(points, faces);
        // + create internal mesh data
        mesh.initialize_mesh_data(points, faces);
        geodesic::GeodesicAlgorithmExact exact_algorithm(&mesh);

        //for (i=0; i<num_of_meshpoint; i++){
            // write to tri*mshd
        //fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", vertex_i, vertex[vertex_i].x,vertex[vertex_i].y,vertex[vertex_i].z);  //Print per calculation
        for(d=0;d<num_of_defects;d++){
            source_vertex_index = defect[d].closest_v;
            target_vertex_index = vertex_i;
            geodesic::SurfacePoint source (&mesh.vertices()[source_vertex_index]);
            geodesic::SurfacePoint target (&mesh.vertices()[target_vertex_index]);
            exact_algorithm.geodesic(source, target, path);
            distance_to_defect = geodesic::length(path);
            //printf("\ni %ld   d %ld  dist=  %.4f\n", i, d, distance_to_defect);

            r2 = distance_to_defect*distance_to_defect;
            vertex[vertex_i].distToDefect[d] = distance_to_defect;   //LLLL
            //fprintf(f_mshd,"\t%.11f", distance_to_defect);                //Print per calculation
            fprintf(f_dist,"%ld\t%ld\t%.11f\n", vertex_i, d, distance_to_defect);
            progress_bar(d, num_of_defects);
            progress_bar_count = progress_bar_on_ou(d, num_of_defects, progress_bar_count);
        }
        fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", vertex_i, vertex[vertex_i].x,vertex[vertex_i].y,vertex[vertex_i].z);   //print all
        for(d=0;d<num_of_defects;d++){              //print all
            fprintf(f_mshd,"\t%.11f", vertex[vertex_i].distToDefect[d]);        //print all
        }       //print all
        fprintf(f_mshd,"\n");
        fflush(f_mshd);


        //}
        //write_triMshd_tail(f_mshd);       //per vertex, do by hand
        break;
    }
    case 5:
    {
        printf("\n **Using euclidean distances**\n");
        //for (i=0; i<num_of_meshpoint; i++){
        fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", vertex_i, vertex[vertex_i].x,vertex[vertex_i].y,vertex[vertex_i].z);
        for(d=0;d<num_of_defects;d++){
            r[0] = vertex[vertex_i].x-defect[d].x;
            r[1] = vertex[vertex_i].y-defect[d].y;
            r[2] = vertex[vertex_i].z-defect[d].z;
            r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
            distance_to_defect = sqrt(r2);
            vertex[vertex_i].distToDefect[d] = distance_to_defect;
            fprintf(f_mshd,"\t%.11f", distance_to_defect);
            fprintf(f_dist,"%ld\t%ld\t%.11f\n", vertex_i, d, distance_to_defect);
        }
        //}
        //write_triMshd(f_mshd);    //per vertex, do by hand
        break;
    }
    }

    fclose(f_mshd);

}


/*******************************************************************/
double great_circle_distance(long i, long d){

    double r[3], acosArg, l_arc;

  // If defects are at exact meshpoint, return l_arc = 0.
    if (displace) {
        if (i==defect[d].closest_v){
            return 0.;
        }
    }

  //Distance from vertex to defect core (still, relatively small core radius, shortest dist on surface)
    r[0] = vertex[i].x*defect[d].x;
    r[1] = vertex[i].y*defect[d].y;
    r[2] = vertex[i].z*defect[d].z;
    acosArg = (r[0]+r[1]+r[2])/pow(calculated_radius,2.);

  // correct numerical precision errors
    if (acosArg <-1) acosArg = -1.;
    else if (acosArg>1) acosArg = 1.;

    l_arc = acos(acosArg)*calculated_radius;
    return l_arc;
}

/*******************************************************************/
template<class Points, class Faces>
void read_mesh_from_current_data(Points& points, Faces& faces){

    for (long i = 0; i<num_of_meshpoint; i++){
        points.push_back(vertex[i].x);
        points.push_back(vertex[i].y);
        points.push_back(vertex[i].z);
    }
    //std::cout << points.size()/3 << std::endl;


    for (long t = 0; t<num_of_triangles; t++){
        faces.push_back(triangle[t].v1);
        faces.push_back(triangle[t].v2);
        faces.push_back(triangle[t].v3);
    }
    //std::cout << faces.size()/3 << std::endl;
    unsigned maxE = *max_element(faces.begin(),faces.end());
    unsigned minE = *min_element(faces.begin(),faces.end());
    std::cout << "min and max vertex on a face "<< maxE <<" "<< minE << std::endl;

}


/*******************************************************************/
/*void write_ge(FILE *f_ou)
{
    long i,j;

	for (i=0; i<num_of_meshpoint; i++){
		fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10f\t%.10f\t%.10f\t",
			i,
			vertex[i].x,
			vertex[i].y,
			vertex[i].z,
			vertex[i].area,
			vertex[i].h2,
			vertex[i].kg,
            vertex[i].S,
            vertex[i].hx,
            vertex[i].hy,
            vertex[i].hz,
            vertex[i].nx,
            vertex[i].ny,
            vertex[i].nz);

        //LLLL

        for (j=0; j<12; j++){
            fprintf(f_ou,"%.8lg\t",vertex[i].distToDefect[j]);
        }           //LLLLLL*/
/*
		for (j=0; j<vertex[i].num_of_neighbors; j++){
			fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);
		}

		fprintf(f_ou,"\n");
	}
}*/



/*******************************************************************/

/*******************************************************************/
/*          // ALREADY in common.cpp
void progress_bar(long t, long total)
{
	long i, ticks, percent, num_of_ticks=20;
	char progress[32];

	percent = 100*t/total;
	ticks = num_of_ticks*percent/100;

	for (i=0; i<20; i++){
		progress[i] = ((i<=ticks)?'#':' ');
	}
	progress[20] = '\0';

	printf("Progress: [%s] %ld%%\r",progress,percent);
	fflush(stdout);
}

*/


/*******************************************************************/

void end_calculate_distances()
{
   // Save used defect coordinates and other info - only for first vertex point
    if(vertex_i==0){
        void write_defect_file(FILE *);
        FILE *f_def;
        std::string f_def_saveName = "defs_ico"+std::to_string(mesh_number)+"_"+std::to_string(distance_method)+".dat";
        f_def = fopen(f_def_saveName.c_str(),"w");
        write_defect_file(f_def);                       // 3.4.1
        fclose(f_def);
    }

	//int label_components();   // to CPP
	CPU_Time cpu_time;

	//FILE *f_ou;

	time(&t2);
  	get_time(t2-t1,&cpu_time);


	printf("\n\n");
	//printf("\tNumber of vertices %ld\n",num_of_meshpoint);
	printf("\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
  	printf("\n");

  	fprintf(f_ou,"\n\nDists = { ");
  	for(int d=0;d<num_of_defects;d++){              //print all
        fprintf(f_ou,"%.11f  ", vertex[vertex_i].distToDefect[d]);        //print all
    }       //print all
    fprintf(f_ou,"}\n");

    fprintf(f_ou,"\n\tCPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);

	 fprintf(f_ou,"\nDONE");

/*
	f_ou = fopen("ou.dat","w");
	fprintf(f_ou, "Mesh file: %s\n",f_name); //+++
	fprintf(f_ou, "Defects file: %s\n",fD_name); //+++
	fprintf(f_ou,"Number of vertices %ld\n",num_of_meshpoint);
	fprintf(f_ou,"Number of obtuse triangles %ld -> %.3f%%\n",num_of_obtuse,num_of_obtuse*100./num_of_triangles);  //+++
	fprintf(f_ou,"Number of domains %d\n\n",num_of_domains);
	fprintf(f_ou,"Integral of Eta-K is = %lg\n",total_disclinicity);
	fprintf(f_ou,"Total area %lg\n",total_area);
	fprintf(f_ou,"Total volume %lg\n",total_volume);    //++++
	fprintf(f_ou,"Radial variance asphericity is %.10g\n",asphericity);    //++++
	fprintf(f_ou,"Max dl = %.4g\n", dl_max);   //++++
	fprintf(f_ou,"Min dl = %.4g\n", dl_min);   //++++
	fprintf(f_ou,"Run time %g\n",run_time);
	fprintf(f_ou,"Num of steps %ld\n",num_of_iteration);
	fprintf(f_ou,"Time step %g\n",DT);
	fprintf(f_ou,"Total initial stress %lg\n",init_stress);
	fprintf(f_ou,"Range of local stress %lg\n",phi_range);
	fprintf(f_ou,"Total stress %lg\n",total_stress);
	fprintf(f_ou,"Total RHS^2 %lg\n",convergence);
	fprintf(f_ou,"CPU Time %ld:%ld:%ld:%ld\n",
	 cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
	fclose(f_ou);*/

}
/*******************************************************/
/*************************************************/
//  write_defect_file:
//       save defs.dat with defect information, including cost core radius
void write_defect_file(FILE *f_def)
{
    long d;
    fprintf(f_def,"$ Crystal coordination $\n%d\n",p_crystal);
    fprintf(f_def,"$ Number of defects $\n%ld\n",num_of_defects);
    fprintf(f_def,"$ Defect info $\nx\t\ty\t\tz\t\tqDens\t\tqOff\t\tcore_r\t\tClosest vertex");
    for (d=0; d<num_of_defects; d++){
        fprintf(f_def,"\n%.10f\t%.10f\t%.10f\t\t%.8f\t\t%G\t\t%.10f\t\t%ld",
            defect[d].x,
            defect[d].y,
            defect[d].z,
            defect[d].q,
            defect[d].qOff,
            defect[d].core_r,
            defect[d].closest_v);
    }
}

/*******************************************************************/

void write_triMshd(FILE *f_mshd)
{
	long i,d,t;
	fprintf(f_mshd,"N_vertices\n");
    fprintf(f_mshd,"%ld\n", num_of_meshpoint);
    fprintf(f_mshd,"N_defects\n");
    fprintf(f_mshd,"%ld\n", num_of_defects);
    fprintf(f_mshd,"i\tx\ty\tz\td1...N_defects\n");

    for(i=0; i<num_of_meshpoint; i++){
        fprintf(f_mshd,"%ld\t%.11f\t%.11f\t%.11f ", i, vertex[i].x,vertex[i].y,vertex[i].z);
    // print vertices info and idx (starting at 0)
        for(d=0; d<num_of_defects; d++){
        // read distances
            fprintf(f_mshd,"\t%.11f", vertex[i].distToDefect[d]);
        }
        fprintf(f_mshd,"\n");
    }
    fprintf(f_mshd,"N_triangles\n");
    fprintf(f_mshd,"%ld\n", num_of_triangles);
    for(t=0; t<num_of_triangles; t++){
        // print triangle vertices indx (starting at 0)
        fprintf(f_mshd,"%ld\t%ld\t%ld\n", triangle[t].v1,triangle[t].v2,triangle[t].v3);
    }

}

/*******************************************************************/

void write_triMshd_head(FILE *f_mshd)
{
	fprintf(f_mshd,"N_vertices\n");
    fprintf(f_mshd,"%ld\n", num_of_meshpoint);
    fprintf(f_mshd,"N_defects\n");
    fprintf(f_mshd,"%ld\n", num_of_defects);
    fprintf(f_mshd,"i\tx\ty\tz\td1...N_defects\n");
}
void write_triMshd_tail(FILE *f_mshd)
{
    long t;

    fprintf(f_mshd,"N_triangles\n");
    fprintf(f_mshd,"%ld\n", num_of_triangles);
    for(t=0; t<num_of_triangles; t++){
        // print triangle vertices indx (starting at 0)
        fprintf(f_mshd,"%ld\t%ld\t%ld\n", triangle[t].v1,triangle[t].v2,triangle[t].v3);
    }
}


/*******************************************************************/
/*          // ALREADY in common.cpp
void get_time(time_t t, CPU_Time *cpu_time)
{
  long sec_in_day, sec_in_ora, sec_in_min;

  sec_in_day = 86400;
  sec_in_ora = 3600;
  sec_in_min = 60;

  cpu_time->d = t/sec_in_day; t = t%sec_in_day;
  cpu_time->h = t/sec_in_ora; t = t%sec_in_ora;
  cpu_time->m = t/sec_in_min; t = t%sec_in_min;
  cpu_time->s = t;
}*/

/*******************************************************************/
