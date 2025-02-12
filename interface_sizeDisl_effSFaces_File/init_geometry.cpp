///// geodesic_algorithm  ---
#include <iostream>
#include <fstream>
#include <string>

#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"
///// ---- geodesic_algorithm
#include "common.h"


char f_name[132];            // path+name of .msh file
char f_msh_name[132];            // name of .msh file
char fD_name[132];           // name of defect info file

int distance_method;        // flag for distance calculation algorithm or reading from file

Triangle triangle[MAX_SIZE];
Defect defect[MAX_SIZE];

long num_of_defects;
long num_of_triangles;

int displace = 1;           // Default:Locate each defect (x,y,z) at the mesh vertex closest to the initial coordinates.

double total_KInt;
double dropletRadius;        // get droplet constant volume through initial sphere radius
double lattice_constant;     // fix the lattice constant value given in inputParameters
double core_ratio_cst;       // assume all defects have the same relative size, read size from input
int p_crystal;
double burgerNumber;        // model parameter, strenght of dislocation vorticity field
double dislocExtend;        // model parameter, extend of dislocation vorticity field from defect

double R0_mesh;                // radius or input mesh to be rescaled to dropletRadius
double sizeScaling_factor;    // paramSize to meshSize factor to rescale coordinates
double calculated_radius;     // after rescaling to dropletRadius, check geo radius with calculated volume

void read_input_parameters(std::string, FILE *);             // 3.1
void import_mesh_and_distances(char *, FILE *);             // 3.2
void find_neighbors();                                      // 3.2.1
void import_mesh(char *);                                   // 3.3
void read_defects(char *);                                  // 3.4
void get_geometry();                                        // 3.5
double calculate_volume();                                  // 3.5.0    Added on 28-01-20
void add_defect_deformation();                              // 3.6
void write_ge_file(FILE *);                                 // 3.7

template<class Points, class Faces>
    void read_mesh_from_current_data(Points&, Faces&);      // 3.6.6.

/************************************************************************/

// 3. init_geometry
void init_geometry()
{
    printf("\n ****** Get Input and Geometry ******\n");
    read_input_parameters("inputParameters.txt", f_ou);     // 3.1
    fflush(f_ou);

    if (distance_method==4){
      printf("\n  Getting mesh from :  %s\n",f_name);
      import_mesh_and_distances(f_name, f_ou);                    // 3.2
      printf("\n  Distances have been calculated previously and now read from  %s too\n",f_msh_name);
      //read_defects(fD_name);         // To get crystal coordination   // now after geometry
	} else {
        printf("\n  Getting mesh from :  %s\n",f_name);
        import_mesh(f_name);                                // 3.3
        //printf("\n  Reading defects from :  %s\n",fD_name);
        //read_defects(fD_name);                              // 3.4            //now after geometry
    }
    fflush(f_ou);

    printf("\n  Calculating geometry \n");
    get_geometry();                                         // 3.5
                                                                                                                                                                                                                                                                                                                  fflush(f_ou);

    printf("\n  Reading defects from :  %s\n",fD_name);
    read_defects(fD_name);

    add_defect_deformation();                               // 3.6
    fflush(f_ou);

    // Save geometry values
	FILE *f_ge;
	f_ge = fopen("ge.dat","w");
	write_ge_file(f_ge);                                    // 3.7
    fclose(f_ge);
}

/************************************************************************/
/*******************************************************************/
/*************************************************************/

// 3.1. read_input_parameters:
//              get file information and run parameters from input_parameters.txt
void read_input_parameters(std::string f_in_name, FILE *f_out)
{
  // Read inputParameters.txt file
    FILE *fileStream;
    fileStream=fopen(f_in_name.c_str(),"r");
    fscanf(fileStream,"MeshFile - %s\n", f_name);
    fscanf(fileStream,"MeshName - %s\n", f_msh_name);
    fscanf(fileStream,"DefectFile - %s\n", fD_name);
    fscanf(fileStream,"%*s - %d\n", &distance_method);
    fscanf(fileStream,"%*[^\n] \n");        // skip this line
    fscanf(fileStream,"%*s . %*s - %lf . %lf\n", &dropletRadius, &R0_mesh);
    fscanf(fileStream,"%*s - %lf\n", &lattice_constant);
    fscanf(fileStream,"%*s - %lf\n", &core_ratio_cst);
    fscanf(fileStream,"%*s . %*s - %lf . %lf\n", &burgerNumber, &dislocExtend);
    fscanf(fileStream,"%*s . %*s - %lf . %lf\n", &init_stress, &phi_range);
    fscanf(fileStream,"%*s - %ld\n", &seed);
    //fscanf(fileStream,"%*s . %*s - %ld . %lG\n", &num_of_iteration, &run_time);       // change to read DT**
    fscanf(fileStream,"%*s . %*s - %ld . %lf\n", &num_of_iteration, &time_step);
    fscanf(fileStream,"%*s - %d\n", &integration_method);
    fscanf(fileStream,"%*[^\n] \n");        // skip this line
    fscanf(fileStream,"%*s _ %lf\n", &convergence_stop);

  // Fx global variables
    // sizeScaling_factor = dropletRadius/R0_mesh;    // Removed once volume calculation was included prior get_geometry() on 28-01-20
    burgerNumber = burgerNumber/lattice_constant;

  // Save read input in ou.dat
    fprintf(f_out, "<<< Input Parameters >>>\n");
    fprintf(f_out,"MeshFile  %s\n", f_name);
    fprintf(f_out,"MeshName  %s\n", f_msh_name);
    fprintf(f_out,"MeshRadius  %.3lf    --(with Evolver, this might be inconsistent, better read out measured)\n", R0_mesh);
    //fprintf(f_out,"ScalingFactor  %.7lG\n", sizeScaling_factor);
    fprintf(f_out,"DefectFile  %s\n", fD_name);
    fprintf(f_out,"DistanceMethod  %d\n", distance_method);
    fprintf(f_out,"InitialRadius[10mu]  %lG\n", dropletRadius);
    fprintf(f_out,"LatticeConstant[10mu]  %.7lG\n", lattice_constant);
    fprintf(f_out,"CoreRatio (CHECK on def.tmp) %.7lG\n", core_ratio_cst);
    fprintf(f_out,"CoreRadius (CHECK on def.dat)  %.7lG\n", dropletRadius*core_ratio_cst);
    fprintf(f_out,"DislocationBurger  %.7lG\n", burgerNumber);
    fprintf(f_out,"DislocationExtend  %.7lG\n", dislocExtend);
    fprintf(f_out,"InitStress    Range  %.4lG    %.4lG\n", init_stress, phi_range);
    fprintf(f_out,"Seed %ld\n", seed);
    //fprintf(f_out,"N_Steps      T_run   %ld     %lG\n", num_of_iteration, run_time);     // change to read DT**
    fprintf(f_out,"N_Steps      D_t   %ld     %lG\n", num_of_iteration, time_step);
    fprintf(f_out,"Integration  %d\n", integration_method);
    fprintf(f_out,"ConvergenceStop  %.3lG\n", convergence_stop);


}

/*******************************************************/

// 3.2. import_mesh_and_distances:
//              get mesh and distance to defect from the same file
void import_mesh_and_distances(char *f_name, FILE *f_out)
{
	long i, d, t, chi, dummy, num_of_edges=0;
	//std::string dummystring;
	char line[LINESIZE];


	FILE *f_in = fopen(f_name,"r");

	fgets(line,LINESIZE-1,f_in);    //"N_vertices"
    fscanf(f_in,"%ld \n", &num_of_meshpoint);
    if (num_of_meshpoint>MAX_SIZE){
		printf("Error: the number of vertices (%ld) exceeds MAX_SIZE (%d)\n",num_of_meshpoint,MAX_SIZE);
		exit(0);
	}
    fgets(line,LINESIZE-1,f_in);    //"N_defects"
    fscanf(f_in,"%ld \n", &num_of_defects);
    if (num_of_defects>MAX_SIZE){
		printf("Error: the number of defects (%ld) exceeds MAX_SIZE (%d)\n",num_of_defects,MAX_SIZE);
		exit(0);
	}
	fgets(line,LINESIZE-1,f_in);    // "i  x  y  z  d1...dNdefecs"

    //for(i=0; i<4; i++){
    for(i=0; i<num_of_meshpoint; i++){
        fscanf(f_in,"%ld %lg %lg %lg ", &dummy, &vertex[i].x,&vertex[i].y,&vertex[i].z);
    // read coordinates
        for(d=0; d<num_of_defects; d++){
        // read distances
            fscanf(f_in,"%lg ", &vertex[i].distToDefect[d]);

            //printf("\n\n\t THE distance for %ld to %ld is %lg \n\t The size scaling factor is %lg\n",i,d,vertex[i].distToDefect[d],sizeScaling_factor);
        /////////****** RESCALE distance read only after scalingFactor has been calculated
            //vertex[i].distToDefect[d] *= sizeScaling_factor;
        }
    }
    fscanf(f_in,"%*s\n");     //"N_triangles"
    fscanf(f_in,"%ld \n", &num_of_triangles);
    if (num_of_triangles>MAX_SIZE){
		printf("Error: the number of triangles (%ld) exceeds MAX_SIZE (%d)\n",num_of_triangles,MAX_SIZE);
		exit(0);
	}
    for(t=0; t<num_of_triangles; t++){
        // read triangle vestices
        fscanf(f_in,"%ld %ld %ld \n", &triangle[t].v1,&triangle[t].v2,&triangle[t].v3);
    }

	fclose(f_in);

	find_neighbors();                                       // 3.2.1

	for (i=0; i<num_of_meshpoint; i++){
		num_of_edges += vertex[i].num_of_neighbors;
	}

	if (!num_of_edges%2){
		printf("Error: bad triangulation, 2E = %ld\n",num_of_edges);
		exit(0);
	}
    // Eurler characteristic
	chi = num_of_meshpoint-num_of_edges/2+num_of_triangles;
	printf("\tEuler characteristic = %ld\n",chi);         //+++ verf
	//printf("meshpoint = %ld\n",num_of_meshpoint);         //+++ verf
	//printf("triangles = %ld\n",num_of_triangles);         //+++ verf
	//printf("edge= %ld\n",num_of_edges);         //+++ verf

    // change for topologies with more handles
	if (chi!=2){
		printf(" _____ Error: bad triangulation, chi = %ld______\n",chi);
		exit(0);
	}

  // Write mesh properties
    fprintf(f_out, "\n<<< Mesh Info >>>\n");
    fprintf(f_out,"N_vertices  %ld\n", num_of_meshpoint);
    fprintf(f_out,"N_triangles  %ld\n", num_of_triangles);
    fprintf(f_out,"EulerCharacteristic  %ld\n", chi);


}


/*******************************************************/
/*************************************************/
// 3.2.1 find_neighbors:
//              assign to each vertex element, which vertices index are neighbors

void find_neighbors()
{
	int is_neighbor(long,long);             // 3.2.1.1
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
/*************************************************/
/*****************************************/
// 3.2.1.1 is_neighbor:
//              function for find_neighbors() only

int is_neighbor(long i, long j)
{
	long k;

	for (k=0; k<vertex[i].num_of_neighbors; k++){
		if (j==vertex[i].neighbor[k]) return 1;
	}

	return 0;
}

/*****************************************/
/*************************************************/
/*******************************************************************/

// 3.3. import_mesh:
//              get mesh from the same file and calculate distances in this program

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

	find_neighbors();                   // 3.2.1

	for (i=0; i<num_of_meshpoint; i++){
		num_of_edges += vertex[i].num_of_neighbors;
	}

	if (!num_of_edges%2){
		printf("Error: bad triangulation, 2E = %ld\n",num_of_edges);
		exit(0);
	}
    // Eurler characteristic
	chi = num_of_meshpoint-num_of_edges/2+num_of_triangles;
	printf("\tEuler characteristic = %ld\n",chi);         //+++ verf

    // change for topologies with more handles
	if (chi!=2){
		printf("Error: bad triangulation, chi = %ld\n",chi);
		exit(0);
	}
	  // Write mesh properties
    fprintf(f_ou, "\n<<< Mesh Info >>>\n");
    fprintf(f_ou,"N_vertices  %ld\n", num_of_meshpoint);
    fprintf(f_ou,"N_triangles  %ld\n", num_of_triangles);
    fprintf(f_ou,"EulerCharacteristic  %ld\n", chi);


}

/*******************************************************/

// 3.5.0. calculate_volume: [[ADDED on 28-01-20]]
//                From current meshpoint positions, calculate volume of shape (for rescaling or checking purposes)

double calculate_volume()
{
    long t, i, j, k, n;
	double theta_i,theta_j,theta_k,cot_i,cot_j,cot_k;
	double x[3],y[3],z[3], dx, dy, base, height,side_i2,side_j2;
	double norm, area_t;

   //// Geometry info
	double areaTot = 0;
	double VolTot =0;

	// Initialize vector quantities
	for (i=0; i<num_of_meshpoint; i++){
        vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].nx = 0;
		vertex[i].ny = 0;
		vertex[i].nz = 0;
		vertex[i].area = 0;
	}

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
        areaTot += area_t;       //+++ CHK

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

        // Sum area contribution to vertices
        if (theta_i >= PI/2 ){
            vertex[i].area += area_t/2;
            vertex[j].area += area_t/4;
            vertex[k].area += area_t/4;
        } else if (theta_j >= PI/2){
            vertex[i].area += area_t/4;
            vertex[j].area += area_t/2;
            vertex[k].area += area_t/4;
        } else if (theta_k >= PI/2){
            vertex[i].area += area_t/4;
            vertex[j].area += area_t/4;
            vertex[k].area += area_t/2;
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

    for (i=0;i<num_of_meshpoint;i++){
      // Area
        total_area += vertex[i].area;

      // Mean Curvature
        vertex[i].hx /= vertex[i].area;
        vertex[i].hy /= vertex[i].area;
        vertex[i].hz /= vertex[i].area;

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
		VolTot += (vertex[i].x*vertex[i].nx+vertex[i].y*vertex[i].ny+vertex[i].z*vertex[i].nz)*vertex[i].area/3.;
    }

    double vol_sphere = 4*PI*pow(dropletRadius,3)/3.;
	printf("\tTotal volume = %g   (sphere of R=1 should be ~ %.4g)\n",VolTot, vol_sphere);
	//printf("\n");

    return VolTot;

}

// 3.5. get_geometry:
//              calculate geometry related quantities, including weights for integration

void get_geometry()
{
	long t, i, j, k, n, d;

	double theta_i,theta_j,theta_k,cot_i,cot_j,cot_k;
	double x[3],y[3],z[3], dx, dy, base, height,side_i2,side_j2;
	double norm, area_t;

   //// Control obtuse angles
    long num_of_obtuse = 0;
	long obtuse_t[MAX_SIZE];
	float obtuse_tolerance = OBTUSE_tolerance;    //+++ admit only this % of obtuse angles in msh
	float obtuse_percent;

   //// Geometry info
	total_area = 0;
	total_volume =0;
	double asphericity = 0;
	double max_dZ = 0.0;

    //// Control
	double total_K = 0.;
	total_KInt = 0.;
    dl_min = 1E5;                          // extern variable
	double dl_max = 0.0;

  //// Asphericity
	double mean_radius = 0;
	double radial_deviation = 0;

  //// Rescaling
    printf("\n\tCHECKING INITIAL VOLUME READ FROM INPUT MESH\n");
	double tmpVol = calculate_volume();
	calculated_radius = pow(tmpVol*3./(4*PI),1./3);
	printf("\tcalculated Radius %.7lf",calculated_radius);      //+++++CHK
    printf("\n\tsize difference(set-calculated) = %.4lG%%  (should be ~0)\n",(dropletRadius-calculated_radius)*100./dropletRadius);
	sizeScaling_factor = dropletRadius/calculated_radius;


    // Rescale meshpoint coordinates to given dropletRadius
    for (i=0; i<num_of_meshpoint; i++){
        vertex[i].x *= sizeScaling_factor;
        vertex[i].y *= sizeScaling_factor;
        vertex[i].z *= sizeScaling_factor;
    }
   // Rescale distances read from tri_msh generated elsewhere; for the other methods, defect positions are instead rescaled
    if(distance_method==4){
        for (i=0; i<num_of_meshpoint; i++){
            for(d=0;d<num_of_defects;d++){
                vertex[i].distToDefect[d] *= sizeScaling_factor;
            }
        }
    }

    printf("\n\n\tFinding out volume after....\n");
	tmpVol = calculate_volume();
	calculated_radius = pow(tmpVol*3./(4*PI),1./3);
	printf("\tcalculated Radius %.7lf",calculated_radius);      //+++++CHK
    printf("\n\tsize difference(set-calculated) = %.4lG%%  (should be ~0)\n\n",(dropletRadius-calculated_radius)*100./dropletRadius);

  /////  Initialize
	// Initialize vector quantities
	for (i=0; i<num_of_meshpoint; i++){
        vertex[i].hx = 0;
		vertex[i].hy = 0;
		vertex[i].hz = 0;
		vertex[i].h2 = 0;
		vertex[i].kg = 2*PI;
        vertex[i].nx = 0;
		vertex[i].ny = 0;
		vertex[i].nz = 0;
		vertex[i].area = 0;
		vertex[i].S = 0;
        vertex[i].vorticity = 0;
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

        // Substract theta angles to corresponding Gaussian kg
        vertex[i].kg -= theta_i;
        vertex[j].kg -= theta_j;
        vertex[k].kg -= theta_k;

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


        // Sum weight contribution
        // vertex i, find position of vertices j,k
        for (n=0; n<vertex[i].num_of_neighbors; n++){
            if (vertex[i].neighbor[n]==j){vertex[i].weight[n] += cot_k/2;}
            if (vertex[i].neighbor[n]==k){vertex[i].weight[n] += cot_j/2;}
		}
		// vertex j
		for (n=0; n<vertex[j].num_of_neighbors; n++){
            if (vertex[j].neighbor[n]==i){vertex[j].weight[n] += cot_k/2;}
            if (vertex[j].neighbor[n]==k){vertex[j].weight[n] += cot_i/2;}
		}
		// vertex k
		for (n=0; n<vertex[k].num_of_neighbors; n++){
            if (vertex[k].neighbor[n]==j){vertex[k].weight[n] += cot_i/2;}
            if (vertex[k].neighbor[n]==i){vertex[k].weight[n] += cot_j/2;}
		}
	}

	//printf("\tTotal surface area from triangles area %g\n",total_area); //++++ chk
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

        // Weights
        for (j=0; j<vertex[i].num_of_neighbors; j++){
			vertex[i].weight[j] /= vertex[i].area;
		}

        // Gaussian Curvature
        total_K += vertex[i].kg;
        vertex[i].kg /= vertex[i].area;
        total_KInt += vertex[i].kg*vertex[i].area;          //+++ BCK

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

    // Number of surfactant lattice points approximating all crystal cells as hexagonal to an avg lattice constant by input
    double area_crystal_cell = sqrt(3.)*lattice_constant*lattice_constant/2.;
    double number_latticePoints  = long(total_area/area_crystal_cell)*1.;

    // Geometry Z height
    max_dZ  = max_z-min_z;
    //double geo_halfHeight = max_dZ/2;
    calculated_radius = pow(total_volume*3./(4*PI),1./3);
    double vol_sphere = 4*PI*pow(dropletRadius,3)/3.;
    double area_sphere = 4*PI*pow(dropletRadius,2);

	//printf("\tTotal volume with (x,y,z)/3 = %.10g   (sphere should be ~ %.6g)\n",total_volume, vol_sphere);    //+++ CHK
	double tst_volume = 0;
    for(i=0; i<num_of_meshpoint; i++){
		tst_volume += vertex[i].x*vertex[i].nx*vertex[i].area;
        }//printf("\tTotal volume with (x,0,0) = %.10g   (sphere should be ~ %.6g)\n",total_volume, vol_sphere);
    tst_volume = 0;
    for(i=0; i<num_of_meshpoint; i++){
		tst_volume += vertex[i].y*vertex[i].ny*vertex[i].area;
        }//printf("\tTotal volume with (0,y,0) = %.10g   (sphere should be ~ %.6g)\n",total_volume, vol_sphere);
    tst_volume = 0;
    for(i=0; i<num_of_meshpoint; i++){
        tst_volume += vertex[i].z*vertex[i].nz*vertex[i].area;
        }//printf("\tTotal volume with (0,0,z) = %.10g   (sphere should be ~ %.6g)\n",total_volume, vol_sphere);
        ////////+++++ */


    obtuse_percent = 100.0*num_of_obtuse/num_of_triangles;

	//printf("\n");
	printf("\tVertices %ld\n",num_of_meshpoint);
	printf("\tTriangles %ld\n",num_of_triangles);
	printf("\tObtuse triangles %ld  -> %.3f%%\n",num_of_obtuse,obtuse_percent);
	printf("\tTotal surface area = %g  (sphere should be ~ %.5g)\n",total_area,area_sphere);
	printf("\tTotal volume = %g   (sphere should be ~ %.4g)\n",total_volume, vol_sphere);
	printf("\tApprox number crystal cells = %.1le  \n",number_latticePoints);
	printf("\tAsphericity by radial variance = %.7G   (sphere should be 0)\n",asphericity);
	//printf("\tDisclinicity over surface = %.9lG  (should be 0)\n",total_disclinicity);  //+++     //RM
	printf("\tMax dl = %.4G, Min dl = %.4G\n", dl_max, dl_min);   //++++
	printf("\tDt < %.3lG  (if cut .45)", DT_tolerance*dl_min*dl_min); //+++
	printf("\n");

	// Check is percentage of obstuse angles in msh less than tolerance
	if (obtuse_percent > obtuse_tolerance)
	{
        printf("Too many obtuse triangles\nPROBLEM\n");
        exit(1);
	}

	// CHECK RESCALING
    // check read vs calculated:
    //printf("\n\t  set %.7lf",dropletRadius);      //+++++CHK
    printf("\n\t  calculated Radius %.7lf",calculated_radius);      //+++++CHK
    printf("\n\t size difference(set-calculated) = %.4lG%%  (should be ~0)\n",(dropletRadius-calculated_radius)*100./dropletRadius);

  // Write geometry properties
    fprintf(f_ou, "\n<<< Geometry Info >>>\n");
    fprintf(f_ou,"N_obtuse  %ld  (%.3f%%)\n", num_of_obtuse, obtuse_percent);
    fprintf(f_ou,"Scaling Factor  %.6lG\n", sizeScaling_factor);
    fprintf(f_ou,"Area  %.6lG\n", total_area);
    fprintf(f_ou,"Volume  %.6lG\n", total_volume);
    fprintf(f_ou,"CalculatedRadius  %.5lG\n", calculated_radius);
    fprintf(f_ou,"N_Molecules = %.1le\n",number_latticePoints);
    fprintf(f_ou,"Asphericity  %.6lG\n", asphericity);
    fprintf(f_ou,"Height  %.6lG\n", max_z);
    fprintf(f_ou,"Max/Min_edge_length  %.4lG  %.4lG\n", dl_max, dl_min);
}

/*******************************************************/

// 3.4. read_defects:
//              read defect information (but take core radius from inputParameters)

void read_defects(char *fD_name)
{
    double charge_screening(long);                              // 3.4.2
    long d;

    FILE *f_in = fopen(fD_name,"r");
    char line[LINESIZE+1];


    // Read crystal ideal coordination number
    fgets(line,LINESIZE-1,f_in);    // $Crystal coordination $
    fscanf(f_in,"%d\n",&p_crystal);
    fprintf(f_ou, "CrystalCoordination  %d\n", p_crystal);
    printf("\tCrystalCoordination  %d\n", p_crystal);

    // Read number of defects
    //fgets(line,LINESIZE-1,f_in);    // newline
    fgets(line,LINESIZE-1,f_in);    // $ Number of defects $
    //fscanf(f_in,"%ld",&num_of_defects);
    if(distance_method==4){fgets(line,LINESIZE-1,f_in);} else {fscanf(f_in,"%ld\n",&num_of_defects);}
    if (num_of_defects>MAX_SIZE){
		printf("Error: the number of defects (%ld) exceeds MAX_SIZE (%d)\n",num_of_defects,MAX_SIZE);
		exit(0);
	}

     // Read defect information
    //fgets(line,LINESIZE-1,f_in);    // newline
    fgets(line,LINESIZE-1,f_in);    // $ Defect info $
    fgets(line,LINESIZE-1,f_in);    // x   y  z  q  qOff  core_r

    //    double core_size_cst = core_ratio_cst*dropletRadius;
    double dummy;
    double core_ratio, core_size;
    for (d =0; d<num_of_defects; d++){
        fscanf(f_in,"%lg %lg %lg %lg %g %lg",
        &defect[d].x,
        &defect[d].y,
        &defect[d].z,
        &defect[d].q,
        &defect[d].qOff,        // -1 for dislocation cloud at vertices, or 1>qOff>0 for faces
        &core_ratio);                              // All defects SAME core radius
        //&defect[d].core_r);                     // DIFFERENT radii for each defect
        core_size = core_ratio*dropletRadius;
        defect[d].core_r = core_size;       //  Consider all defects to have same core radius
    }


    fclose(f_in);

    fprintf(f_ou,"N_defects  %ld\n", num_of_defects);
    printf("\tN_defects  %ld\n", num_of_defects);
    printf("\tCore Radius  = %.3lG\n", core_size);

    // write configuration used in ou.dat
    fprintf(f_ou,"PairConfigList  [ ");
    printf("\tPairConfigList = [ ");
    for (d=0; d<num_of_defects; d++){
        fprintf(f_ou, "%G, ",defect[d].qOff);
        printf("%G, ",defect[d].qOff);
    }
    fprintf(f_ou,"]\n");
    printf("]\n");

    // Defect disclinicity
    for (d=0; d<num_of_defects; d++){
        defect[d].q *= 2*PI/p_crystal;       //*** GOBACK
        //printf("for defect %ld, charge density is %lg\n",d,defect[d].q);    //+++ chk
    }

    // Effective charge of screened disclination
    for (d=0; d<num_of_defects; d++){
        defect[d].q += defect[d].qOff*charge_screening(d);
        //printf("\n");printf("for defect %ld, charge density is %lg\n",d,defect[d].q);    //+++ chk
    }

  // Rescale defect coordinates to given dropletRadius
    for(d=0;d<num_of_defects;d++){
        defect[d].x *= sizeScaling_factor;
        defect[d].y *= sizeScaling_factor;
        defect[d].z *= sizeScaling_factor;
    }

     // Locate each defect (x,y,z) at the mesh vertex closest to the initial coordinates.
    // False if defect centers at exact tips of inscribed icosahedron
    //if(distance_method==1) displace =0;           ///+++++++++++CHECK VOLVER

    double r[3],r2;
    for(d=0;d<num_of_defects;d++){
        double min_dist = 1E6;
        long i, i_min = 0;
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
    } else{printf("\tDefects NOT displaced\n");}


  // Save used defect coordinates and other info
    void write_defect_file(FILE *);
    FILE *f_def;
    f_def = fopen("defs.dat","w");
    write_defect_file(f_def);                       // 3.4.1
    fclose(f_def);

}

/*******************************************************/
/*************************************************/
// 3.4.1. write_defect_file:
//                      save defs.dat with defect information, including cost core radius
void write_defect_file(FILE *f_def)
{
    long d;
    fprintf(f_def,"$ Crystal coordination $\n%d\n",p_crystal);
    fprintf(f_def,"$ Number of defects $\n%ld\n",num_of_defects);
    fprintf(f_def,"$ Defect info $\nx\t\ty\t\tz\t\tqDens\t\tqOff\t\tcore_radius\t\tClosest vertex");
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

/*************************************************/
// 3.4.2 charge_screening:
//              function to calculate the effective charge of a screened disclination
double charge_screening(long d)
{
    double screening = lattice_constant/calculated_radius;     // ****GOBACK
    double offset = burgerNumber*screening;
    return offset;

}
/*************************************************/
/*******************************************************/
/*******************************************************************/

// 3.6. add_defect_deformation:
//              read defect core radius and (depending on the distance method), calculate distances
void add_defect_deformation()
{
  // Internat functions forward declarations
    double great_circle_distance(long, long );          // 3.6.1.
    double delta_approx(long, double);                  // 3.6.2.
    double vorticity_approx(long ,double);               // 3.6.2b
    void write_triMshd(FILE *);                         // 3.6.4
    void write_triMshd_head(FILE *);                    // 3.6.5
    void write_triMshd_tail(FILE *);                    // 3.6.7

    long i;
    double r[3],distance_to_defect, r2;
    double discl_chk = 0;
    double defEff = 0.;
    double dislEff = 0.;

    double total_disclinicity = 0;
    double dislocation_sum = 0;

/*
    // Locate each defect (x,y,z) at the mesh vertex closest to the initial coordinates.
    // False if defect centers at exact tips of inscribed icosahedron
    int displace = 1;
    if(distance_method==1) displace =0;

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
    } else{printf("\tDefects NOT displaced\n");}

*/                  //RM


  // Write .mshd file with distances and triangle information
    std::string f_mshd_name;
    f_mshd_name = "tri_"+std::string (f_msh_name)+"d";
    FILE *f_mshd = fopen(f_mshd_name.c_str(), "w");


    // Variables needed for geodesic_algorithm
    std::vector<double> points;     // keep types because it is read by mesh.initialize_mesh_data
    std::vector<unsigned> faces;
    std::vector<geodesic::SurfacePoint> path;
    geodesic::Mesh mesh;
    unsigned source_vertex_index, target_vertex_index;
    // Choose how to calculate distance on surface
        // 1: great-circle distance on a sphere, 2:dijkstra or 3: exact geodesic (mesh) algortihm.
        // 4: read distances from file, 5: euclidean
    switch(distance_method){
    case 1: // circle
    {
        printf("\n  ** Using great-circle distance**\n");

      //  Include progress bar on ou.dat to check remotely
        int progress_bar_count = 0;
        fprintf(f_ou, "\n ---- Distance Progress: [ 10 .");
        fflush(f_ou);

        printf("\n  \n");

        for (long i=0; i<num_of_meshpoint; i++){
          // write to tri*mshd
            for(long d=0;d<num_of_defects;d++){
                distance_to_defect = great_circle_distance(i,d);        // 3.6.1
                r2 = distance_to_defect*distance_to_defect;
                defEff = defect[d].q*delta_approx(d,r2);                // 3.6.2
                vertex[i].S += defEff;
                dislEff = burgerNumber*vorticity_approx(d,distance_to_defect);          // 3.6.2b
                //vertex[i].vorticity +=dislEff;
                vertex[i].vorticity +=defEff;
                vertex[i].distToDefect[d] = distance_to_defect;   //LLLL
            }
            progress_bar(i, num_of_meshpoint);                          // 3.6.3  IN: common.cpp
            total_disclinicity += vertex[i].S*vertex[i].area;
            dislocation_sum += vertex[i].vorticity*vertex[i].area;
            progress_bar_count = progress_bar_on_ou(i, num_of_meshpoint, progress_bar_count);            // 3.6.3b (too)  IN: common.cpp
        }
        write_triMshd(f_mshd);
        fprintf(f_ou, " 0! ]  ;) \n "); fflush(f_ou);                                        // 3.6.4
        break;
    }

    case 2:  // dijkstra
    {
        printf("\n  ** Using Dijkstra algorithm to calculate geodesic distances**\n");

        // SAVE tri*mshd file Line-by-line
        write_triMshd_head(f_mshd);                                     // 3.6.5

      //  Include progress bar on ou.dat to check remotely
        int progress_bar_count = 0;
        fprintf(f_ou, "\n ---- Distance Progress: [ 10 .");
        fflush(f_ou);                             // 3.6.5

        read_mesh_from_current_data(points, faces);                     // 3.6.6
        // + create internal mesh data
        mesh.initialize_mesh_data(points, faces);
        geodesic::GeodesicAlgorithmDijkstra dijkstra_algorithm(&mesh);

        for (long i=0; i<num_of_meshpoint; i++){
          // write to tri*mshd
            fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", i, vertex[i].x,vertex[i].y,vertex[i].z);
            for(long d=0;d<num_of_defects;d++){
                source_vertex_index = defect[d].closest_v;
                target_vertex_index = i;
                geodesic::SurfacePoint source (&mesh.vertices()[source_vertex_index]);
                geodesic::SurfacePoint target (&mesh.vertices()[target_vertex_index]);
                dijkstra_algorithm.geodesic(source, target, path);
                distance_to_defect = geodesic::length(path);

                r2 = distance_to_defect*distance_to_defect;
                defEff = defect[d].q*delta_approx(d,r2);
                vertex[i].S += defEff;
                dislEff = burgerNumber*vorticity_approx(d,distance_to_defect);          // 3.6.2b
                //vertex[i].vorticity +=dislEff;
                vertex[i].vorticity +=defEff;
                vertex[i].distToDefect[d] = distance_to_defect;   //LLLL
                fprintf(f_mshd,"\t%.11f", distance_to_defect);
            }
            progress_bar(i, num_of_meshpoint);
            fprintf(f_mshd,"\n");
            total_disclinicity += vertex[i].S*vertex[i].area;
            dislocation_sum += vertex[i].vorticity*vertex[i].area;
            progress_bar_count = progress_bar_on_ou(i, num_of_meshpoint, progress_bar_count);            // 3.6.3b (too)  IN: common.cpp
        }
        write_triMshd_tail(f_mshd);
        fprintf(f_ou, " 0! ]  ;) \n "); fflush(f_ou);                                 // 3.6.7
        break;
    }

    case 3:  // exact geodesic_algortihm
    {
        printf("\n  ** Using Exact algorithm to calculate geodesic distances**\n");

        // SAVE tri*mshd file Line-by-line
        write_triMshd_head(f_mshd);                                         // 3.6.5

      //  Include progress bar on ou.dat to check remotely
        int progress_bar_count = 0;
        fprintf(f_ou, "\n ---- Distance Progress: [ 10 .");
        fflush(f_ou);

        read_mesh_from_current_data(points, faces);
        // + create internal mesh data
        mesh.initialize_mesh_data(points, faces);
        geodesic::GeodesicAlgorithmExact exact_algorithm(&mesh);

        for (long i=0; i<num_of_meshpoint; i++){
          // write to tri*mshd
            fprintf(f_mshd, "%ld\t%.11f\t%.11f\t%.11f ", i, vertex[i].x,vertex[i].y,vertex[i].z);
            for(long d=0;d<num_of_defects;d++){
                source_vertex_index = defect[d].closest_v;
                target_vertex_index = i;
                geodesic::SurfacePoint source (&mesh.vertices()[source_vertex_index]);
                geodesic::SurfacePoint target (&mesh.vertices()[target_vertex_index]);
                exact_algorithm.geodesic(source, target, path);
                distance_to_defect = geodesic::length(path);

                r2 = distance_to_defect*distance_to_defect;
                defEff = defect[d].q*delta_approx(d,r2);
                vertex[i].S += defEff;
                dislEff = burgerNumber*vorticity_approx(d,distance_to_defect);          // 3.6.2b
                //vertex[i].vorticity +=dislEff;
                vertex[i].vorticity +=defEff;
                vertex[i].distToDefect[d] = distance_to_defect;   //LLLL
            }
            progress_bar(i, num_of_meshpoint);
            fprintf(f_mshd,"\n");
            total_disclinicity += vertex[i].S*vertex[i].area;
            dislocation_sum += vertex[i].vorticity*vertex[i].area;
            progress_bar_count = progress_bar_on_ou(i, num_of_meshpoint, progress_bar_count);            // 3.6.3b (too)  IN: common.cpp
        }
        write_triMshd_tail(f_mshd);
        fprintf(f_ou, " 0! ]  ;) \n "); fflush(f_ou);
        break;
    }
    case 4:
    {
        printf("\n   ** Distances have been read from .mshd file **\n");
        for (long i=0; i<num_of_meshpoint; i++){
            for(long d=0;d<num_of_defects;d++){
                distance_to_defect = vertex[i].distToDefect[d];
              // DISCLINATION
                defEff = defect[d].q*delta_approx(d,distance_to_defect*distance_to_defect);
                vertex[i].S += defEff;
                // CHECK S per vertex:
                //vertex[i].distToDefect[d] = defEff;        //++++++ CHK

              // DISLOCATION -  explicit
                dislEff = burgerNumber*vorticity_approx(d,distance_to_defect);          // 3.6.2b
                //vertex[i].vorticity +=dislEff;
                vertex[i].vorticity +=defEff;
            }

            total_disclinicity += vertex[i].S*vertex[i].area;
            dislocation_sum += vertex[i].vorticity*vertex[i].area;
        }
        break;
    }
    case 5:
    {
        printf("\n   ** Using euclidean distances**\n");
        for (long i=0; i<num_of_meshpoint; i++){
            for(long d=0;d<num_of_defects;d++){
                r[0] = vertex[i].x-defect[d].x;
                r[1] = vertex[i].y-defect[d].y;
                r[2] = vertex[i].z-defect[d].z;
                r2 = r[0]*r[0]+r[1]*r[1]+r[2]*r[2];
                distance_to_defect = sqrt(r2);
                defEff = defect[d].q*delta_approx(d,r2);
                vertex[i].S += defEff;
                vertex[i].distToDefect[d] = distance_to_defect;
                dislEff = burgerNumber*vorticity_approx(d,distance_to_defect);          // 3.6.2b
                //vertex[i].vorticity +=dislEff;
                vertex[i].vorticity +=defEff;
            }
            total_disclinicity += vertex[i].S*vertex[i].area;
            dislocation_sum += vertex[i].vorticity*vertex[i].area;
        }
        write_triMshd(f_mshd);                                              // 3.6.4
        break;
    }

    }
    fclose(f_mshd);
    if(distance_method==4){remove(f_mshd_name.c_str());}

    // CHECK DIstances  +++DistError
    /*printf("\n\n\n");
    printf("Distance from d%d to i %ld exp =%.8lG\tdouble%lf  (dN%ld)",0,18686,vertex[18686].distToDefect[0],vertex[18686].distToDefect[0],defect[0].closest_v);
    printf("\nDistance from d%d to i %ld exp =%.8lG\tdouble%lf  (dN%ld)",6,18690,vertex[18690].distToDefect[6],vertex[18690].distToDefect[6],defect[6].closest_v);
    printf("\nDistance from d%d to i %ld exp =%.8lG\tdouble%lf  (dN%ld)",28,5457,vertex[5457].distToDefect[28],vertex[5457].distToDefect[28],defect[28].closest_v);
    printf("\n\n\n");*/


    // The integral of S should be 4Pi, so normalize accordingly
    printf("\n\tThe integrated defect density before correction for the given delta approximation is:\n");
    printf("\trho/(2PI) = %lG  (should be = 2)\n",total_disclinicity/2./PI);
    printf("\t -> Correcting by normalizing to K/rho factor\n");

  // Correct total disclination charge equal to int K, AND int vorticity equal zero
    double dislocation_check_before = dislocation_sum;       //+++TMP
    double dislocation_check_after = 0.;       //+++TMP

    for (i=0; i<num_of_meshpoint; i++){
        vertex[i].S = vertex[i].S*total_KInt/total_disclinicity;        //+++ BCK total_K   //LLLLL comment
        vertex[i].vorticity = vertex[i].vorticity-dislocation_sum/total_area;

        // chck sum now
        discl_chk += vertex[i].S*vertex[i].area;        // ++++ chk
        dislocation_check_after += vertex[i].vorticity*vertex[i].area;      //+++TMP
    }
    printf("\tCHECK dislocation_sum before = %lg\n",dislocation_check_before);  //+++TMP
    printf("\tCHECK dislocation_sum after = %lg\n",dislocation_check_after);       //+++TMP

    /*
    if (discl_chk!=2){printf("\t Integral of defect density rho/{2PI} is = %lg, the diff with 2 is = %lg\n",discl_chk/2/PI,2.-discl_chk/2/PI);}
    else {printf("\n\tnumerical EUREKA\n");} */     // ++++ chk

    // Add gaussian curvature screening and corrected vorticity
    for (i=0; i<num_of_meshpoint; i++){
        vertex[i].S -= vertex[i].kg;
        //vertex[i].S += vertex[i].vorticity;     //+++** no explicit dislocation density, only screened charge
    }
    total_disclinicity = 0;
    dislocation_sum = 0;
    for (i=0; i<num_of_meshpoint; i++){
        total_disclinicity += vertex[i].S*vertex[i].area;
        dislocation_sum += vertex[i].vorticity*vertex[i].area;
    }
    printf("\tIntegral of rho-K now = %lg\n",total_disclinicity);
    printf("\tIntegral of omega now = %lg\n",dislocation_sum);
    fprintf(f_ou,"DisclinicitySum  %.4lG\t\t(should be ~0)\n", total_disclinicity);
    fprintf(f_ou,"DislocationSum  %.4lG\t\t(should be ~0)\n", dislocation_sum);

}

/*******************************************************************/
/*************************************************/
// 3.6.1. great_circle_distance:
//              calculate geodesic distance on a sphere as the arc of the great circle connecting the two points
double great_circle_distance(long i, long d){

    double r[3], acosArg, l_arc;

  // If defects are at exact meshpoint, return l_arc = 0.
    if (displace) {
        if (i==defect[d].closest_v){
            //printf("\n\n\nJACKPOT\n\n\n");      ++Chk DistError
            return 0.;
        }
    }
  //Distance from vertex to defect core (still, relatively small core radius, shortest dist on surface)
    r[0] = vertex[i].x*defect[d].x;
    r[1] = vertex[i].y*defect[d].y;
    r[2] = vertex[i].z*defect[d].z;
    acosArg = (r[0]+r[1]+r[2])/pow(calculated_radius,2.);

    // CHECK        ++Chk DistError
    /*if (i==18690 && d==6) {
        printf("\n\n\n\n(x,y,z)i-(x,y,z)d %.23lG %.23lG %.23lG (should be 0)",vertex[i].x,vertex[i].y,vertex[i].z);
        printf("\n x,y,z)i-(x,y,z)d %.23lG %.23lG %.23lG (should be 0)",defect[d].x,defect[d].y,defect[d].z);
        printf("GOTCHA\n\t %lf (should be 1)\n\n\n",acosArg);
    }
    if (i==19126 && d==2){
     printf("\n\n\n\n(x,y,z)i-(x,y,z)d %.23lG %.23lG %.23lG (should be 0)",vertex[i].x-defect[d].x,vertex[i].y-defect[d].y,vertex[i].z-defect[d].z);
        printf("GOTCHA   OK\n\t %lf (should be 1)\n\n\n",acosArg);}*/

  // correct numerical precision errors
    if (acosArg <-1) acosArg = -1.;
    else if (acosArg>1) acosArg = 1.;

    l_arc = acos(acosArg)*calculated_radius;
    /*if (i==18690 && d==6) printf("\n\n\n\nDIST\n\t %lf (should be 0)\n\n\n",l_arc);
    if (i==19126 && d==2) printf("\n\n\n\nDIST  OK\n\t %lf (should be 0)\n\n\n",l_arc);*/
    return l_arc;
}

/*******************************************************************/
/*************************************************/
// 3.6.2. delta_approx:
//              function to approximate DISCLINATION density
double delta_approx(long d, double dist2)
{
    //double delta = 1/pow(sqrt(2*PI)*defect[d].core_r,3);
    double delta = 1.0;
    delta = delta*exp(-dist2/pow(defect[d].core_r,2)/2);
    return delta;
}

/*
/*******************************************************************/
/*/*************************************************/
// 3.6.2b vorticity_approx:
//              function to approximate DISLOCATION density through burger vector field vorticity
double vorticity_approx(long d,double r)
{
    //assume burger same for all d, otherwise burger = defect[d].burger
    double c_extend = dislocExtend*calculated_radius;
    // ensure a/R^3 scaling of the dislocation charge density
    double dislocationScaling = lattice_constant/pow(calculated_radius,3);
    double omega = -dislocationScaling*exp(-r/c_extend);
    //return omega;
    return 0.0;
}
/*******************************************************************/
/*************************************************/
// 3.6.4. write_trisMshd:
//              for fast distance calculation method, write the tri_mshd file at once
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
/*************************************************/
// 3.6.5. write_trisMshd_head:
// 3.6.7. write_triMshd_tail:
//              for slow distance calculation method, write the tri_mshd file on the fly
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
/*************************************************/
// 3.6.6. read_mesh_from_current_data:
//              create Mesh object for geodesic_algorithm
template<class Points, class Faces>
void read_mesh_from_current_data(Points& points, Faces& faces)
{
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
/*************************************************/
/*******************************************************************/
// 3.7. write_ge_file:
//              write ge.dat file containing all geo info by vertex i
void write_ge_file(FILE *f_ou)
{
    long i,j;

	for (i=0; i<num_of_meshpoint; i++){
		fprintf(f_ou,"%ld\t%.10f\t%.10f\t%.10f\t%.10f\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10lg\t%.10f\t%.10f\t%.10f\t",
			i,
			vertex[i].x,
			vertex[i].y,
			vertex[i].z,
			vertex[i].area,
			vertex[i].h2,
			vertex[i].kg,
            vertex[i].S,
            vertex[i].vorticity,
            vertex[i].hx,
            vertex[i].hy,
            vertex[i].hz,
            vertex[i].nx,
            vertex[i].ny,
            vertex[i].nz);

        //LLLL

        for (j=0; j<num_of_defects; j++){
            fprintf(f_ou,"%.8lg\t",vertex[i].distToDefect[j]);
        }           //LLLLLL*/

		for (j=0; j<vertex[i].num_of_neighbors; j++){
			fprintf(f_ou,"%.10f\t",vertex[i].weight[j]);
		}

		fprintf(f_ou,"\n");
	}
}
/*******************************************************************/
/************************************************************************/
