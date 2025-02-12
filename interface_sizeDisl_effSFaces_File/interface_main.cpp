/******************************************************************
**                                                               **
**                           INTERFACE_EFFSFACES                 **
**                                                               **
** Program: interface_effSFaces                                  **
** Version: 3.6 (disclination AND Faces with effective charge)   **
**                                                               **
** (From interface_effS)                                         **
** Last Change: Aug 24 2018                                      **
**                                                               **
*******************************************************************/

#include "common.h"

/*  Other files:
   Sources :
    - common.cpp
    - init_geometry.cpp
    - run_field_integration.cpp

   Headers:
    - common.h
   (and from geodesic_algorithm, in init_geometry.cpp):
    - geodesic_algorithm_base.h
    - geodesic_algorithm_dijkstra.h
    - geodesic_algorithm_dijkstra_alternative.h
    - geodesic_algortihm_exact.h
    - geodesic_algorithm_elements.h
    - geodesic_algorithm_graph_base.h
    - geodesic_algorithm_subdivision.h
    - geodesic_constants_and_simple_functions.h
    - geodesic_memory.h
    - geodesic_mesh.h
    - geodesic_mesh_elements.h
*/


/*******************************************************************/
/**************************        MAIN        *********************/
/*******************************************************************/

main()
{
    printf("\n");
    time_t t0, t1;
    time(&t0);                  // 1. Get general initial time

	f_ou = fopen("ou.dat","w"); // 2. Keep ou.dat open to write output throughout program

	init_geometry();            // 3. Get geometry input and calculate geo values

    time(&t1);                  // 4. Get integration initial time
	CPU_Time cpu_time;
	get_time(t1-t0,&cpu_time);  // 6.1
    printf("\n CPU Time on Geometry %ld:%ld:%ld:%ld\n",
     cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);
    fprintf(f_ou,"\nCPU_Geometry_Time %ld:%ld:%ld:%ld\n",
     cpu_time.d,cpu_time.h,
	 cpu_time.m,cpu_time.s);

    run_field_integration();    // 5. Integrate stress field equation
    fflush(f_ou);
	end(t0, t1);                // 6. Print final quantities     in common.cpp

    fclose(f_ou);

    printf("\nEND\n");
  	return EXIT_SUCCESS;
}

/*******************************************************************/
