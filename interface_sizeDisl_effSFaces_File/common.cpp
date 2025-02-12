#include "common.h"

// *** Mesh variables used throughout the code **
long num_of_meshpoint;
Vertex vertex[MAX_SIZE];
double total_area;
double total_volume;
double dl_min;              // minimal edge length in the mesh -> DT_max calculation

// *** File variables used throughout the code **
FILE *f_ou;

// *** Integration parameters read from input **
double init_stress;                 // initial total pre-stress
double phi_range;                   // set the initial stress at each point randomly within this range
long seed;                          // for random fx to set initial stress value
long num_of_iteration;
//double run_time;                    // integration time step set accordingly        // change to read DT**
double time_step;                   // integration time step
int integration_method;             // 1. Euler   2. Runge-Kutta (2 steps)  3. Runge-Kutta (4 steps) -> [bug]
double convergence_stop;            // stop integration when num_iterations completed or reaching this error

// *** integration variables used throughout the code **
double convergence;



/*******************************************************************/
/*******************************************************************/

// **** Common functions   **

/************************************************************************/
// 3.6.3. progress_bar:
//          print on-screen percentage bar of a process at t from total
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
/****************************************************/
// 3.6.4. progress_bar_on_ou:
//          print on ou.dat percentage bar of a process at t from total
int progress_bar_on_ou(long t, long total, int prog_count)
{
	float percent = 100*t/total;

	int current_count = percent/10;
	if (current_count > prog_count)
	{
        fprintf(f_ou," %d .",10-current_count);
        fflush(f_ou);
    }
    return current_count;

}
/************************************************************************/
/*******************************************************************/
// 6. end:
//      Get end times and print (on screen and in ou.dat) final outputs
void end(time_t t0, time_t t1)
{
    CPU_Time cpu_time_total,cpu_time_geo,cpu_time_run;
    time_t t2;
	time(&t2);
  	get_time(t2-t1,&cpu_time_run);                  // 6.1
    get_time(t2-t0, &cpu_time_total);
    get_time(t1-t0, &cpu_time_geo);

	export_conf(-1,-1);         //save last.dat file        // 5.2.2 in run_field_integration

    double total_stretch=0;
    for(long i=0; i<num_of_meshpoint; i++){
        total_stretch += vertex[i].phi*vertex[i].phi*vertex[i].area;
    }

  // Write final values of interest in ou.dat
    fprintf(f_ou, "\n<<< Final Variables >>>\n");
    fprintf(f_ou,"TotalStretch  %.8lG\n", total_stretch);
    fprintf(f_ou,"RHS^2  %.4lG\n", convergence);
    /*fprintf(f_ou,"CPU_Geometry_Time %ld:%ld:%ld:%ld\n",
     cpu_time_geo.d,cpu_time_geo.h,
	 cpu_time_geo.m,cpu_time_geo.s);*/  //Print out before
    fprintf(f_ou,"CPU_Run_Time %ld:%ld:%ld:%ld\n",
     cpu_time_run.d,cpu_time_run.h,
	 cpu_time_run.m,cpu_time_run.s);
    fprintf(f_ou,"CPU_Total_Time %ld:%ld:%ld:%ld\n",
	 cpu_time_total.d,cpu_time_total.h,
	 cpu_time_total.m,cpu_time_total.s);

  // Print final screen output
    printf("\n\n ****** End Of Integration Steps ******\n");
    printf("\tTotalStretch  %.8lG\n", total_stretch);
    printf("\tRHS^2  %.4lG\n", convergence);
    printf("\n CPU_Run_Time %ld:%ld:%ld:%ld\n",
     cpu_time_run.d,cpu_time_run.h,
	 cpu_time_run.m,cpu_time_run.s);
	 printf("-------------------------------\n");
    printf(" CPU_Total_Time %ld:%ld:%ld:%ld\n\n",
     cpu_time_total.d,cpu_time_total.h,cpu_time_total.m,cpu_time_total.s);

    /*    UNCOMMENT to write .m (mathematica) file
	f_ou = fopen("last.m","w");
	export_graphic_complex(f_ou);
	fclose(f_ou);
	*/



}
/*******************************************************************/
// 6. end:
//    write CPU time in hour:minute:second format
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
}
/*******************************************************************/
/************************************************************************/
