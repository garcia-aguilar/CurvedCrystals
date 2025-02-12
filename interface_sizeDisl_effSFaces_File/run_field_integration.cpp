#include "common.h"

double DT;
double rhs_at_i[MAX_SIZE];

void init_random();                     // 5.1.
void run();                             // 5.2

/************************************************************************/
// 5. run_field_integration
void run_field_integration()
{
    printf("\n ****** Run field integration ******\n");
  // initial random stress values
    init_random();                      // 5.1.


  // time step
	//DT = run_time/num_of_iteration;     // change to read DT**
	DT = time_step;
	printf("\tTotal Run_time = %.2lG\n",DT*num_of_iteration);
	fprintf(f_ou,"Total T_run = %.2lG\n",DT*num_of_iteration);
	// check against DT from lmin (Euler)
	double max_DT;
	max_DT = DT_tolerance*dl_min*dl_min;
	if(DT > max_DT){
        printf("\n\tDT from input (%G) > Dt_max (%.4G)", DT, max_DT);
        DT = max_DT;
        printf("\n\tDT changed to DT = %.4G", DT);
        printf("\n\tT_run now = %.4lg\n",DT*num_of_iteration);
	}
	fprintf(f_ou,"DT  %.4lG\n", DT);
	fprintf(f_ou,"Max_DT  %.4lG\n", max_DT);
	fprintf(f_ou,"T_run_now  %.3lG\n", DT*num_of_iteration);

	printf("\n");
	run();                              // 5.2


}
/*******************************************************************/
// 5.1. init_random:
//                  assing random initial values of the field phi to each point (constrained to total_stress)
void init_random()
{
	double ran2(long *),mid_phi;                    // 5.1.1.
	long i;
    double total_stress = 0.;

	mid_phi = init_stress/total_area;

	for (i=0; i<num_of_meshpoint; i++){
        vertex[i].phi = mid_phi+phi_range*(1-2*ran2(&seed))/2;
		//vertex[i].phi = mid_phi;                                    //+++ uniform init
		total_stress += vertex[i].phi*vertex[i].area;
	}

	printf("\tStress input = %.11lg, stress total = %.11lg\n",init_stress,total_stress);
	printf("\t\tDifference = %lg\n",init_stress-total_stress);       //+++ chk
	printf("\t -> Correcting to stress input\n");

	//check total_stress = input
    double phi_corr = (init_stress-total_stress)/total_area;
    total_stress = 0.;
    for (i=0; i<num_of_meshpoint; i++){
//        double newS = vertex[i].phi+phi_corr;
        vertex[i].phi += phi_corr;
		total_stress += vertex[i].phi*vertex[i].area;
	}
    printf("\tStress input = %.11lG, shifted stress = %.8lG\n",init_stress,total_stress);
	printf("\t\tDifference now = %.4lG\n",init_stress-total_stress);

  // Write run variables
    fprintf(f_ou, "\n<<< Integration Info >>>\n");
    fprintf(f_ou,"InitStressCalculated  %.8G\t\t(should be ~0)\n", total_stress);

}

/*******************************************************************/
/*************************************************************/
// 5.1.1. ran2:
//          get a random number
double ran2(long *idum)
{
  	int j;
  	long k;
  	static long idum2=123456789;
  	static long iy=0;
  	static long iv[NTAB];
  	double temp;

  	if (*idum <= 0) {
    	if (-(*idum) < 1) *idum=1;
    	else *idum = -(*idum);
    	idum2 = (*idum);
    	for (j=NTAB+7; j>=0; j--) {
      		k = (*idum)/IQ1;
      		*idum = IA1*(*idum-k*IQ1)-k*IR1;
      		if (*idum < 0) *idum += IM1;
      		if (j < NTAB) iv[j] = *idum;
    	}
    	iy = iv[0];
  	}
  	k = (*idum)/IQ1;
  	*idum = IA1*(*idum-k*IQ1)-k*IR1;
  	if (*idum < 0) *idum += IM1;
  	k = idum2/IQ2;
  	idum2 = IA2*(idum2-k*IQ2)-k*IR2;
  	if (idum2 < 0) idum2 += IM2;
  	j = iy/NDIV;
  	iy = iv[j]-idum2;
  	iv[j] = *idum;
  	if (iy < 1) iy += IMM1;
  	if ((temp = AM*iy) > RNMX) return RNMX;
  	else return temp;
}
/*************************************************************/
/*******************************************************************/
// 5.2. run:
//           run integration for N steps and save intermediate configuration and time data
void run()
{
    void get_rhs(double *);                             // 5.2.1
    void write_hi_line(FILE *, long);                   // 5.2.3
    void one_step();                                    // 5.2.4


	long t;
	long save_time_data = floor(num_of_iteration/TIME_DATA_POINTS);        // save only TIME_DATA_POINTS number in hi.dat
	long Nexport = floor(num_of_iteration/FIELD_CONF_FILES);        // export only FIELD_CONF_FILES number of s-t files

	FILE *f_hi;

	f_hi = fopen("hi.dat","w");

	get_rhs(rhs_at_i);                                       // 5.2.1

  //  Include progress bar on ou.dat to check remotely
    int progress_bar_count = 0;
	fprintf(f_ou, "\n --- Progress: [ 10 .");
	fflush(f_ou);

	for (t=0; t<num_of_iteration; t++){
		if (t%Nexport==0){
			//sprintf(f_na,"gc_%ld.m",t);       //+++ .m -> .dat
			//f_ou = fopen(f_na,"w");           //+++ .m-> .dat
			export_conf(t,Nexport);                         // 5.2.2
			//export_graphic_complex(f_ou); // ++++ check no .m files
			//fclose(f_ou); //++++ .m -> .dat

          //**// Add divergence stop
			double total_stretch=0;
            for(long i=0; i<num_of_meshpoint; i++){
                total_stretch += vertex[i].phi*vertex[i].phi*vertex[i].area;
            }
            if (total_stretch > DIVERGENCE_STOP){
                printf("   Energy  > %.2lG  -> STOP\n", DIVERGENCE_STOP);
                double divStepsPerCent = t*100/num_of_iteration;
                printf("\t divergent before %ld steps  (%.3lf %%)\n", t, divStepsPerCent);
                fprintf(f_ou, "\n\n Energy  > %.2lG  -> STOP\n\n________________EXIT 1____________________", DIVERGENCE_STOP);
                printf("\n________________EXIT 1____________________\n\n");
                exit (EXIT_FAILURE);
            }

		}
		if (t%save_time_data==0){
            write_hi_line(f_hi,t);                              // 5.2.3
            progress_bar_count = progress_bar_on_ou(t, num_of_iteration, progress_bar_count);            // 3.6.3b (too)  IN: common.cpp
          //**// Add convergence_stop
            if (convergence < convergence_stop){
                fprintf(f_ou, " .... ]  o.O \n ");
                printf("   Numerical RHS^2 < %.2lG  -> Convergence\n", convergence_stop);
                fprintf(f_ou,"\n\t  >>>>>> Numerical RHS^2 < %.2lG  -> Convergence\n", convergence_stop);
                double convStepsPerCent = t*100/num_of_iteration;
                printf("\t convergence reached after %ld steps  (%.1lf %%)\n", t, convStepsPerCent);
                break;
            }
        }
		progress_bar(t, num_of_iteration);                      // 3.6.3 (too)  IN: common.cpp
		one_step();                                             // 5.2.4
		get_rhs(rhs_at_i);                                       // 5.2.1
	}

	if(t==num_of_iteration){
        fprintf(f_ou, " 0! ]  ;) \n "); fflush(f_ou);
        printf("   End of iterations with RHS^2 = %.2lG < %.2lG  (**Check convergence**)\n", convergence, convergence_stop);
        fprintf(f_ou, "\n\t   >>>>>>   End of iterations with RHS^2 = %.2lG < %.2lG  (*^^^*Check convergence*^^^*)\n", convergence, convergence_stop);
    }


	fclose(f_hi);
}
/*******************************************************************/
/*************************************************************/
// 5.2.1. get_rhs:
//          Calculate the right-hand side of the Allen-Cahn equation and current total convergence
void get_rhs(double *rhs)
{
	double laplace(long);                   // 5.2.1.1
	double sigma = 1.0;
	long i;
	convergence = 0.0;

	// Calculate the right-hand side of the Allen-Cahn equation
	for (i=0; i<num_of_meshpoint; i++){
		rhs[i] = sigma*laplace(i)-vertex[i].S;
		convergence += rhs[i]*rhs[i]*vertex[i].area;
    }
}
/*************************************************************/
/****************************************************/
// 5.2.1.1 laplace:
//          calculate laplacian for a triangular mesh
double laplace(long i)
{
	double lb=0;
	long j;

	for (j=0; j<vertex[i].num_of_neighbors; j++){
		lb += vertex[i].weight[j]*(vertex[vertex[i].neighbor[j]].phi-vertex[i].phi);
	}
	return lb;
}
/****************************************************/
/*************************************************************/
// 5.2.2. export_conf:
//              save a s-t(N_step).dat file with field information and local convergence at each mesh point
void export_conf(long t, long period)
{
	char f_na[32];
	long i;


	FILE *f_ou;

	if (t%period!=0) return;
	sprintf(f_na,"s-t%ld.dat",t);
	if (t<0) sprintf(f_na,"last.dat");

	f_ou = fopen(f_na,"w");

	for (i=0; i<num_of_meshpoint; i++){
        //double conv = (sigma*laplace(i)-vertex[i].S)*(sigma*laplace(i)-vertex[i].S)*vertex[i].area;
        double conv = rhs_at_i[i]*rhs_at_i[i]*vertex[i].area;
		fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.13f\t%.10lg\n",
		vertex[i].x,
		vertex[i].y,
		vertex[i].z,
		vertex[i].phi,
		conv);
	}
	fclose(f_ou);
}
/*************************************************************/
// 5.2.3. write_hi_line:
//          write one line at T=t(steps) of the hi.dat file
void write_hi_line(FILE *f_ou, long t)
{
	long i,some_v;
	some_v =num_of_meshpoint/3;

	double total_stress = 0.0;
	double total_stretch = 0.0;

	for (i=0; i<num_of_meshpoint; i++){
        total_stress += vertex[i].phi*vertex[i].area;
        total_stretch += vertex[i].phi*vertex[i].phi*vertex[i].area;
	}
	fprintf(f_ou,"%.10f\t%.10f\t%.10f\t%.10f\t%.10g\t%.10g\n",t*DT,vertex[0].phi,vertex[some_v].phi,total_stress,total_stretch,convergence);
}
/*************************************************************/
// 5.2.4. one_step:
//
void one_step()
{
    void euler();               // 5.2.4.1
    void rk2();                 // 5.2.4.2
    void rk4();                 // 5.2.4.3

	switch(integration_method){
  	case 1:
		euler();
    	break;
  	case 2:
		rk2();
    	break;
  	case 3:
		rk4();
    	break;
	}
}
/*************************************************************/
/****************************************************/
// 5.2.4.1 euler:
//          Euler integration
void euler()
{
	double h=DT;
	double rhs[MAX_SIZE];
	long i;

	get_rhs(rhs);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi += h*rhs[i];
	}
}
/****************************************************/
// 5.2.4.2 rk2:
//          Runge-Kutta 2-step integration
void rk2()
{
	double h=DT/2;
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	long i;

	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi;
	}

	get_rhs(rhs1);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*rhs1[i];
	}

	get_rhs(rhs2);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+h*(rhs1[i]+rhs2[i]);
	}
}
/****************************************************/
// 5.2.4.3 rk4:
//          Runge-Kutta 4-step integration
void rk4()
{
	double h=DT/6;
	double rhs0[MAX_SIZE];
	double rhs1[MAX_SIZE];
	double rhs2[MAX_SIZE];
	double rhs3[MAX_SIZE];
	double rhs4[MAX_SIZE];
	long i;

	for (i=0; i<num_of_meshpoint; i++){
		rhs0[i] = vertex[i].phi;
	}

	get_rhs(rhs1);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs1[i];
	}

	get_rhs(rhs2);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+0.5*DT*rhs2[i];
	}

	get_rhs(rhs3);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+DT*rhs3[i];
	}

	get_rhs(rhs4);

	for (i=0; i<num_of_meshpoint; i++){
		vertex[i].phi = rhs0[i]+h*(rhs1[i]+2*rhs2[i]+2*rhs3[i]+rhs4[i]);
	}
}
/****************************************************/
/*************************************************************/
/*******************************************************************/
