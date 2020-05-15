/* BROWNIAN DYNAMICS OF SELF-PROPELLED PARTICLES IN A ONE-DIMENSIONAL OFF-LATTICE SYSTEM */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <signal.h>

#define MOD(a,b) ((((a)%(b))+(b))%(b))

/* POTENTIAL FUNCTION */
float LJ(float epsilon, int sigma, float Vrange, float x) {
  float result;
  if (x<Vrange) {result = 4*epsilon*( pow((sigma/x),12) - pow((sigma/x),6) - pow((1/(float)2.5),12) + pow((1/(float)2.5),6) );}
  else {result=0;}
  return result;
}

/* GRADIENT OF THE POTENTIAL FUNCTION */
float LJprime(float epsilon, int sigma, float Vrange, float x) {
  float result;
  if (x<Vrange) {result = 24*epsilon*pow(sigma,6)*(-2*pow(sigma,6)+pow(x,6))/pow(x,13) - 24*epsilon*pow(sigma,6)*(-2*pow(sigma,6)+pow(Vrange,6))/pow(Vrange,13);}
  else {result=0;}
  return result;
}

/* Mod (remainder) for floats. */
float nfmod(float a,float b)
{
    return a - b * floor(a / b);
}

int main(int argc, char **argv){

	srand(time(NULL));					/* Initialize random number generator */

	/* PARAMETERS */
	int DEBUG;

		// Debug variable. 
		// 	-1 = no info
		// 	 0 = sim progress
		// 	 1 = sim status
		// 	 2 = particle dynamics
		// 	 3 = particle dynamics (more)
		// 	 6 = collective dynamics
		// 	 7 = "ONLY MEASURE SIZE..."
		// 	 8 = SPATIAL CORRELATIONS

	char simName[100];    // Name of the simulation.
	char output_dir[192];                 /* Output directory */
	int L;						/* Lattice size */
	float alpha;		/* Tumbling rate: initial value, interval, final value */
	float fi;			/* Density: initial value, interval, final value */
	int N;						/* # of particles */
	float Tmax;						/* Total simulation time */
	float Tint;						/* Snapshot measurement interval */
	float Vrange;						/* Range of the LJ potential */
	float epsilon;	/* LJ potential intensity parameter */
	int sigma;						/* LJ potential particle diameters parameter */
	float Fp;			/* Propulsion force magnitude */
	float beta;		/* 1/KbT: initial value, interval, final value */
	float mu;			/* Friction coefficient: initial value, interval, final value */
	float Dt;						/* Translational diffusion coefficient*/
	int CMOB;			/* Cluster mobility selector: 1 = motile clusters - 0 = still clusters */
	int INIT_STATE;					/* Initial state of the system selector: 0=random 1=coarsened 2=gas */
	float cluster_cutoff;					/* Cutoff distance for considering particles belong to same cluster */
	float dt;           					/* Time step */
	int measures;		/* Number of measures */ 


	/* READ PARAMETERS */
	scanf("Simulation name: %s\n", simName);
	scanf("Output path: %s\n", output_dir);
	scanf("Debug: %d\n", & DEBUG);
	scanf("Lattice size (L): %d\n", & L);
	scanf("Tumbling rate (alpha): %f\n", & alpha);
	scanf("Density (phi): %f\n", & fi);
	scanf("Simulation time (time units): %f\n", & Tmax);
	scanf("Number of measures: %d\n", & measures);
	scanf("LJ pot. range (cutoff): %f\n", & Vrange);
	scanf("LJ pot. intensity (epsilon): %f\n", & epsilon);
	scanf("LJ pot.  diameter (sigma): %d\n", & sigma);
	scanf("Propulsion force (v, Fp): %f\n", & Fp);
	scanf("Friction coeficient (mu): %f\n", & mu);
	scanf("Temperature (beta): %f\n", & beta);
	scanf("Cluster mobility (CMOB): %d\n", & CMOB);
	scanf("Cluster cutoff: %f\n", & cluster_cutoff);
	scanf("Initial state: %d\n", & INIT_STATE);
	scanf("Time Step: %f\n", & dt);
	scanf("Snap measuring interval: %f\n", & Tint);

	printf("Simulation name: %s\n", simName);
	printf("Output path: %s\n", output_dir);
	printf("Debug: %d\n",  DEBUG);
	printf("Lattice size (L): %d\n",  L);
	printf("Tumbling rate (alpha): %f\n",  alpha);
	printf("Density (phi): %f\n",  fi);
	printf("Simulation time (time units): %f\n",  Tmax);
	printf("Number of measures: %d\n", measures);
	printf("LJ pot. range (cutoff): %f\n",  Vrange);
	printf("LJ pot. intensity (epsilon): %f\n",  epsilon);
	printf("LJ pot.  diameter (sigma): %d\n",  sigma);
	printf("Propulsion force (v, Fp): %f\n",  Fp);
	printf("Friction coeficient (mu): %f\n",  mu);
	printf("Temperature (beta): %f\n",  beta);
	printf("Cluster mobility (CMOB): %d\n",  CMOB);
	printf("Cluster cutoff: %f\n",  cluster_cutoff);
	printf("Initial state: %d\n",  INIT_STATE);
	printf("Time Step: %f\n", dt);
	printf("Snap measuring interval: %f\n", Tint);


		/* DEBUG */ if (DEBUG >= 1) {printf("\n\n\n|||||||||||||||||||| SIMULATION START ||||||||||||||||||||\n\n\n");}


	epsilon=epsilon/24;

	/* VARIABLES */
	N = L*fi;
	Dt = mu/beta;
	float r_potmin = 1.123149; 			// This is for the shifted force for sig=1 and potcut=2.5, for unshifted: sigma * pow(2,1/6);
	float positions[N];           /* Particles' positions */
	int partOrder[N];				/* Array idx = particle ID. Array element = particles order. */
	int orderedPart[N];			/* Array idx = particles order. Array element = particle ID. */
	int directions[N];            /* Particles' swimming directions (1=+x; 0=-x;) */
	float trueDirections[N];			/* Particles' resulting velocities */
	int i, j;					            /* Counter */
	float V[N];					          /* Potential on each particle*/
	float Vprime[N];				      /* Gradient of the potential on each particle*/
	int step;					            /* Time step counter */
	int clusters[N];				      /* Particles' clusters */
	float LBclusters[N/2];			  /* Vector with the left boundaries of the clusters */
	float RBclusters[N/2];			  /* Vector with the right boundaries of the clusters */
	int Nsize[N];					        /* Vector with the number of clusters of each size */
	int Ndir[3]={0,0,0};				  /* Vector with the number of clusters in each direction : +, -, 0 */
	int Ndir_swim[2]={0,0};			  /* Vector with the number of particles swimming in each direction : +, - */
	int Ndir_res[3]={0,0,0};			/* Vector with the number of particles moving in each direction : +, -, 0 */
	float C[L/2];					        /* Correlations */
	float E=0;					          /* Energy of the system (normalized by the number of particles) */
	clock_t tt;					          /* Real time */
	clock_t t_0=clock();				  /* Real starting time */

	// Checks.
	if (sigma != 1 && Vrange != 2.5)
	{
		printf("Bad r_potmin! Recompute force zero (r_potmin) for given values of sigma and cutoff.\n");
		exit(1);
	}

	if (fi > 1/(r_potmin))
	{
		printf("Density too high! Particle overlap! Maximum density is: %f\n", 1/r_potmin);
		exit(1);
	}

	/* FILE PARAMETERS */

	// If it doesn't exist, create it.
	struct stat st = {0};
	if (stat(output_dir, & st) == -1) {mkdir(output_dir, 0777);}

	char newdir[192];									/* Directory for the new simulation */
	sprintf(newdir, "%s/sim_a%.3f_f%.3f_t%010.2f_L%.5d_D%05.2f_Fp%.2f_eps%07.4f_CMOB%d_IS%d_tint%06.1f_dt%7.5f", output_dir, alpha, fi, Tmax, L, Dt, Fp, epsilon, CMOB, INIT_STATE,Tint,dt);	

	// If it doesn't exist, create it.
	struct stat st_bis = {0};
	if (stat(newdir, & st_bis) == -1) {mkdir(newdir, 0777);}

	FILE *sc;										/* Saves the number of clusters of the same size overall */
	char buffersc[192];
	sprintf(buffersc, "%s/sizedistr.dat", newdir);

	FILE *scn;										/* Saves the number of clusters of the same normalized size overall */
	char bufferscn[192];
	sprintf(bufferscn, "%s/sizedistr_norm.dat", newdir);

	FILE *dc;										/* Saves the number of clusters for each direction overall */
	char bufferdc[192];
	sprintf(bufferdc, "%s/dirdistr.dat", newdir);

	FILE *pdc;										/* Saves the number of particles swimming in each direction overall */
	char bufferpdc[192];
	sprintf(bufferpdc, "%s/swim_dirdistr.dat", newdir);

	FILE *rpdc;										/* Saves the number of particles moving (result of LJ) in each direction overall */
	char bufferrpdc[192];
	sprintf(bufferrpdc, "%s/result_dirdistr.dat", newdir);

	FILE *nrj;										/* Saves the energy of the system */
	char buffernrj[192];
	sprintf(buffernrj, "%s/energy.dat", newdir);
	nrj=fopen(buffernrj, "wb");

	FILE *nc;										/* Saves the total number of clusters */
	char buffernc[192];
	sprintf(buffernc, "%s/nclusters.dat", newdir);
	nc=fopen(buffernc, "wb");
	
	FILE *snap;										/* Saves positions and directions of the particles for initial times */
	char buffersnap[192];
	sprintf(buffersnap, "%s/system_beg.dat", newdir);
	snap=fopen(buffersnap, "wb");
	
	FILE *snapbis;									/* Saves positions and directions of the particles for the steady state */
	char buffersnapbis[192];
	sprintf(buffersnapbis, "%s/system_ss.dat", newdir);
	snapbis=fopen(buffersnapbis, "wb");
	
	FILE *corr;										/* Saves site-to-site occupancy correlations */
	char buffercorr[192];
	sprintf(buffercorr, "%s/corr.dat", newdir);
	corr=fopen(buffercorr, "wb");

	FILE * pars;
	char bufferspars[192];								/* Saves the parameters of this simulation */
	sprintf(bufferspars, "%s/parameters.dat", newdir);

	FILE *logfile;									/* Log file */
	char bufferlog[192];
	sprintf(bufferlog, "%s/out.log", newdir);
	logfile  = fopen(bufferlog, "w");

	pars=fopen(bufferspars, "wb");
	fprintf(pars, "Debug variable: %d (0 = no info printed - 1 = some info printed - 2 = all info printed)\n", DEBUG);
	fprintf(pars, "Lattice size: %d\n", L);
	fprintf(pars, "Tumbling rate: %f\n", alpha);
	fprintf(pars, "Particle density: %f\n", fi);
	fprintf(pars, "Number of particles: %d\n", (int)N);
	fprintf(pars, "Total simulation time: %f\n", Tmax);
	fprintf(pars, "Number of measures: %d\n", measures);
	fprintf(pars, "Range of the LJ potential: %f\n", Vrange);
	fprintf(pars, "LJ potential intensity parameter (epsilon): %f\n", epsilon);
	fprintf(pars, "LJ potential particle diameters parameter (sigma): %d\n", sigma);
	fprintf(pars, "Propulsion force: %f\n", Fp);
	fprintf(pars, "Beta - Inverse of the temperature energy: %f\n", beta);
	fprintf(pars, "Mu - Friction coefficient / mobility: %f\n", mu);
	fprintf(pars, "Tranlational diffusivity: %f\n", Dt);
	fprintf(pars, "Mobility of the clusters: %d\n", CMOB);
	fprintf(pars, "Cluster cutoff distance: %f\n", cluster_cutoff);
	fprintf(pars, "Initial state of the system: %d (0 = random - 1 = gas - 2 = coarsened)\n", INIT_STATE);
	fprintf(pars, "Time step: %f\n", dt);
	fprintf(pars, "Snap measuring interval: %f\n", Tint);
	fclose(pars);

	if (DEBUG>0) {
	printf("\nSize of the system: %d\n# of particles: %d\nTumbling rate: %.3f\nTotal simulation time: %f\nMeasuring interval: %f\n",L,N,alpha,Tmax,Tmax/measures );
	}

	/* INITIAL STATE OF THE SYSTEM */

	/* DEBUG */ if (DEBUG>=1) 
					{
						printf("\n\n\n==========> INITIALIZING THE SYSTEM\n");
					}
	/* DEBUG */	if (DEBUG==3) 
					{
						if (INIT_STATE==0) {printf(" - RANDOM\n");}
						else if (INIT_STATE==2) {printf(" - SINGLE CLUSTER\n");}
						else if (INIT_STATE==1) {printf(" - DILUTE GAS\n");}
					}

	/* DEBUG INITIAL STATE: positions and directions by hand. */
	if (INIT_STATE==-1)
	{
		positions[0] = L/2;
		positions[1] = L/2 + 1.5;
		positions[2] = L/10;
		positions[3] = L-L/10;		
		directions[0] = -1;
		directions[1] = 1;
		directions[2] = -1;
		directions[3] = -1;
		trueDirections[0] = directions[0];
		trueDirections[1] = directions[1];
		trueDirections[2] = directions[2];
		trueDirections[3] = directions[3];
	}


	/* RANDOM INITIAL STATE */
	if (INIT_STATE==0) {
		for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
			float r1=floor(rand()/(float)RAND_MAX*(L-0.000001));				/* Random position [0, L) */

				if (DEBUG==3) {printf("Position %d: %.3f\n",i+1,r1);}

			int busy=0;
			for (j=0; j<i; j++) {
				float ddist=r1-positions[j];
				if (ddist<0) {ddist=-ddist;}
				if (ddist>L/2) {ddist=L-ddist;}

				    if (DEBUG==3) {printf("  Position %d: %.3f - Distance: %.3f\n",j+1,positions[j],ddist);}

				if (ddist<sigma) {busy=1; continue;}
			}

				if (DEBUG==3) {printf("Busy: %d\n",busy);}

			if (busy==1) {i=i-1; continue;}
			positions[i]=r1;

			/* Random direction (1=+x; 0=-x;) */
			int r2=rand()%2;
			int r3;
			if (r2==1) {r3=1;}
			else if (r2==0) {r2=-1; r3=-1;}
			directions[i]=r2;
			trueDirections[i]=r3;

		}
	}

	/* GAS INITIAL STATE */
	else if (INIT_STATE==1) {
		for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
			float interval=L/(float)N;
			positions[i]=i*interval;						/* Ordered particles in a dilute gas at constant distance */
			  
			/* Random direction (1=+x; 0=-x;) */
			int r2=rand()%2;							
			int r3;
			if (r2==1) {r3=1;}
			else if (r2==0) {r2=-1; r3=-1;}
			directions[i]=r2;
			trueDirections[i]=r3;
		}
	}
	
	/* COARSENED INITIAL STATE */
	else if (INIT_STATE==2) 
		{			
		for (i=0; i<N; i++) 
			{						/* Particles are N, ordered from 1 to N */
		  		float start_pos=L/(float)2-N/(float)2;
		  		positions[i] = start_pos+i;						/* Ordered particles in a single cluster in the middle of the system */

				/* Random direction (1=+x; 0=-x;) */
				int r2=rand()%2;							
		  		int r3;
		  		if (r2==1) {r3=1;}
		  		else if (r2==0) {r2=-1; r3=-1;}
				directions[i]=r2;
		  		trueDirections[i]=r3;
			}
		}

	/* EQUIDISTANT PARTICLES */
	else if (INIT_STATE==3) 
		{
		for (int i = 0; i < N; i++) 
			{						
		  		positions[i] = (float)i/fi;						

				/* Random direction (1=+x; 0=-x;) */
				int r2=rand()%2;
		  		int r3;
		  		if (r2==1) {r3=1;}
		  		else if (r2==0) {r2=-1; r3=-1;}
				directions[i]=r2;
		  		trueDirections[i]=r3;
			}
		}

	/* Initiaize cluster arrays */
	for (int i = 0; i < N; ++i)
		{
			clusters[i]=0;
			Nsize[i]=0;
		}

		/*DEBUG*/ 	if (DEBUG > 1)
						{
							printf("\nIntial state of the system:\n");
							for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
							printf("\n");
							for (i=0; i<N; i++) {printf("%d ",directions[i]);}
							printf("\n");
						}


	/* Compute particle order */

	/* DEBUG */ if (DEBUG >= 1)
				{
					printf("\n\n========> Compute particle order");
				}

	// Initialize order arrays.
	float orderedPos[N];
	int orderedParticles[N];
	int particlesOrder[N];
	for (int i = 0; i < N; ++i)
	{
		orderedPos[i] = positions[i];	// This is not used directly.
      	orderedParticles[i]=i;			// The ID of the sorted particles, from closer to 0 to further..
      	particlesOrder[i]=i;			// The order of the particle i.
	}

	// Sort.
	for(int i = 0; i < N-1; i++)
	{
	  int Imin = i;
	  for(int j = i+1; j < N; j++)
	  {
	     if( orderedPos[j] < orderedPos[Imin] )
	     {
	        Imin = j;
	     }
	  }
	  float temp = orderedPos[Imin];
	  orderedPos[Imin] = orderedPos[i];
	  orderedPos[i] = temp;

	  float tempI = orderedParticles[Imin];
	  orderedParticles[Imin] = orderedParticles[i];
	  orderedParticles[i] = tempI;
	}

	for (int i = 0; i < N; ++i)
	{
	  particlesOrder[orderedParticles[i]]=i;
	}

	// Neighbour list: each element contains de neightbours of each particle. 
		// neighbours[i][0] the left one.
		// neighbours[i][1] the right one.

	/* DEBUG */ if (DEBUG >= 1)
				{
					printf("\n\n========> Neighbour list");
				}

	// Initialize neighbour array: if -1 no neighbours.
	int neighbours[N][2];
	for (int i = 0; i < N; ++i)
	{
		neighbours[i][0]=-1;
		neighbours[i][1]=-1;
	}

	for (int i = 0; i < N; ++i)
	{
		if (i != orderedParticles[MOD(particlesOrder[i]-1,N)])
		{
			neighbours[i][0] = orderedParticles[MOD(particlesOrder[i]-1,N)];
		}

		if (i != orderedParticles[MOD(particlesOrder[i]+1,N)])
		{
			neighbours[i][1] = orderedParticles[MOD(particlesOrder[i]+1,N)];
		}
	}

		/* DEBUG */	if (DEBUG > 2)
					{
						printf("\nPositions: \n");
						for(int n = 0 ; n < N; n++ ) {   
						  printf("%f ", positions[n]);
						}

						printf("\nOrdered Positions: \n");
						for(int n = 0 ; n < N; n++ ) {   
						  printf("%f ", orderedPos[n]);
						}

						printf("\nThe particles order: \n");
						for(int n = 0 ; n < N; n++ ) {   
						  printf("%d ", particlesOrder[n]);
						}

						printf("\nThe ordered particles: \n");
						for(int n = 0 ; n < N; n++ ) {   
						  printf("%d ", orderedParticles[n]);
						}
						printf("\n");
						printf("\nThe particles neighbours: \n");
						for(int n = 0 ; n < N; n++ ) {   
							printf("Particle %d. Neighbours: %d and %d \n", n, neighbours[n][0],neighbours[n][1]);
							printf("Position %f. Neighbours: %f and %f \n", positions[n], positions[neighbours[n][0]], positions[neighbours[n][1]]);
							printf("\n");
						}
						printf("\n");
					}
	




	/* TIME LOOP - DYNAMICS*/						

	for (step=0; step < (int)ceil((Tmax/dt)); step++) {

		// For only output of the last iteration, uncomment this.
		//logfile  = freopen(bufferlog, "w", logfile);


			/* DEBUG */ if (DEBUG >= 1) {fprintf(logfile, "\n\n\n==========> TIME LOOP - t = %d\n",step+1);}


		/* INDIVIDUAL DYNAMICS */					

			/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n========> PARTICLE DYNAMICS");}

		/* Reset the particle-to-cluster tagger */
		for (i=0; i<N; i++) {clusters[i]=0;}

		/* Random order of movement. */				

			/* DEBUG */ if (DEBUG > 2) {fprintf(logfile, "\nMvmnt order (1st=0): \n");}

		int orderofmovement[N];	
		for (i=0; i<N; i++) 
		{
			/* Random particle to move at turn i */
			int par = rand()%N+1;

			/* Save particle par in turn i (even if repetition happens) */
		  	orderofmovement[i] = par;
			
			/* DEBUG */ if (DEBUG >= 2) {fprintf(logfile, "%d ", par-1);}
		}





	/* PARICLE LOOP. Change each particle in the above order. */

	int k;
	int ptm = 0;
	for (k=0; k < N; k++) 
	{
		/* Particle to move now */
	  	ptm = orderofmovement[k]-1;

	  	// Check particle is inside the box.
	  	if (positions[ptm] < 0 || positions[ptm] >= L)
		{
			printf("\n OUT OF THE BOX! \npositions[%d] = %f\n", ptm, positions[ptm]);
			exit(1);
		}

	  		/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n========> PATICLE LOOP - ptm = %d\n",ptm);}
	  	
	  		/* DEBUG */ if (DEBUG > 2) 
	  					{
	  						fprintf(logfile, "\n\nParticle (1st=0): %d --------------------------------------\nNeighs: %d and %d\n\n",ptm,neighbours[ptm][0],neighbours[ptm][1]);
	  					}


	  	/* TUMBLE */

			/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n======> 1. Tumble\n");}

	  	/* Random number to compare to alpha */
		int prob=rand();
		if (prob<=alpha*RAND_MAX) 
		{
			int tmbl=rand();				/* Random number to redefine direction */
			if (tmbl<=0.5*RAND_MAX) { directions[ptm] = -directions[ptm]; }
		}

			/* DEBUG */ if ( DEBUG > 2 )
					{
						// Initialize.
					  	int orderedParticlesConc[2*N];

					  	for (int i = 0; i < 2*N; ++i)
					  	{
							orderedParticlesConc[i] = orderedParticles[i%N];
					  	}

						float orderedPos2[N];
						int orderedParticles2[N];
						int particlesOrder2[N];
						for (int i = 0; i < N; ++i)
						{
							orderedPos2[i] = positions[i];	// This is not used directly.
					      	orderedParticles2[i]=i;			// The ID of the sorted particles, from closer to 0 to further..
					      	particlesOrder2[i]=i;			// The order of the particle i.
						}

						// Sort.
						for(int i = 0; i < N-1; i++)
						{
						  int Imin = i;
						  for(int j = i+1; j < N; j++)
						  {
						     if( orderedPos2[j] < orderedPos2[Imin] )
						     {
						        Imin = j;
						     }
						  }
						  float temp = orderedPos2[Imin];
						  orderedPos2[Imin] = orderedPos2[i];
						  orderedPos2[i] = temp;

						  float tempI = orderedParticles2[Imin];
						  orderedParticles2[Imin] = orderedParticles2[i];
						  orderedParticles2[i] = tempI;
						}

						int j=0;
						for (int i = 0; i < 2*N; ++i)
						{
							if (j < N)
							{
								if (orderedParticlesConc[i] == orderedParticles2[j]){ j++; }
								else { j=0; }
							}
							else if (j == N)
							{
								fprintf(logfile, "No Q-jumps.\n");
								break;
							}
						}

						if (j != N)
						{
							fprintf(logfile, "\n¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ QUANTUM JUMP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
							fprintf(logfile, "\nThe ordered particles: \n");
							for(int n = 0 ; n < 2*N; n++ ) {   
							  fprintf(logfile, "%d ", orderedParticlesConc[n]);
							}
							fprintf(logfile, "\nThe ordered particles 2: \n");
							for(int n = 0 ; n < N; n++ ) {   
							  fprintf(logfile, "%d ", orderedParticles2[n]);
							}
							fprintf(logfile, "\n");
							fprintf(logfile, "\nState: \n");
							for (i=0; i<N; i++) {fprintf(logfile, "%.3f ", positions[i]);}
							fprintf(logfile, "\n");

							//exit(0);
						}
					}

			/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n======> 2. Potential\n");}

		/* POTENTIAL (NEW) */
		float pos = positions[ptm];		/* Position on the array of the particle */
		V[ptm] = 0;						/* Reset the potential on that particle */
		Vprime[ptm] = 0;				/* Reset the gradient of the potential on that particle */

		if (epsilon != 0) 
		{
			int ptmaux = ptm;	
			bool alone = false;
			float distaux = 0;
			int numbLeftNeighs = 0;

			// Left neighbours.
			while ( !alone && numbLeftNeighs < N) 
			{

				// Distance to vth left neighbour.
				distaux = fabsf( pos - positions[neighbours[ptmaux][0]] );

				// PBCs
				if (distaux > L/2) { distaux = L - distaux; }

				// If vth left neighbour is in range, add its contribution to the potential. And go to the (v+1)th left neighbour.
				if ( distaux < Vrange )
				{
					numbLeftNeighs++;
					// printf("\nV(left): %f\n", Vprime[ptm]);
					V[ptm] += LJ(epsilon, sigma, Vrange, distaux);
					Vprime[ptm] += LJprime(epsilon, sigma, Vrange, distaux);

					// Checks
					if (isnan(Vprime[ptm]) || isinf(Vprime[ptm]) ) 
						{
							printf("\n\nLEFT\n");
							printf("ptm: %d\n", ptm); 
							printf("ptmaux: %d\n", ptmaux); 
							printf("pos: %f\n", pos);
							printf("posLneigh[ptmaux]: %f\n", positions[neighbours[ptmaux][0]]);
							printf("distaux: %f\n", distaux);
							printf("Vrange: %f\n", Vrange);
							printf("eps: %f\n", epsilon);
							printf("sigma: %d\n", sigma);
							printf("Vprime[%d -> %d]: %f\n", neighbours[ptmaux][0], ptm, LJprime(epsilon, sigma, Vrange, distaux));
							printf("VprimeTOT[%d]: %f\n\n", ptm, Vprime[ptm]);
							printf("Number of left beighbours: %d\n", numbLeftNeighs);
							printf("1st neighbour\n");
							printf("positionsL[%d]: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
							printf("Rest of neighbours\n");
							printf("positions[ptmaux = %d]: %f\n", ptmaux, positions[neighbours[ptmaux][0]]);
							printf("\nState: \n");
		      				for (i=0; i<N; i++) {printf("%.3f ", positions[i]);}
		      				printf("\n\n");
							exit(1);
						}
					if (ptm == neighbours[ptmaux][0])
						{
							fprintf(logfile, "¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ EH, AUTOREFERENCIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
							//exit(0);
						}

					ptmaux = neighbours[ptmaux][0];
				}

				// If not, the particle is alone :(
				else {alone = true;}

					/*DEBUG */  if (DEBUG > 2) { fprintf(logfile, "Alone: %d\n", alone); }
			}
			
				/*DEBUG */  if (DEBUG > 2) { fprintf(logfile, "No more left neighbours.-------------------\n"); }

			// Repeat for right neighbours.
			ptmaux = ptm;	
			alone = false;
			distaux = 0;
			int numbRightNeighs = 0;

			while ( !alone && numbRightNeighs < N) 
			{
				// Distance to vth right neighbour.
				distaux = fabsf( pos - positions[neighbours[ptmaux][1]] );

				// PBCs
				if (distaux > (float)L/2) { distaux = L - distaux; }

				// If vth right neighbour is in range, add (substract, since its on the right) its contribution to the potential. And go to the (v+1)th right neighbour.
				if (distaux < Vrange)
				{
					numbRightNeighs++;
					// printf("\nV(right): %f\n", Vprime[ptm]);
					V[ptm] += LJ(epsilon, sigma, Vrange, distaux);
					Vprime[ptm] -= LJprime(epsilon, sigma, Vrange, distaux);

					// Checks
					if (isnan(Vprime[ptm]) || isinf(Vprime[ptm]) ) 
						{
							printf("\n\nRIGHT\n");
							printf("ptm: %d\n", ptm); 
							printf("ptmaux: %d\n", ptmaux); 
							printf("pos: %f\n", pos);
							printf("posRneigh[ptmaux]: %f\n", positions[neighbours[ptmaux][1]]);
							printf("distaux: %f\n", distaux);
							printf("Vrange: %f\n", Vrange);
							printf("eps: %f\n", epsilon);
							printf("sigma: %d\n", sigma);
							printf("Vprime[%d -> %d]: %f\n", neighbours[ptmaux][1], ptm, LJprime(epsilon, sigma, Vrange, distaux));
							printf("VprimeTOT[%d]: %f\n\n", ptm, Vprime[ptm]);
							printf("Number of right beighbours: %d\n", numbRightNeighs);
							printf("1st neighbour\n");
							printf("positionsR[%d]: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
							printf("Rest of neighbours\n");
							printf("positions[ptmaux = %d]: %f\n", ptmaux, positions[neighbours[ptmaux][1]]);
							printf("\nState: \n");
		      				for (i=0; i<N; i++) {printf("%.3f ", positions[i]);}
		      				printf("\n\n");

							exit(1);
						}
					if (ptm == neighbours[ptmaux][1])
						{
							fprintf(logfile, "¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡ EH, AUTOREFERENCIA !!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
							exit(1);
						}	

					ptmaux = neighbours[ptmaux][1];

				}

				// If not, the particle is alone :(
				else {alone = true;}

					/*DEBUG */  if (DEBUG > 2) { fprintf(logfile, "Alone: %d\n", alone); }
			}

				/*DEBUG */  if (DEBUG > 2) { fprintf(logfile, "No more right neighbours.-------------------\n"); }
		}

	  		/*DEBUG */ 	if (DEBUG > 2) 
							{
								fprintf(logfile, "\nForce felt by particle %d.\nF: %.10f\n", ptm, -Vprime[ptm]); 
								fprintf(logfile, "V: %.10f\n",V[ptm]);
							}


		/* NEW POSITION */ 

			/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n======> 3. New position\n");}

		float npos;							/* New position */
		int vel=directions[ptm];				/* Velocity of the particle: 1 site/time-step in the direction */

		/* GAUSSIAN WHITE NOISE */
		float U1=(rand()+1)/(float)RAND_MAX;			/* Uniformly distributed RV */
		//if (DEBUG==3) {fprintf(logfile, "U1: %.3f\n",U1);}
		float U2=(rand()+1)/(float)RAND_MAX;			/* Uniformly distributed RV */
		//if (DEBUG==3) {fprintf(logfile, "U2: %.3f\n",U2);}
		float Z1=sqrt(-2*log(U1))*cos(2*M_PI*U2);		/* Normally distributed RV */
		//if (DEBUG==3) {fprintf(logfile, "Z1: %.3f\n",Z1);}
		float Z2=sqrt(-2*log(U1))*sin(2*M_PI*U2);		/* Normally distributed RV */
		//if (DEBUG==3) {fprintf(logfile, "Z2: %.3f\n",Z2);}
		float eta=Z1;
		if (DEBUG==3) {fprintf(logfile, "eta: %.3f\n",eta);}


		// -> Compute new position <- //
		npos = pos + beta*Dt*dt*( Fp*vel - Vprime[ptm] ) + sqrt(2*Dt*dt)*eta;


			/* DEBUG */ if (DEBUG >= 2) 
					{
				        fprintf(logfile, "Old position: %.3f\n",pos);
				        fprintf(logfile, "Velocity: %d\n",vel);
				        fprintf(logfile, "Expected new position: %.3f\n",npos);
					}

		/* Resulting direction */
		trueDirections[ptm] = npos - pos;

			/* DEBUG */ if (DEBUG >= 2) {fprintf(logfile, "True velocity: %.3f\n",trueDirections[ptm]);}



		/* FIND IF MOVEMENT IS ALLOWED (NEW) */
		
			/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n======> 4. Movmnt allowed?\n");}

		/* PBCs (pos = 0 is included, pos = L is excluded). Changed to nfmod, now the jump can be bigger than the box. */
			// if (npos < 0.f) {npos = npos + L;}
			// else if (npos >= (float)L) {npos = npos - L;}
		npos = nfmod(npos, (float)L);
		if (npos == (float)L){npos = 0.f;}

			/* DEBUG */ if (DEBUG >= 2) {fprintf(logfile, "Expected new position (with PBCs): %.3f\n",npos);}

	  	// Compute distance to neighbours and jump distance.
	    float distL = fabsf(pos - positions[neighbours[ptm][0]]);
	    float distR = fabsf(pos - positions[neighbours[ptm][1]]);
	    float distJump = fabsf(npos - pos);

	    	/* DEBUG */ if ( DEBUG > 2)
				    {
				    	fprintf(logfile, "\n");
						fprintf(logfile, "\nState: \n");
						for (i=0; i<N; i++) {fprintf(logfile, "%.3f ", positions[i]);}
						fprintf(logfile, "\n");
				    	fprintf(logfile, "\n\nnpos: %f\ndistJump: %f\ndistL: %f\ndistL-sigma: %f\n", npos, distJump, distL, distL-sigma);
						fprintf(logfile, "distJump: %f\ndistR: %f\ndistR-sigma: %f\n",distJump, distR, distR-sigma);
				    }

		// Check new position is inside box.
		if (npos >= (float)L || npos < 0.f || positions[ptm] >= (float)L || positions[ptm] < 0.f || distR < 1.f || distL < 1.f)
			{
				printf("\n\nBEFORE\n");
				printf("eta: %.3f\n",eta);
				printf("Vprime: %f\n", Vprime[ptm]);
				printf("positions[%d]: %f\n", ptm, positions[ptm]);
				printf("directions: %d\n",vel);
				printf("trueDirections: %f\n", trueDirections[ptm]);
				printf("position[%d] L: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
				printf("position[%d] R: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
				printf("pos: %f\nnpos: %f\ndistJump: %f\ndistL: %f\ndistL-potmin: %f\n", pos, npos, distJump, distL, distL-r_potmin);
				printf("distR: %f\ndistR-potmin: %f\n", distR, distR-r_potmin);
			}	

	    // If moving to the left.
	    if (trueDirections[ptm] < 0.f) 
		    {
		    		/* DEBUG */ if ( DEBUG > 2) {fprintf(logfile, "\nMOVING LEFT\n");}

		    	// Ahora calculamos bien las distancias dado que puede ocurrir que la distancia calculada sea mayor de L/2 pero no cruce las PBCs).
		    	// Por tanto if (dist > L/2) {dist = L-dist;} no es suficiente.

		    	// If the left neighbour is at the other side of the PBC frontier.
		    	if (pos - distL < 0.f)
		    		{ distL = L - distL; }

		    	if (distL < 1)
		    		{
						printf("distL = %f < 1\n", distL);		    	
					}

		    	// If the new position is at the other side of the PBC frontier.
		    	if (pos - distJump < 0.f)
		    		{ distJump = L - distJump; }

		    	// UPDATE POSITION, If the jump is smaller than the distance to its neighbor.
		    	if (distJump <= distL - r_potmin)
		    		{ 
		    			if (DEBUG > 1){	printf("1L. jump < distNeigh\n"); } 

		    			positions[ptm] = npos; 
		    		}

		    	// If the jump is bigger, the particle sticks onto its neighbor only
			    	// IF IT HAS SPACE ON THE RIGHT, if not, it moves to an equidistant distance from its neighs.
			    	// This occurs when distR < r_potmin and distR + distL < 2*r_potmin (as in the initial case) 
		    		// the particle will move increasing distR and if there is a particle on the left it will 
		    		// overlap, so we have to take care of this.		    	
		    	else if (distJump > distL - r_potmin && distL + distR >= 2*r_potmin)
		    		{
		    			if (DEBUG > 1){ printf("2L. jump > distNeigh, space\n"); }

		    			// If r_potmin of the particle is at the other side of the PBC frontier.		    		
		    			if (positions[neighbours[ptm][0]] + r_potmin > (float)L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + r_potmin - (float)L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][0]] + r_potmin < (float)L && positions[neighbours[ptm][0]] + r_potmin >= 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + r_potmin;
			    			}
		    			// This can occur.
		    			else if (positions[neighbours[ptm][0]] + r_potmin == (float)L)
			    			{
			    				positions[ptm] = 0.f;
			    			}
		    			else
		    				{
		    					printf("ERROR! - L\n");
		    					exit(1);
		    				}
		    		}

		    	// If the particle doesn't have space, it moves to a equidistant distance from its neighbours.
		    	else if (distJump > distL - r_potmin && distL + distR < 2*r_potmin)
			    	{
						if (DEBUG > 1){	printf("2L. jump > distNeigh, no space\n"); }
		    			// If (distR+distL)/2 is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][0]] + (distR+distL)/2 > (float)L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + (distR+distL)/2 - (float)L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][0]] + (distR+distL)/2 < (float)L && positions[neighbours[ptm][0]] + (distR+distL)/2 >= 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + (distR+distL)/2;
			    			}
		    			// This can occur.
		    			else if (positions[neighbours[ptm][0]] + (distR+distL)/2 == (float)L)
			    			{
			    				positions[ptm] = 0.f;
			    			}
		    			else
		    				{
		    					printf("ERROR! - L\n");
		    					exit(1);
		    				}			    		
			    	}
		    	else 
		    		{
		    			fprintf(logfile, "\nERROR!\ndistJump: %f\ndistR: %f\ndistL: %f",distJump, distR, distL);
		    			exit(1);
		    		}
		    }

	    // If moving to the right.
	    if (trueDirections[ptm] > 0.f) 
		    {
		    		/* DEBUG */ if ( DEBUG > 2) {fprintf(logfile, "\nMOVING RIGHT\n");}

		    	// If the right neighbour is at the other side of the PBC frontier.
		    	if (pos + distR > (float)L)
		    		{ distR = L - distR; }

		    	if (distR < 1)
		    		{
						printf("distR = %f < 1\n", distR);		    	
					}

		    	// If the new position is at the other side of the PBC frontier.
		    	if (pos + distJump > (float)L)
		    		{ distJump = L - distJump; }

		    	// UPDATE POSITION, If the jump is smaller than the distance to its neighbor.
		    	if (distJump <= distR - r_potmin)
		    		{ 
		    			if (DEBUG > 1){	printf("1L. jump < distNeigh\n"); } 
		    			positions[ptm] = npos; 
		    		}

		    	// If the jump is bigger, the particle sticks onto its neighbor only
			    	// IF IT HAS SPACE ON THE LEFT, if not, it sits there.
			    	// This occurs when distR < r_potmin and distR + distL < 2*r_potmin (as in the initial case) 
		    		// the particle will move increasing distR and if there is a particle on the left it will 
		    		// overlap, so we have to take care of this.
		    	else if (distJump > distR - r_potmin && distL + distR >= 2*r_potmin)
		    		{
		    			if (DEBUG > 1){ printf("2L. jump > distNeigh, space\n"); }
		    			// If r_potmin of the particle is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][1]] - r_potmin < 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - r_potmin + (float)L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][1]] - r_potmin >= 0.f && positions[neighbours[ptm][1]] - r_potmin < (float)L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - r_potmin;
			    			}
		    			// This shouldn't occur.
		    			else if (positions[neighbours[ptm][1]] - r_potmin == (float)L)
			    			{
			    				positions[ptm] = 0.f;
			    				printf("OJO!\n");
			    			}
		    			else
		    				{
		    					printf("ERROR! - R\n");
		    					exit(1);
		    				}		    			

		    		}
		    	
		    	// If the particle doesn't have space, it moves to a equidistant distance from its neighbours.
		    	else if (distJump > distR - r_potmin && distL + distR < 2*r_potmin)
			    	{
						if (DEBUG > 1){	printf("2R. jump > distNeigh, no space\n"); }

		    			// If (distR+distL)/2 is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][1]] - (distR+distL)/2 < 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - (distR+distL)/2 + (float)L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][1]] - (distR+distL)/2 >= 0.f && positions[neighbours[ptm][1]] - (distR+distL)/2 < (float)L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - (distR+distL)/2;
			    			}
		    			// This shouldn't occur.
		    			else if (positions[neighbours[ptm][1]] - (distR+distL)/2 == (float)L)
			    			{
			    				positions[ptm] = 0.f;
			    				printf("OJO!\n");
			    			}
			    	}

		    	else 
		    		{
		    			fprintf(logfile, "\nERROR!\ndistJump: %f\ndistR: %f\ndistL: %f",distJump, distR, distL);
		    			exit(1);
		    		}
		    }

		// Due to floating point precision (e.g. 500 - 0.000013 = 500.000000) sometimes positions[ptm] = 500 occur, so we need to:
		if (positions[ptm] == (float)L){ positions[ptm] = 0.f; }

	    // Check box again.
		if (positions[ptm] >= L || positions[ptm] < 0 || distR < 1.f || distL < 1.f)
			{
				printf("\n\nAFTER\n");
				printf("\n\npositions[%d]: %f\n", ptm, positions[ptm]);
				printf("trueDirections: %f\n", trueDirections[ptm]);
				printf("position[%d] L: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
				printf("position[%d] R: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
				printf("pos: %f\nnpos: %f\ndistJump: %f\ndistL: %f\ndistL-potmin: %f\n", pos, npos, distJump, distL, distL-r_potmin);
				printf("distR: %f\ndistR-potmin: %f\n", distR, distR-r_potmin);
				exit(1);
				//raise(SIGSEGV);
			}	    		

	    	/* DEBUG */ if ( DEBUG > 2)
				    {
					    fprintf(logfile, "\n");
						fprintf(logfile, "\nState: \n");
						for (i=0; i<N; i++) {fprintf(logfile, "%.3f ", positions[i]);}
						fprintf(logfile, "\n");
					    fprintf(logfile, "\n\ndistJump: %f\ndistL: %f\ndistL-r_potmin: %f\n",distJump, distL, distL-r_potmin);
						fprintf(logfile, "distJump: %f\ndistR: %f\ndistR-r_potmin: %f\n",distJump, distR, distR-r_potmin);
					}

		/* Save the direction of the particle: swimming and resulting (every Tmax/measures timesteps) */
		if ((step+1)%(int)ceil((Tmax/(measures*dt)))==0 && step>0) 
		{
			if (directions[ptm]==-1) {Ndir_swim[0]++;}
			else {Ndir_swim[1]++;}
			if (trueDirections[ptm]<0) {Ndir_res[0]++;}
			else if (trueDirections[ptm]>0) {Ndir_res[2]++;}
			else {Ndir_res[1]++;}
		}

	} // End of particle loop.

		/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n========> PARTICLE LOOP (exit)\n");}

		/* DEBUG */ if (DEBUG > 1)
						{
							fprintf(logfile, "\nState of the system after individual dynamics:\n");
							for (i=0; i<N; i++) {fprintf(logfile, "%.3f ", positions[i]);}
							fprintf(logfile, "\n");
							for (i=0; i<N; i++) {fprintf(logfile, "%d ", directions[i]);}
							fprintf(logfile, "\n");
						}





	/* COLLECTIVE DYNAMICS */

	/* EVALUATE CLUSTERS */
	if (DEBUG >= 1) {fprintf(logfile, "\n\n\n==========> CLUSTER DYNAMICS - t = %d\n",step+1);}

	float R_x=L, L_x=0, o_R_x, limit_broken=0;			/* Right-boundary, left-boundary and old right-boundary positions */
	int R_part, L_part, size=1;					/* Right-boundary particle, displacement for broken clusters and size*/
	int broken=0;						/* Broken cluster variable */
	int incluster=1, incluster_L=1;				/* Inside cluster variables */
	int c=1;							/* Cluster counter */
	float prv, nxt;						/* Incluster limits */
	for (i=0; i<N; i++) {
	  if (positions[i]<R_x) {R_x=positions[i]; R_part=i;}	/* Find most-left particle - starting right boundary */
	  if (positions[i]>L_x) {L_x=positions[i]; L_part=i;}	/* Find most-right particle - starting left boundary */
	}
	if (DEBUG == 6) {fprintf(logfile, "Left-most particle: %d @ position %.5f\nRight-most particle: %d @ position %.5f\n",R_part+1,R_x,L_part+1,L_x);}
	prv=R_x-cluster_cutoff;
	if (prv<0 && L_x>prv+L) {						/* If there most-right particle is less than cluster_cutoff away (w/ PBCs) from most-left particle... */
	  broken=c;							  /* cluster is broken */
	  clusters[R_part]=c;
	  clusters[L_part]=c;					  /* both R_part and L_part are in the cluster (1st one) */
	  size=2;							  /* size of this cluster is now 2 */
	}
	if (DEBUG == 6) {fprintf(logfile, "Broken? (1=yes - 0=no): %d\n",broken);}
	if (broken==1) {						/* If first cluster is broken... */ /* Compute its size and limits to the left */
	  while (incluster_L==1) {					  /* as long as we're still in this cluster to the left... */
	    prv=L_x-cluster_cutoff;					    /* maximum distance for a particle to be in the same cluster */
	    if (DEBUG == 6) {fprintf(logfile, "    New limit for cluster inclusion: %.5f\n", prv);}
	    for (i=0; i<N; i++) {
	      incluster_L=0;					    /* by default assume it's the end of the cluster -> we exit this cluster */
	      if (positions[i]>prv && positions[i]<L_x && i!=L_part) {  /* if there is a particle to the left of the current boundary less than cluster_cutoff away... */
	        L_x=positions[i];					      /* update boundary to new particle */
	        L_part=i;
	        size++;						      /* increase size */
	        clusters[i]=c;					      /* tag cluster */
	        incluster_L=1;					      /* we're still inside the cluster! */
	        if (DEBUG == 6) {fprintf(logfile, "    Particle %d @ position %.5f -> cluster #%d of size %d\n", i+1, positions[i], c, size);}
	        break;						      /* exit the particle loop and start over with new boundary */
	      }
	    }
	    if (DEBUG == 6 && incluster==0) {fprintf(logfile, "  End of cluster reached\n");}
	  }
	  if (DEBUG == 6) {fprintf(logfile, "Left-boundary of broken cluster is particle %d @ position %.5f\n", L_part+1, L_x);}
	}
	if (DEBUG == 6) {fprintf(logfile, "\nMeasure clusters towards the right...\n");}
	while (R_x<L_x) {						/* While new right-boundary is to the left of the left-boundary (we haven't probed the whole system yet) */
	  while (incluster==1) {					  /* as long as we're still in this cluster */
	    float nxt=R_x+cluster_cutoff;				    /* maximum distance for a particle to be in the same cluster */
	    if (DEBUG == 6) {fprintf(logfile, "  New limit for cluster inclusion: %.5f\n", nxt);}
	    for (i=0; i<N; i++) {
	      incluster=0;						    /* by default assume it's the end of the cluster -> we exit this cluster */
	      if (positions[i]<nxt && positions[i]>R_x && i!=R_part) {  /* if there is a particle to the right of the current boundary less than cluster_cutoff away... */
	        R_x=positions[i];					      /* update boundary to new particle */
	        size++;						      /* increase size */
	        if (size==2) {clusters[R_part]=c;}			      /* tag cluster of previous boundary if this is the second particle in the cluster */
	        R_part=i;
	        clusters[i]=c;					      /* tag cluster */
	        incluster=1;					      /* we're still inside the cluster! */
	        if (DEBUG == 6) {fprintf(logfile, "    Particle %d @ position %.5f -> cluster #%d of size %d\n", i+1, positions[i], c, size);}
	        break;						      /* exit the particle loop and start over with new boundary */
	      }
	    }
	    if (DEBUG == 6 && incluster==0) {fprintf(logfile, "  End of cluster reached\n");}
	  }
	  if (c==1 && broken==1)  {limit_broken=R_x;}		  /* save the right-boundary of broken cluster */
	  if (DEBUG == 6 && broken==1) {fprintf(logfile, "Right-limit of the broken cluster: %.5f\n",limit_broken);}
	  if (size>1) {c++;}					  /* change cluster because we've reached the end */
	  size=0;
	  o_R_x=R_x;
	  R_x=L_x;							  /* reset R_x to the left-most position that is still to the right of old R_x */
	  for (i=0; i<N; i++) {					  /* Find most-left particle - new cluster right boundary */
	    if (positions[i]<R_x && positions[i]>o_R_x && i!=R_part) {
	      R_x=positions[i]; R_part=i; incluster=1; size=1;
	    }
	  }
	  if (DEBUG == 6) {fprintf(logfile, "\nNew cluster #%d starting at position %.5f with size %d\n", c, R_x, size);}
	}

	/* save # of clusters */
	if ((step+1)%(int)ceil((Tmax/(measures*dt))) == 0 && step>0) 
		{
			fprintf(nc, "%f	%d\n", (step+1)*dt, c-1);
		}			

	/*DEBUG*/	if (DEBUG == 6) 
				{
	      			fprintf(logfile, "\n");
	      			for (i=0; i<N; i++) {fprintf(logfile, "%d ",clusters[i]);}
	      			fprintf(logfile, "\n");
	      			fprintf(logfile, "Number of clusters: %d\n",c-1);
	    		}





	/* CLUSTERS MVMNT */

	if (DEBUG >= 1) {fprintf(logfile, "\n\n======> Cluster movement");}

	if (CMOB==1) {						/* If moving clusters are chosen do collective dynamics */

	  /* CHANGE POSITION OF CLUSTERS */
	  if (DEBUG==6) {fprintf(logfile, "\nCOLLECTIVE DYNAMICS - t=%f\n",(step+1)*dt);}
	  int dir, size;						/* Cluster direction and size variables */
	  float clusterprob;					/* Cluster jumping probability variable */
	  float clRB, clLB;						/* Cluster boundaries */
	  float d=limit_broken+1;					/* Displacement for broken cluster */
	  int p;
	  for (i=1; i<c; i++) 
	  {
	    dir=0, size=0;						/* Reset direction and size for new cluster */
	    if (DEBUG==6) {fprintf(logfile, "\nCluster: %d -- Zeroth boundaries: %.3f-%.3f",i,clLB,clRB);}
	    for (p=0; p<N; p++) 
	    {					/* Compute the direction, size and boundaries of the cluster */
	      if (clusters[p]==i) 
	      {
	        float posit=positions[p];				/* Position */
	        if (i==broken) 
	        {
	          posit=posit-d;					/* Displacement of broken cluster */
	          if (posit<0) {posit=posit+L;}			/* PBCs */
	        }
	        if (size==0) {clRB=posit; clLB=posit;}		/* Reset boundaries for new cluster */
	        else 
	        {
	          if (posit>clRB) {clRB=posit;}
	          if (posit<clLB) {clLB=posit;}
	        }
	        size++;						/* Increment in size */
	        dir=dir+directions[p];				/* Compute resulting direction */
	        if (DEBUG==6) 
	        {
	          if (i==broken) {fprintf(logfile, "\n  Position of particle %d (displaced cluster by %.3f): %.3f -- New boundaries: %.3f-%.3f",size,d,posit,clLB,clRB);}
	          else {fprintf(logfile, "\n  Position of particle %d: %.3f -- New boundaries: %.3f-%.3f",size,posit,clLB,clRB);}
	        }
	      }
	    }
	    if (i==broken) {clRB+=d; clLB+=d;}			/* Undo displacement for broken cluster */
	    if (clRB>=L) {clRB-=L;}					/* PBCs */
	    LBclusters[i-1]=clLB;					/* Save the boundaries of the cluster */
	    RBclusters[i-1]=clRB;
	    clusterprob=fabsf(dir)/(float)size;			/* Define the jumping probability for the cluster (cambiado abs() por fabsf())*/
	    int rclust=rand();
	    if (DEBUG==6) 
	    {
	      fprintf(logfile, "\nPosition: %.3f-%.3f\nSize: %d\nDirection: %d\nProbability: %.3f\nRandom number: %.3f\n",clLB,clRB,size,dir,clusterprob,rclust/(float)RAND_MAX);
	    }
	    if (rclust<=clusterprob*RAND_MAX) 
		    {			/* Jumping probability met */
		      float nclRB=clRB+1;					/* Limit for excluded volume to the right */
		      if (nclRB>=L) {nclRB-=L;}				/* PBCs */
		      float pclLB=clLB-1;					/* Limit for excluded volume to the left */
		      if (pclLB<0) {pclLB+=L;}				/* PBCs */
		      int allowed=1;
		      if (dir>0) 
			      {						/* If positive direction... */
			        if (DEBUG==6) {fprintf(logfile, "Right-moving cluster. New limit if jump: %.3f\n",nclRB);}
			        for (p=0; p<N; p++) 
			        {				/* Check movement of cluster is allowed */
			          if (clusters[p]==i) {continue;}
			          float dd=fabsf(positions[p]-nclRB);			/* Distance to other particles */
			          // if (dd<0) {dd=-dd;}				/* Absolute value (lo he puesto en la linea de arriba)*/
			          if (dd>L/2) {dd=L-dd;}				/* PBCs */
			          if (dd<sigma) 
			          {
			            allowed=0;
			            if (DEBUG==6) {fprintf(logfile, "  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
			            break;
			          }							/* Movement isn't allowed if particles are close to it */
			          if (DEBUG==6) {fprintf(logfile, "  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
			        }
			        if (DEBUG==6) {fprintf(logfile, "Allowed: %d\n",allowed);}
			        if (allowed==1) 
			        {					/* ...and movement is allowed */
			          for (p=0; p<N; p++) 
			          {
			            if (clusters[p]==i) 
			            {
				            /* Move cluster particles up */	
				            positions[p]=positions[p]+Fp*dt;
				            if (positions[p]>=L) {positions[p]=positions[p]-L;}/* PBCs */
				            }
			          }
			        }
			      }
		      else if (dir<0) 
			      {					/* If negative direction... */
			        if (DEBUG==6) {fprintf(logfile, "Left-moving cluster. New limit if jump: %.3f\n",pclLB);}
			        for (p=0; p<N; p++) 
			        {				/* Check movement of cluster is allowed */
			          if (clusters[p]==i) {continue;}
			          float dd=fabsf(positions[p]-pclLB);			/* Distance to other particles */
			          // if (dd<0) {dd=-dd;}				/* Absolute value (lo he puesto en la linea de arriba)*/
			          if (dd>L/2) {dd=L-dd;}				/* PBCs */
			          if (dd<sigma) 
			          {
			            allowed=0;
			            if (DEBUG==6) {fprintf(logfile, "  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
			            break;
			          }							/* Movement isn't allowed if particles are close to it */
			          if (DEBUG==6) {fprintf(logfile, "  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
			        }
			        if (DEBUG==6) {fprintf(logfile, "Allowed: %d\n",allowed);}
			        if (allowed==1) 
			        {					/* ...and movement is allowed */
			          for (p=0; p<N; p++) 
			          {
			            if (clusters[p]==i) 
			            {				  /* Move cluster particles down */
			              positions[p]=positions[p]-dt*Fp;
			              if (positions[p]<0) {positions[p]=positions[p]+L;}/* PBCs */
			            }
			          }
			        }
			      }
		      if (DEBUG==6) 
			      {
			        int ii;
			        fprintf(logfile, "State of the system:\n");
			        for (ii=0; ii<N; ii++) {fprintf(logfile, "%.3f ",positions[ii]);}
			        fprintf(logfile, "\n");
			        for (ii=0; ii<N; ii++) {fprintf(logfile, "%d ",directions[ii]);}	
			        fprintf(logfile, "\n");
			        for (ii=0; ii<N; ii++) {fprintf(logfile, "%d ",clusters[ii]);}
			        fprintf(logfile, "\n");
			      }
		    }

	    /* Only at selected time steps - MEASURING THE SYSTEM */
	    if ((step+1)%(int)ceil((Tmax/(measures*dt))) ==0 && step>0) 
		    {				
		      Nsize[size-1]++;					/* Histogram of the cluster sizes at this time */
		      if (DEBUG==6) {fprintf(logfile, "Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
		      if (dir<0) {Ndir[0]++;}				/* Histogram of the cluster directions at this time */
		      else if (dir==0) {Ndir[1]++;}
		      else if (dir>0) {Ndir[2]++;}
		    }
	  }

	}	/* End of CLUSTER MVMNT */





	// MEASUREMENT, SAVING & OUTPUT

	if (DEBUG >= 1) {fprintf(logfile, "\n\n\n==========> MEASUREMENT & OUTPUT\n");}

	/* ONLY MEASURE SIZES AND DIRECTIONS (si CMOB = 1 esta parte esta metida mas arriba). */

	if (CMOB==0) 
	{		
	  if (DEBUG==7) {fprintf(logfile, "\n\nMEASURING CLUSTERS - t=%d\n\n",step+1);}
	  int dir, size;						/* Cluster direction and size variables */
	  int p;
	  for (i=1; i<c; i++) 
	  {
	    dir=0, size=0;						/* Reset direction and size for new cluster */
	    for (p=0; p<N; p++) 
	    {					/* Compute the direction, size and boundaries of the cluster */
	      if (clusters[p]==i) 
	      {
	        size++;						/* Increment in size */
	        dir=dir+directions[p];				/* Compute resulting direction */
	      }
	    }
	    if (DEBUG==7) {fprintf(logfile, "\nCluster: %d -- Size: %d",i,size);}

	    /* Only at selected time steps - MEASURING THE SYSTEM */
	    if ((step+1)%(int)ceil((Tmax/(measures*dt))) ==0 && step>0) 
	    {				
	      Nsize[size-1]++;					/* Histogram of the cluster sizes at this time */
	      if (DEBUG>=2) {fprintf(logfile, "Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
	      if (dir<0) {Ndir[0]++;}				/* Histogram of the cluster directions at this time */
	      else if (dir==0) {Ndir[1]++;}
	      else if (dir>0) {Ndir[2]++;}
	    }
	  }
	}


	/* SPATIAL CORRELATIONS */

	int tag1, tag2, idist;
	float dist;

	/* Only at selected time steps */
	if ((step+1)%(int)ceil((Tmax/(measures*dt)))  == 0 && step>0) 
	{

	  	/* DEBUG */ if (DEBUG==8)
				  {
				    fprintf(logfile, "\n\n\nState of the system:\n");
				    for (tag1=0; tag1<N; tag1++) {fprintf(logfile, "%.3f ",positions[tag1]);}
				    fprintf(logfile, "\n");
				    fprintf(logfile, "\n\n\nCOMPUTING CORRELATIONS...\n\n");
				  }

		/* Reset correlations */
		for (tag1=0; tag1<L/2; tag1++) {C[tag1]=0;}

		/* DEBUG */ if (DEBUG==8) {fprintf(logfile, "Correlations reset\n\n");}

		/* Sum over particles */	
		for (tag1=0; tag1<N; tag1++) 
		{	

			/* Sum over all other particles */			
			for (tag2=tag1; tag2<N; tag2++) 
			{
				dist = positions[tag1] - positions[tag2];

				if (dist<0) {dist=-dist;}
				if (dist>(L/2)) {dist=L-dist;}

				idist = (int)ceil(dist);
				C[idist] += 1/(float)N;

				/* DEBUG */ if (DEBUG==8) 
					      {
					        fprintf(logfile, "Particles %d and %d\nDistance: %.3f\nInteger distance: %d\nC[%d]=%.3f\n\n",tag1+1,tag2+1,dist,idist,idist,C[idist]);
					      }
			}
	  }

	  /* Save correlations */
	  fprintf(corr, "%f	", (step+1)*dt);
	  for (idist=0; idist<L/2; idist++) {fprintf(corr, "%.10f	", C[idist]);}
	  fprintf(corr, "\n");
	  if (DEBUG==8) {fprintf(logfile, "\n\nFINISHED COMPUTING CORRELATIONS\n\n\n");}

	}

	/* Save total number of clusters and measure the Energy of the system */
	int sn;
	if ((step+1)%(int)ceil((Tmax/(measures*dt)))  == 0 && step>0)
	{
		/* DEBUG */ if (DEBUG==8) {fprintf(logfile, "\nComputing energy...\n");}

		E=0;							/* reset energy */
	  	for (sn=0; sn<N; sn++) 
	  	{
	    	E += 0.5*pow(trueDirections[sn],2)+0.5*V[sn];	

	    	/*DEBUG*/ if (DEBUG==8) {fprintf(logfile, "  Particle %d - v=%.3f - V=%.3f - E=%.3f\n",sn,trueDirections[sn],V[sn],0.5*pow(trueDirections[sn],2)+V[sn]);}
	  	}

	  	/* compute and normalize E */
	  	E=E/(float)N;									

	  	/*DEBUG*/ if (DEBUG==8) {fprintf(logfile, "\nEnergy = %.5f\n",E);}

	  	/* save E */
	  	fprintf(nrj, "%f	%f\n", (step+1)*dt, E);
	}


		//Save system snapshot/directions.
		/* snap initial times */

		int saveStep = (int)Tint/(dt*measures);
		if (saveStep < 1){saveStep = 1;}

		if ( step % saveStep == 0 && (float)step <= Tint/dt)
		{
			fprintf(snap, "%f ", (step+1)*dt);
			for (sn=0; sn<N; sn++) 
			{
	  			fprintf(snap, "%.10f %d ", positions[sn], directions[sn]);
			}
			fprintf(snap, "\n");
		}

		/* snap final times */

		else if ( step % saveStep == 0 && (float)step >= (Tmax - (float)Tint) / dt )
		{
			fprintf(snapbis, "%f ", (step+1)*dt);
			for (sn=0; sn<N; sn++) 
			{
	  			fprintf(snapbis, "%.10f %d ", positions[sn], directions[sn]);
			}
			fprintf(snapbis, "\n");
		}

		/* Print to terminal */
		if (DEBUG>1) {
		  fprintf(logfile, "\nFinal state of the system @ t=%d:\n",step+1);
		  for (i=0; i<N; i++) {fprintf(logfile, "%.3f ",positions[i]);}
		  fprintf(logfile, "\n");
		  for (i=0; i<N; i++) {fprintf(logfile, "%d ",directions[i]);}
		  fprintf(logfile, "\n");
		  for (i=0; i<N; i++) {fprintf(logfile, "%d ",clusters[i]);}
		  fprintf(logfile, "\n");
		  if ((step+1)%measures==0 && step>0) {
		    for (i=0; i<N; i++) {fprintf(logfile, "%d ",Nsize[i]);}
		    fprintf(logfile, "\n");
		  }
		}
		else if (DEBUG==0) 
		{
		  tt = clock() - t_0;
		  double seconds_taken=((double)tt)/CLOCKS_PER_SEC;
		  int minutes_taken=(int)floor(seconds_taken/(float)60);
		  seconds_taken-=minutes_taken*60;
		  int hours_taken=(int)floor(minutes_taken/(float)60);
		  minutes_taken-=hours_taken*60;
		  char message[96];
		  sprintf(message,"%.2f %% elapsed after %d hours %d minutes and %.2f seconds",(step+1)/(double)(Tmax/dt)*100,hours_taken,minutes_taken,seconds_taken);
		  if (step==0) 
		  {
		    printf("\nL=%d phi=%.3f alpha=%.3f epsilon=%.5f v=%.1f beta=%.3f D=%.3f CMOB=%d IS=%d dt=%.5f\n",L,fi,alpha,epsilon,Fp,beta,Dt,CMOB,INIT_STATE,dt);
		    printf("%s",message);
		  }
		  else if ( step < (ceil(Tmax/dt) - 1.f) )
		  {
		    int mm;
		    for (mm=0; mm<96; mm++) {printf("\b");}
		    printf("%s",message);
		  }
		  else {printf("\n");}
		}

	} /* END OF TIME LOOP */
	
	/* DEBUG */ if (DEBUG >= 1) {printf("\n\n==========> TIME LOOP (ending)\n\n\n");}

	fclose(logfile);
	




	/*SAVE CLUSTER SIZE DISTRIBUTION (NSIZE)*/
	sc=fopen(buffersc, "wb");
	scn=fopen(bufferscn, "wb");
	int norm=0;
	for (i=0; i<N; i++) {norm+=Nsize[i];}
	for (i=0; i<N; i++) 
	{
		if (Nsize[i]!=0) 
		{
			fprintf(sc, "%d	%.10f\n", i+1, Nsize[i]/(float)norm);
			fprintf(scn, "%f	%.10f\n", (i+1)/(float)N, Nsize[i]/(float)norm);
		}
	}

	fclose(sc);
	fclose(scn);

	/*SAVE CLUSTER DIRECTION DISTRIBUTION (NDIR)*/
	dc=fopen(bufferdc, "wb");
	for (i=0; i<3; i++) {if (Ndir[i]!=0) {fprintf(dc, "%d	%.10f\n", i-1, Ndir[i]/(float)(Ndir[0]+Ndir[1]+Ndir[2]));}}
	fclose(dc);

	/*SAVE PARTICLE DIRECTION DISTRIBUTION (NDIR_swim)*/
	pdc=fopen(bufferpdc, "wb");
	for (i=0; i<2; i++) {if (Ndir_swim[i]!=0) {fprintf(pdc, "%d	%.10f\n", (int)pow(-1,i), Ndir_swim[i]/(float)(Ndir_swim[0]+Ndir_swim[1]));}}
	fclose(pdc);

	/*SAVE PARTICLE RESULTING DIRECTION DISTRIBUTION (NDIR_res)*/
	rpdc=fopen(bufferrpdc, "wb");
	for (i=0; i<3; i++) {if (Ndir_res[i]!=0) {fprintf(rpdc, "%d	%.10f\n", i-1, Ndir_res[i]/(float)(Ndir_res[0]+Ndir_res[1]+Ndir_res[2]));}} 
	fclose(rpdc);

	fclose(nrj);
	fclose(nc);
	fclose(snap);	
	fclose(snapbis);	
	fclose(corr);

	epsilon=epsilon*24;

	/* DEBUG */ if (DEBUG >= 1) {printf("\n\n\n|||||||||||||||||||| SIMULATION END ||||||||||||||||||||\n\n\n");}

return 0;
}								/* END OF MAIN */
