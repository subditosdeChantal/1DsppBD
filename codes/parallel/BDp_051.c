/* BROWNIAN DYNAMICS OF SELF-PROPELLED PARTICLES IN A ONE-DIMENSIONAL OFF-LATTICE SYSTEM */
/* Solved saving cluserts duplicately (in clustarr) Changed visited in cluster evaluation: visited[particlesOrder[auxL]] to visited[particlesOrder[auxL]]. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdbool.h>
#include <signal.h>
#include <string.h>

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
	float L;						/* Lattice size */
	float alpha;		/* Tumbling rate: initial value, interval, final value */
	float fi;			/* Density: initial value, interval, final value */
	int N;						/* # of particles */
	float Tmax;						/* Total simulation time */
	float T0 = 0;						/* Total simulation time */
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
	scanf("Number of particles (N): %d\n", & N);
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
	printf("Number of particles (N): %d\n", N);
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

	L = N/fi;
	Dt = mu/beta;
	epsilon=epsilon/24;

	/* FILE PARAMETERS */

	// If it doesn't exist, create it.
	struct stat st = {0};
	if (stat(output_dir, & st) == -1) {mkdir(output_dir, 0777);}

	char newdir[192];									/* Directory for the new simulation */
	sprintf(newdir, "%s/sim_a%.6f_f%.3f_t%010.2f_N%.5d_D%05.2f_Fp%.2f_eps%07.4f_CMOB%d_IS%d_tint%06.1f_dt%7.5f", output_dir, alpha, fi, Tmax, N, Dt, Fp, epsilon, CMOB, INIT_STATE,Tint,dt);	

	// If it doesn't exist, create it.
	struct stat st_bis = {0};
	if (stat(newdir, & st_bis) == -1) {mkdir(newdir, 0777);}



	FILE *logfile;									/* Log file */
	char bufferlog[192];
	sprintf(bufferlog, "%s/out.log", newdir);
	logfile  = fopen(bufferlog, "w");

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
	
	FILE *freearr;										/* Saves free particles at 1/10 of the measure times */
	char bufferfreearr[192];
	sprintf(bufferfreearr, "%s/freearr.dat", newdir);
	freearr=fopen(bufferfreearr, "wb");
	
	FILE *clustarr;										/* Saves cluster arrays at 1/10 of the measure times */
	char bufferclustarr[192];
	sprintf(bufferclustarr, "%s/clustarr.dat", newdir);
	clustarr=fopen(bufferclustarr, "wb");
	
	FILE *clustdyn;										/* Saves cluster dynamics at initial times (like snapshots) */
	char bufferclustdyn[192];
	sprintf(bufferclustdyn, "%s/clustdyn_beg.dat", newdir);
	clustdyn=fopen(bufferclustdyn, "wb");
	
	FILE *clustdynbis;										/* Saves cluster dynamics at steady state (like snapshots) */
	char bufferclustdynbis[192];
	sprintf(bufferclustdynbis, "%s/clustdyn_ss.dat", newdir);
	clustdynbis=fopen(bufferclustdynbis, "wb");
	
	FILE *csdt;											/* Saves CSD over time (at measures) */
	char buffercsdt[192];
	sprintf(buffercsdt, "%s/csd_t.dat", newdir);
	csdt=fopen(buffercsdt, "wb");

	FILE *pars;
	char bufferspars[192];								/* Saves the parameters of this simulation */
	sprintf(bufferspars, "%s/parameters.dat", newdir);

	if (DEBUG>0) {
	printf("\nSize of the system: %f\n# of particles: %d\nTumbling rate: %.3f\nTotal simulation time: %f\nMeasuring interval: %f\n",L,N,alpha,Tmax,Tmax/measures );
	}

		/* DEBUG */ if (DEBUG >= 1) {fprintf(logfile,"\n\n\n|||||||||||||||||||| SIMULATION START ||||||||||||||||||||\n\n\n");}


	/* VARIABLES */
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
	int clusterArray[(int)floor(N/2)][N];
	float LBclusters[N/2];			  /* Vector with the left boundaries of the clusters */
	float RBclusters[N/2];			  /* Vector with the right boundaries of the clusters */
	int Nsize[N];					        /* Vector with the number of clusters of each size */
	int Ndir[3]={0,0,0};				  /* Vector with the number of clusters in each direction : +, -, 0 */
	int Ndir_swim[2]={0,0};			  /* Vector with the number of particles swimming in each direction : +, - */
	int Ndir_res[3]={0,0,0};			/* Vector with the number of particles moving in each direction : +, -, 0 */
	float C[(int)ceil(L/2)];					        /* Correlations */
	float E=0;					          /* Energy of the system (normalized by the number of particles) */
	clock_t tt;					          /* Real time */
	clock_t t_0=clock();				  /* Real starting time */

	float CoM[N]; // center of mass
	int sizes[N]; // size of cluster
	float CSD[N]; // cluster size distribution
	float probJump[N]; // jumping prob
	int dirclust[N];  // jumping direction
	int Nclust = 0;

	// Checks.
	if (sigma != 1 && Vrange != 2.5)
	{
		fprintf(logfile,"ERROR Bad r_potmin! Recompute force zero (r_potmin) for given values of sigma and cutoff.\n");
		exit(1);
	}

	if (fi > 1/(r_potmin))
	{
		fprintf(logfile,"ERROR Density too high! Particle overlap! Maximum density is: %f\n", 1/r_potmin);
		exit(1);
	}

	/* INITIAL STATE OF THE SYSTEM */

	/* DEBUG */ if (DEBUG>=1) 
					{
						fprintf(logfile,"\n\n\n==========> INITIALIZING THE SYSTEM\n");
					}
	/* DEBUG */	if (DEBUG==3) 
					{
						if (INIT_STATE==0) {fprintf(logfile," - RANDOM\n");}
						else if (INIT_STATE==2) {fprintf(logfile," - SINGLE CLUSTER\n");}
						else if (INIT_STATE==1) {fprintf(logfile," - DILUTE GAS\n");}
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

				if (DEBUG==3) {fprintf(logfile,"Position %d: %.3f\n",i+1,r1);}

			int busy=0;
			for (j=0; j<i; j++) {
				float ddist=r1-positions[j];
				if (ddist<0) {ddist=-ddist;}
				if (ddist>L/2) {ddist=L-ddist;}

				    if (DEBUG==3) {fprintf(logfile,"  Position %d: %.3f - Distance: %.3f\n",j+1,positions[j],ddist);}

				if (ddist<sigma) {busy=1; continue;}
			}

				if (DEBUG==3) {fprintf(logfile,"Busy: %d\n",busy);}

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

	/* READ FROM FILE */
	else if (INIT_STATE==4) 
		{
		FILE *frame_file;
		char bufferframe0[192];
		sprintf(bufferframe0, "frame0_a%.6f_f%.3f_N%.5d_D%05.2f_Fp%.2f_eps%07.4f_CMOB%d.dat", alpha, fi, N, Dt, Fp, epsilon, CMOB);
		frame_file = fopen(bufferframe0, "r");
		fscanf(frame_file, "%f ", &T0);
		for (int i = 0; i < N; i++) 
			{						
		  		fscanf(frame_file, "%f %d ", &positions[i], &directions[i]);
		  		trueDirections[i] = directions[i];
 			}
 		fclose(frame_file);
		}

	pars=fopen(bufferspars, "wb");
	fprintf(pars, "Debug variable: %d (0 = no info printed - 1 = some info printed - 2 = all info printed)\n", DEBUG);
	fprintf(pars, "Lattice size: %f\n", L);
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
	fprintf(pars, "Initial state of the system: %d (0 = random - 1 = gas - 2 = coarsened - 3 = equidistant - 4 = read from file)\n", INIT_STATE);
	fprintf(pars, "Time step: %f\n", dt);
	fprintf(pars, "Snap measuring interval: %f\n", Tint);
	fprintf(pars, "Initial time: %f\n", T0);
	fclose(pars);

	/* Initiaize cluster arrays */
	for (int i = 0; i < N; ++i)
		{
			clusters[i]=0;
			Nsize[i]=0;
		}

		/*DEBUG*/ 	if (DEBUG > 1)
						{
							fprintf(logfile,"\nIntial state of the system:\n");
							for (i=0; i<N; i++) {fprintf(logfile,"%.3f ",positions[i]);}
							fprintf(logfile,"\n");
							for (i=0; i<N; i++) {fprintf(logfile,"%d ",directions[i]);}
							fprintf(logfile,"\n");
						}


	/* Compute particle order */

	/* DEBUG */ if (DEBUG >= 1)
				{
					fprintf(logfile,"\n\n========> Compute particle order");
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
					fprintf(logfile,"\n\n========> Neighbour list");
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
						fprintf(logfile,"\nPositions: \n");
						for(int n = 0 ; n < N; n++ ) {   
						 fprintf(logfile,"%f ", positions[n]);
						}

						fprintf(logfile,"\nOrdered Positions: \n");
						for(int n = 0 ; n < N; n++ ) {   
						 fprintf(logfile,"%f ", orderedPos[n]);
						}

						fprintf(logfile,"\nThe particles order: \n");
						for(int n = 0 ; n < N; n++ ) {   
						 fprintf(logfile,"%d ", particlesOrder[n]);
						}

						fprintf(logfile,"\nThe ordered particles: \n");
						for(int n = 0 ; n < N; n++ ) {   
						 fprintf(logfile,"%d ", orderedParticles[n]);
						}
						fprintf(logfile,"\n");
						fprintf(logfile,"\nThe particles neighbours: \n");
						for(int n = 0 ; n < N; n++ ) {   
							fprintf(logfile,"Particle %d. Neighbours: %d and %d \n", n, neighbours[n][0],neighbours[n][1]);
							fprintf(logfile,"Position %f. Neighbours: %f and %f \n", positions[n], positions[neighbours[n][0]], positions[neighbours[n][1]]);
							fprintf(logfile,"\n");
						}
						fprintf(logfile,"\n");
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
			fprintf(logfile,"\n ERROR OUT OF THE BOX! \npositions[%d] = %f\n", ptm, positions[ptm]);
			// exit(1);
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
							fprintf(logfile, "\nERROR QUANTUM JUMP. n");
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

							// exit(1);
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
							fprintf(logfile,"\n\nERROR Bad potential.\n");
							fprintf(logfile,"LEFT\n");						
							fprintf(logfile,"ptm: %d\n", ptm); 
							fprintf(logfile,"ptmaux: %d\n", ptmaux); 
							fprintf(logfile,"pos: %f\n", pos);
							fprintf(logfile,"posLneigh[ptmaux]: %f\n", positions[neighbours[ptmaux][0]]);
							fprintf(logfile,"distaux: %f\n", distaux);
							fprintf(logfile,"Vrange: %f\n", Vrange);
							fprintf(logfile,"eps: %f\n", epsilon);
							fprintf(logfile,"sigma: %d\n", sigma);
							fprintf(logfile,"Vprime[%d -> %d]: %f\n", neighbours[ptmaux][0], ptm, LJprime(epsilon, sigma, Vrange, distaux));
							fprintf(logfile,"VprimeTOT[%d]: %f\n\n", ptm, Vprime[ptm]);
							fprintf(logfile,"Number of left beighbours: %d\n", numbLeftNeighs);
							fprintf(logfile,"1st neighbour\n");
							fprintf(logfile,"positionsL[%d]: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
							fprintf(logfile,"Rest of neighbours\n");
							fprintf(logfile,"positions[ptmaux = %d]: %f\n", ptmaux, positions[neighbours[ptmaux][0]]);
							fprintf(logfile,"\nState: \n");
		      				for (i=0; i<N; i++) {fprintf(logfile,"%.3f ", positions[i]);}
		      				fprintf(logfile,"\n\n");
							// exit(1);
						}
					if (ptm == neighbours[ptmaux][0])
						{
							fprintf(logfile, "ERROR SELF REFERENCE.\n");
							// exit(1);
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
							fprintf(logfile,"\n\nERROR Bad potential.\n");
							fprintf(logfile,"RIGHT\n");
							fprintf(logfile,"ptm: %d\n", ptm); 
							fprintf(logfile,"ptmaux: %d\n", ptmaux); 
							fprintf(logfile,"pos: %f\n", pos);
							fprintf(logfile,"posRneigh[ptmaux]: %f\n", positions[neighbours[ptmaux][1]]);
							fprintf(logfile,"distaux: %f\n", distaux);
							fprintf(logfile,"Vrange: %f\n", Vrange);
							fprintf(logfile,"eps: %f\n", epsilon);
							fprintf(logfile,"sigma: %d\n", sigma);
							fprintf(logfile,"Vprime[%d -> %d]: %f\n", neighbours[ptmaux][1], ptm, LJprime(epsilon, sigma, Vrange, distaux));
							fprintf(logfile,"VprimeTOT[%d]: %f\n\n", ptm, Vprime[ptm]);
							fprintf(logfile,"Number of right beighbours: %d\n", numbRightNeighs);
							fprintf(logfile,"1st neighbour\n");
							fprintf(logfile,"positionsR[%d]: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
							fprintf(logfile,"Rest of neighbours\n");
							fprintf(logfile,"positions[ptmaux = %d]: %f\n", ptmaux, positions[neighbours[ptmaux][1]]);
							fprintf(logfile,"\nState: \n");
		      				for (i=0; i<N; i++) {fprintf(logfile,"%.3f ", positions[i]);}
		      				fprintf(logfile,"\n\n");

							// exit(1);
						}
					if (ptm == neighbours[ptmaux][1])
						{
							fprintf(logfile, "ERROR SELF REFERENCE.\n");
							// exit(1);
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
		npos = nfmod(npos, L);
		if (npos == L){npos = 0.f;}

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
		if (npos >= L || npos < 0.f || positions[ptm] >= L || positions[ptm] < 0.f || distR < 1.f || distL < 1.f)
			{
				fprintf(logfile,"\n\nBEFORE\n");
				fprintf(logfile,"eta: %.3f\n",eta);
				fprintf(logfile,"Vprime: %f\n", Vprime[ptm]);
				fprintf(logfile,"positions[%d]: %f\n", ptm, positions[ptm]);
				fprintf(logfile,"directions: %d\n",vel);
				fprintf(logfile,"trueDirections: %f\n", trueDirections[ptm]);
				fprintf(logfile,"position[%d] L: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
				fprintf(logfile,"position[%d] R: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
				fprintf(logfile,"pos: %f\nnpos: %f\ndistJump: %f\ndistL: %f\ndistL-potmin: %f\n", pos, npos, distJump, distL, distL-r_potmin);
				fprintf(logfile,"distR: %f\ndistR-potmin: %f\n", distR, distR-r_potmin);
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
						fprintf(logfile,"ERROR OVERLAP\ndistL = %f < 1\n", distL);		    	
					}

		    	// If the new position is at the other side of the PBC frontier.
		    	if (pos - distJump < 0.f)
		    		{ distJump = L - distJump; }

		    	// UPDATE POSITION, If the jump is smaller than the distance to its neighbor.
		    	if (distJump <= distL - r_potmin)
		    		{ 
		    			if (DEBUG > 1){	fprintf(logfile,"1L. jump < distNeigh\n"); } 

		    			positions[ptm] = npos; 
		    		}

		    	// If the jump is bigger, the particle sticks onto its neighbor only
			    	// IF IT HAS SPACE ON THE RIGHT, if not, it moves to an equidistant distance from its neighs.
			    	// This occurs when distR < r_potmin and distR + distL < 2*r_potmin (as in the initial case) 
		    		// the particle will move increasing distR and if there is a particle on the left it will 
		    		// overlap, so we have to take care of this.		    	
		    	else if (distJump > distL - r_potmin && distL + distR >= 2*r_potmin)
		    		{
		    			if (DEBUG > 1){ fprintf(logfile,"2L. jump > distNeigh, space\n"); }

		    			// If r_potmin of the particle is at the other side of the PBC frontier.		    		
		    			if (positions[neighbours[ptm][0]] + r_potmin > L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + r_potmin - L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][0]] + r_potmin < L && positions[neighbours[ptm][0]] + r_potmin >= 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + r_potmin;
			    			}
		    			// This can occur.
		    			else if (positions[neighbours[ptm][0]] + r_potmin == L)
			    			{
			    				positions[ptm] = 0.f;
			    			}
		    			else
		    				{
		    					fprintf(logfile,"ERROR - 1L\n");
		    					// exit(1);
		    				}
		    		}

		    	// If the particle doesn't have space, it moves to a equidistant distance from its neighbours.
		    	else if (distJump > distL - r_potmin && distL + distR < 2*r_potmin)
			    	{
						if (DEBUG > 1){	fprintf(logfile,"2L. jump > distNeigh, no space\n"); }
		    			// If (distR+distL)/2 is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][0]] + (distR+distL)/2 > L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + (distR+distL)/2 - L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][0]] + (distR+distL)/2 < L && positions[neighbours[ptm][0]] + (distR+distL)/2 >= 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][0]] + (distR+distL)/2;
			    			}
		    			// This can occur.
		    			else if (positions[neighbours[ptm][0]] + (distR+distL)/2 == L)
			    			{
			    				positions[ptm] = 0.f;
			    			}
		    			else
		    				{
		    					fprintf(logfile,"ERROR - 2L\n");
		    					// exit(1);
		    				}			    		
			    	}
		    	else 
		    		{
		    			fprintf(logfile, "\nERROR - 3L\ndistJump: %f\ndistR: %f\ndistL: %f",distJump, distR, distL);
		    			// exit(1);
		    		}
		    }

	    // If moving to the right.
	    if (trueDirections[ptm] > 0.f) 
		    {
		    		/* DEBUG */ if ( DEBUG > 2) {fprintf(logfile, "\nMOVING RIGHT\n");}

		    	// If the right neighbour is at the other side of the PBC frontier.
		    	if (pos + distR > L)
		    		{ distR = L - distR; }

		    	if (distR < 1)
		    		{
						fprintf(logfile,"ERROR OVERLAP\ndistR = %f < 1\n", distR);		    	
					}

		    	// If the new position is at the other side of the PBC frontier.
		    	if (pos + distJump > L)
		    		{ distJump = L - distJump; }

		    	// UPDATE POSITION, If the jump is smaller than the distance to its neighbor.
		    	if (distJump <= distR - r_potmin)
		    		{ 
		    			if (DEBUG > 1){	fprintf(logfile,"1R. jump < distNeigh\n"); } 
		    			positions[ptm] = npos; 
		    		}

		    	// If the jump is bigger, the particle sticks onto its neighbor only
			    	// IF IT HAS SPACE ON THE LEFT, if not, it sits there.
			    	// This occurs when distR < r_potmin and distR + distL < 2*r_potmin (as in the initial case) 
		    		// the particle will move increasing distR and if there is a particle on the left it will 
		    		// overlap, so we have to take care of this.
		    	else if (distJump > distR - r_potmin && distL + distR >= 2*r_potmin)
		    		{
		    			if (DEBUG > 1){ fprintf(logfile,"2R. jump > distNeigh, space\n"); }
		    			// If r_potmin of the particle is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][1]] - r_potmin < 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - r_potmin + L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][1]] - r_potmin >= 0.f && positions[neighbours[ptm][1]] - r_potmin < L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - r_potmin;
			    			}
		    			// This shouldn't occur.
		    			else if (positions[neighbours[ptm][1]] - r_potmin == L)
			    			{
			    				positions[ptm] = 0.f;
			    				fprintf(logfile,"ERROR - OJO!\n");
			    			}
		    			else
		    				{
		    					fprintf(logfile,"ERROR - 1R\n");
		    					// exit(1);
		    				}		    			

		    		}
		    	
		    	// If the particle doesn't have space, it moves to a equidistant distance from its neighbours.
		    	else if (distJump > distR - r_potmin && distL + distR < 2*r_potmin)
			    	{
						if (DEBUG > 1){	fprintf(logfile,"2R. jump > distNeigh, no space\n"); }

		    			// If (distR+distL)/2 is at the other side of the PBC frontier.
		    			if (positions[neighbours[ptm][1]] - (distR+distL)/2 < 0.f)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - (distR+distL)/2 + L;
			    			}
		    			// If it is in the bulk.
		    			else if (positions[neighbours[ptm][1]] - (distR+distL)/2 >= 0.f && positions[neighbours[ptm][1]] - (distR+distL)/2 < L)
			    			{
			    				positions[ptm] = positions[neighbours[ptm][1]] - (distR+distL)/2;
			    			}
		    			// This shouldn't occur.
		    			else if (positions[neighbours[ptm][1]] - (distR+distL)/2 == L)
			    			{
			    				positions[ptm] = 0.f;
			    				fprintf(logfile,"OJO!\n");
			    			}
						else
		    				{
		    					fprintf(logfile,"ERROR - 2R\n");
		    					// exit(1);
		    				}	
			    	}

		    	else 
		    		{
		    			fprintf(logfile, "\nERROR - 3R\ndistJump: %f\ndistR: %f\ndistL: %f",distJump, distR, distL);
		    			// exit(1);
		    		}
		    }

		// Due to floating point precision (e.g. 500 - 0.000013 = 500.000000) sometimes positions[ptm] = 500 occur, so we need to:
		if (positions[ptm] == L){ positions[ptm] = 0.f; }

	    // Check box again.
		if (positions[ptm] >= L || positions[ptm] < 0 || distR < 1.f || distL < 1.f)
			{
				fprintf(logfile,"\n\nERROR Out of BOX or OVERLAP.\n");
				fprintf(logfile,"\n");
				fprintf(logfile,"\n\npositions[%d]: %f\n", ptm, positions[ptm]);
				fprintf(logfile,"trueDirections: %f\n", trueDirections[ptm]);
				fprintf(logfile,"position[%d] L: %f\n", neighbours[ptm][0], positions[neighbours[ptm][0]]);
				fprintf(logfile,"position[%d] R: %f\n", neighbours[ptm][1], positions[neighbours[ptm][1]]);
				fprintf(logfile,"pos: %f\nnpos: %f\ndistJump: %f\ndistL: %f\ndistL-potmin: %f\n", pos, npos, distJump, distL, distL-r_potmin);
				fprintf(logfile,"distR: %f\ndistR-potmin: %f\n", distR, distR-r_potmin);
				// exit(1);
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

		/* DEBUG*/ if (DEBUG >= 1) {fprintf(logfile, "\n\n========> COLLECTIVE DYNAMICS \n");}


    /* EVALUATE CLUSTER (NEW) */

	/*SAVE CLUSTER ARRAY*/
	if ((step+1)%(int)ceil((Tmax/(measures*dt/10))) == 0 && step>0) 
	{
		fprintf(clustarr, "\nTime: %f\n", (step+1)*dt);
	}

	// Initializing variables.
	int i = 1;
	k = 0;
	int l = 0;
	bool visited[N];
    bool inclustL = false;
    bool inclustR = false;
    bool crossL = false;
    bool crossR = false;

	int auxL = 0;
	int auxR = 0;

	memset(visited, false, sizeof(visited));
	memset(CoM, 0, sizeof(CoM));
	memset(sizes, 0, sizeof(sizes)); 
	memset(CSD, 0, sizeof(CSD)); 
	memset(probJump, 0, sizeof(probJump));
	memset(dirclust, 0, sizeof(dirclust)); 
	memset(clusterArray, 0, sizeof(clusterArray));

	Nclust = 0;

	float distLneigh = 0;
	float distRneigh = 0;

	if ((step+1)%(int)ceil((Tmax/(measures*dt/10))) == 0 && step>0) 
	{
		fprintf(freearr, "%f ", step*dt);
	}

	for (int j = 0; j < N; ++j)
	{
			/* DEBUG*/ if (DEBUG >= 3) 
						{
							fprintf(logfile,"\n     -- For j = %d\n", j);
							fprintf(logfile,"        ID of clus i = %d\n", i);		
						}

		// If the particle has not been visited.
		if ( !visited[j] )
		{
				/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"          Entro en el if: visited[%d] = %d\n", j, visited[j]);}

			// Compute neighbor distances.
			distLneigh = fabsf(positions[j] - positions[neighbours[j][0]]);
			distRneigh = fabsf(positions[j] - positions[neighbours[j][1]]);

				/* DEBUG*/ if (DEBUG >= 3) 
							{
								fprintf(logfile,"          distLneigh = %f\n", distLneigh);
								fprintf(logfile,"          distRneigh = %f\n", distRneigh);
							}

			// If neighbors are at the other side of PBC frontier, fix distance.
			if ( positions[j] - distLneigh < 0.f)
				{ 
					distLneigh = L - distLneigh;
						/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"          PBC distLneigh = %f\n", distLneigh);}
				}

			if ( positions[j] + distRneigh > L)
				{ 
					distRneigh = L - distRneigh; 
						/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"          PBC distRneigh = %f\n", distRneigh);}				
				}

			// If the particle doesnt have neighbors, its free.
			if ( distLneigh >= cluster_cutoff && fabsf(positions[j] - positions[neighbours[j][1]]) >= cluster_cutoff )
			{
					/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"\n          1. Free Particle.\n");}			
				clusterArray[0][k] = j;
				visited[j] = true;			

				/*SAVE CLUSTER ARRAY*/
				if ((step+1)%(int)ceil((Tmax/(measures*dt/10))) == 0 && step>0) 
				{
					fprintf(freearr, "%d ", j);
				}
				k++;
				continue;
			}

				/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"\n          2. Particle in cluster.\n");}

	        inclustL = true;
	        inclustR = true;
	        auxL = j;
	        auxR = j;
			clusterArray[i][0] = j;
			visited[j] = true;			
	        CoM[i] += positions[j];
            dirclust[i] += directions[j];
	        l = 1;

	        	/* DEBUG*/ if (DEBUG >= 3) 
	        				{
	        					fprintf(logfile,"           auxL = auxR = j = %d\n", j);
        						fprintf(logfile,"           Visited Array:\n");
    							for (int iii = 0; iii < N; ++iii){ fprintf(logfile,"%d ", visited[iii]); }
    							fprintf(logfile,"\n           inclustL = %d, !visited[auxL = %d] = %d\n", inclustL, auxL, !visited[auxL]);
    						}

	        // Start exploring left cluster neighbors.
	        while(inclustL && !visited[neighbours[auxL][0]])
	        {
					/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"\n       -- 2.1 Left while.\n");}

	        	// Compute neighbor distances.
				distLneigh = fabsf(positions[auxL] - positions[neighbours[auxL][0]]);

					/* DEBUG*/ if (DEBUG >= 3) 
								{
									fprintf(logfile,"           positions[auxL = %d] = %f\n", auxL, positions[auxL]);
									fprintf(logfile,"           positions[Lneigh[auxL] = %d] = %f\n", neighbours[auxL][0], positions[neighbours[auxL][0]]);				
									fprintf(logfile,"           distLneigh = %f\n", distLneigh);				
									fprintf(logfile,"           distLneigh = %f\n", distLneigh);
								}

				// If neighbors are at the other side of PBC frontier, fix distance.
				if ( positions[auxL] - distLneigh < 0.f)
					{
						distLneigh = L - distLneigh;
						crossL = true;

							/* DEBUG*/ if (DEBUG >= 3) 
										{
											fprintf(logfile,"           Neighbor at the other side. crossL =  %d\n", crossL);											
											fprintf(logfile,"           PBC distLneigh = %f\n", distLneigh);
										}
					}

				// If the left neigbor is a cluster neighbor.
	        	if (distLneigh < cluster_cutoff)
	        	{
		        		/* DEBUG*/ if (DEBUG >= 3) 
		        					{
		        						fprintf(logfile,"           l = %d\n", l);
										fprintf(logfile,"           Left neigh is clust neigh (distLneigh = %f < %f).\n", distLneigh, cluster_cutoff);
									}

	        		// Move focus to him.
	        		auxL = neighbours[auxL][0];
	        		/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           auxL = %d\n", auxL);}

	        		// Add him to the cluster array.
					clusterArray[i][l] = auxL;
	        		/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           clusterArray[%d][%d] = %d\n", i, l, auxL);}					

					// Add him to the array of visited particles.
					visited[auxL] = true;
	        		/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           visited[%d] = %d\n", auxL, visited[auxL]);}										
					
					CoM[i] += positions[auxL];
	        		/* DEBUG*/ if (DEBUG >= 3) 
	        					{
	        						fprintf(logfile,"           positions[%d] = %f\n", auxL, positions[auxL]);										
	        						fprintf(logfile,"           CoM[%d] = %f\n", i, CoM[i]);
	        					}

                	if (crossL)
                		{			        					        		
                			CoM[i] -= L;

								/* DEBUG*/ if (DEBUG >= 3) 
                							{
                								fprintf(logfile,"           Cluster crosses left PBC. \n");                			
			        							fprintf(logfile,"           CoM[%d] = %f\n", i, CoM[i]);
			        						}
                		}

	                dirclust[i] += directions[auxL];

		        		/* DEBUG*/ if (DEBUG >= 3) 
	                				{
	                					fprintf(logfile,"           directions[%d] = %d\n", auxL, directions[auxL]);
	        							fprintf(logfile,"           dirclust[%d] = %d\n", i, dirclust[i]);
	        						}

					l++;

		        		/* DEBUG*/ if (DEBUG >= 3) 
									{
										fprintf(logfile,"           Visited Array:\n");
						        		for (int iii = 0; iii < N; ++iii)
						        			{
						        				fprintf(logfile,"%d ", visited[iii]);        				
						        			}

							        	fprintf(logfile,"\n           Cluster Array:\n");
							        	for (int ii = 0; ii < (int)floor(N/2); ++ii)
								        	{
								        		for (int iii = 0; iii < N; ++iii)
								        		{
								        			fprintf(logfile,"%d ", clusterArray[ii][iii]);        				
								        		}
							        			fprintf(logfile,"\n");        				
							    	    	}
							    	}

	        	} 
	        	else 
	        	{ 
	        			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           LEFT END OF CLUSTER\n\n");}

	        		// We reach the left end of the cluster.
	        		inclustL = false;

	        		// We put the left edge particle of the cluster in the 1st positions.
	        		clusterArray[i][l-1] = clusterArray[i][0];
	        		clusterArray[i][0] = auxL;

		        		/* DEBUG*/ if (DEBUG >= 3) 
	        						{
	        							fprintf(logfile,"           Cluster Array (swap edge particle):\n");
							        	for (int ii = 0; ii < (int)floor(N/2); ++ii)
								        	{
								        		for (int iii = 0; iii < N; ++iii)
								        		{
								        			fprintf(logfile,"%d ", clusterArray[ii][iii]);        				
								        		}
							        			fprintf(logfile,"\n");        				
							    	    	}
									}
	        	} 
	        }

    			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"\n           inclustR = %d, !visited[auxR = %d] = %d\n", inclustR, auxR, !visited[auxR]);}


	        // Start exploring right cluster neighbors.
	        while(inclustR && !visited[neighbours[auxR][1]])
	        {
					/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"\n       -- 2.2 Right while.\n");}

	        	// Compute neighbor distances.
				distRneigh = fabsf(positions[auxR] - positions[neighbours[auxR][1]]);

					/* DEBUG*/ if (DEBUG >= 3) 
								{
									fprintf(logfile,"           positions[auxR = %d] = %f\n", auxR, positions[auxR]);
									fprintf(logfile,"           positions[Rneigh[auxR] = %d] = %f\n", neighbours[auxR][1], positions[neighbours[auxR][1]]);								
									fprintf(logfile,"           distRneigh = %f\n", distRneigh);
								}

				// If neighbors are at the other side of PBC frontier, fix distance.
				if ( positions[auxR] + distRneigh > L)
					{
				 		distRneigh = L - distRneigh; 
						crossR = true;

							/* DEBUG*/ if (DEBUG >= 3)
										{
											fprintf(logfile,"           Neighbor at the other side. crossR =  %d\n", crossR);
											fprintf(logfile,"           PBC distRneigh = %f\n", distRneigh);
										}
				 	}

				// If the right neigbor is a cluster neighbor.
	        	if (distRneigh < cluster_cutoff)
	        	{
	        			/* DEBUG*/ if (DEBUG >= 3) 
	        						{
	        							fprintf(logfile,"           l = %d\n", l);
	        							fprintf(logfile,"           Right neigh is clust neigh (distRneigh = %f < %f).\n", distRneigh, cluster_cutoff);
	        						}

	        		// Move focus to him.	        		
	        		auxR = neighbours[auxR][1];
	        			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           auxR = %d\n", auxR);}

	        		// Add him to the cluster array.	        		
					clusterArray[i][l] = auxR;
	        			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           clusterArray[%d][%d] = %d\n", i, l, auxR);}

					// Add him to the array of visited particles.					
					visited[auxR] = true;
	        			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           visited[%d] = %d\n", auxR, visited[auxR]);}

					CoM[i] += positions[auxR];
                	if (crossR)
                		{	        		
                			CoM[i] += L;

			        			/* DEBUG*/ if (DEBUG >= 3) 
                							{
                								fprintf(logfile,"           Cluster crosses right PBC. \n");                			
			        							fprintf(logfile,"           CoM[%d] = %f\n", i, CoM[i]);
			        						}
                		}

	                dirclust[i] += directions[auxR];
	        		
		        		/* DEBUG*/ if (DEBUG >= 3) 
		        					{
		        						fprintf(logfile,"           directions[%d] = %d\n", auxR, directions[auxR]);
		        						fprintf(logfile,"           dirclust[%d] = %d\n", i, dirclust[i]);
		        					}

					l++;

			        	/* DEBUG*/ if (DEBUG >= 3) 
			        				{
			        					fprintf(logfile,"           Visited Array:\n");
						        		for (int iii = 0; iii < N; ++iii)
						        			{
						        				fprintf(logfile,"%d ", visited[iii]);        				
						        			}

							        	fprintf(logfile,"\n           Cluster Array:\n");
							        	for (int ii = 0; ii < (int)floor(N/2); ++ii)
								        	{
								        		for (int iii = 0; iii < N; ++iii)
								        		{
								        			fprintf(logfile,"%d ", clusterArray[ii][iii]);        				
								        		}
							        			fprintf(logfile,"\n");        				
							    	    	}
									}

	        	}
	        	else 
	        	{
	        			/* DEBUG*/ if (DEBUG >= 3) {fprintf(logfile,"           RIGHT END OF CLUSTER\n\n");}       		
	        		// We reach the right end of the cluster.	        		
	        		inclustR = false; 
	        	}       
	        }

			/*SAVE CLUSTER ARRAY*/
			if ((step+1)%(int)ceil((Tmax/(measures*dt/10))) == 0 && step>0) 
			{
				for (int jj = 0; jj < l; ++jj){ fprintf(clustarr, "%d ", clusterArray[i][jj]); }
				fprintf(clustarr, "\n");
			}	        

			sizes[i] = l;
			CoM[i] = nfmod(CoM[i]/sizes[i], L);
			CSD[sizes[i]-1] += 1;
			probJump[i] = abs(dirclust[i])/sizes[i];
			if (dirclust[i] != 0){dirclust[i] = dirclust[i]/abs(dirclust[i]);}
			crossL = false;
			crossR = false;

			if ((step+1)%(int)ceil((Tmax/(measures*dt))) ==0 && step>0) 
			{				
				Nsize[sizes[i]-1]++;								/* Histogram of the cluster sizes at this time */
				if (dirclust[i]<0) {Ndir[0]++;}					/* Histogram of the cluster directions at this time */
				else if (dirclust[i]==0) {Ndir[1]++;}
				else if (dirclust[i]>0) {Ndir[2]++;}
			}

	        i++;
		}
	}

	if ((step+1)%(int)ceil((Tmax/(measures*dt/10))) == 0 && step>0) 
	{
		fprintf(freearr, "\n");
	}

	Nclust = i-1;

	/*SAVE TOTAL NUMBER OF CLUSTERS*/
	if ((step+1)%(int)ceil((Tmax/(measures*dt))) == 0 && step>0) 
	{
		fprintf(nc, "%f	%d\n", (step+1)*dt, Nclust);
	}

	/*SAVE CSD*/
	if ((step+1)%(int)ceil((Tmax/(measures*dt))) == 0 && step>0) 
	{
		fprintf(csdt, "%f ", (step+1)*dt);
		int ii = 0;
		for (ii = 0; ii<N; ii++)
		{
			fprintf(csdt, "%f ", CSD[ii]);
		}
		fprintf(csdt, "\n");
	}


	// 

	/* CLUSTER MOVE (NEW) */

	float dist = 0.f;

	if (CMOB == 1)
	{

		float xR = 0.f;
		float xL = 0.f;
		float nxR = 0.f;
		float nxL = 0.f;
		int rclust = 0;

		for (int i = 1; i < Nclust + 1; ++i)
		{
			// Positions of the edge particles of the cluster.
			xR = positions[clusterArray[i][sizes[i]-1]];
			xL = positions[clusterArray[i][0]];

			// Positions of neighbors of the edge particles.
			nxR = positions[neighbours[clusterArray[i][sizes[i]-1]][1]];
			nxL = positions[neighbours[clusterArray[i][0]][0]];

			dist = 0.f;

			// If moving to the right.
			if (dirclust[i] > 0.f)
			{
	    		// Compute distance taking account of PBCs.
	    		dist = fabsf(nxR - xR);
				if ( xR + dist > L){ dist = L - dist; }
			}

			// Moving to the left.
			else if (dirclust[i] < 0.f)
			{
	    		// Compute distance taking account of PBCs.			
	    		dist = fabsf(xL - nxL);
				if ( xL - dist < 0.f){ dist = L - dist; }
			}

			// If the cluster jumps
			rclust=rand();
			if (rclust <= probJump[i]*RAND_MAX)
			{
				// And the distance to neighbor is larger than the jump distance.
				if (dist >= r_potmin + sigma)
				{
					// The cluster jumps normally.
					// Update center of mass.
					CoM[i] += dirclust[i]*dt*Fp;

					// PBCs
					if (CoM[i] >= L){ CoM[i] -= L; }
					else if (CoM[i] < 0.f){ CoM[i] += L; }
					if (CoM[i] == L){ CoM[i] += 0; }

					// Update positions of all particles in cluster.
					for (int k = 1; k < sizes[i]; ++k)
					{
						positions[clusterArray[i][k]] += dirclust[i]*dt*Fp;

						// PBCs
						if (positions[clusterArray[i][k]] >= L){positions[clusterArray[i][k]] -= L;}
	            		else if (positions[clusterArray[i][k]] < 0.f){positions[clusterArray[i][k]] += L;}
	            		if (positions[clusterArray[i][k]] == L){ positions[clusterArray[i][k]] = 0.f; }
	            	}
	            }
	            // If the distance to neighbor is smaller than the jump distance.
	            else if (dist < r_potmin + sigma)
				{
					// The cluster sticks onto its neighbor.
					CoM[i] += (dist - r_potmin)*dirclust[i]*dt;

					// PBCs
					if (CoM[i] >= L){ CoM[i] -= L; }
					else if (CoM[i] < 0){ CoM[i] += L; }
					if (CoM[i] == L){ CoM[i] += 0; }

					// Update positions of all particles in cluster.
					for (int k = 1; k < sizes[i]; ++k)
					{
	        			positions[clusterArray[i][k]] += (dist - r_potmin)*dirclust[i]*dt;

	        			//PBCs
	        			if (positions[clusterArray[i][k]] >= L){ positions[clusterArray[i][k]] -= L; }
	        			else if (positions[clusterArray[i][k]] < 0){ positions[clusterArray[i][k]] += L; }
	            		if (positions[clusterArray[i][k]] == L){ positions[clusterArray[i][k]] = 0.f; }
					}
				}
			}
		}
	}
	

	/* SPATIAL CORRELATIONS */

	int tag1, tag2, idist;
	dist = 0;

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
			// clusters
			fprintf(clustdyn, "%f ", (step+1)*dt);
			for (sn=1; sn<Nclust+1; sn++) 
			{
	  			fprintf(clustdyn, "%.10f %d %d ", CoM[sn], dirclust[sn], sizes[sn]);
			}
			fprintf(clustdyn, "\n");
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
			// clusters
			fprintf(clustdynbis, "%f ", (step+1)*dt);
			for (sn=1; sn<Nclust+1; sn++) 
			{
	  			fprintf(clustdynbis, "%.10f %d %d ", CoM[sn], dirclust[sn], sizes[sn]);
			}
			fprintf(clustdynbis, "\n");
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
		  if ((step+1)%measures==0 && step>0) 
		  {
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
		    printf("\nL=%f phi=%.3f alpha=%.3f epsilon=%.5f v=%.1f beta=%.3f D=%.3f CMOB=%d IS=%d dt=%.5f\n",L,fi,alpha,epsilon,Fp,beta,Dt,CMOB,INIT_STATE,dt);
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
	
	/* DEBUG */ 
	// if (DEBUG >= 1) {printf("\n\n==========> TIME LOOP (ending)\n\n\n");}

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
	fclose(freearr);
	fclose(clustarr);
	fclose(clustdyn);
	fclose(clustdynbis);
	fclose(csdt);

	epsilon=epsilon*24;

	/* DEBUG */ if (DEBUG >= 1) {fprintf(logfile,"\n\n\n|||||||||||||||||||| SIMULATION END ||||||||||||||||||||\n\n\n");}

return 0;
}								/* END OF MAIN */
