/* BROWNIAN DYNAMICS OF SELF-PROPELLED PARTICLES IN A ONE-DIMENSIONAL OFF-LATTICE SYSTEM */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* POTENTIAL FUNCTION */
float LJ(float epsilon, int sigma, float Vrange, float x) {
  float result;
  if (x<Vrange) {result=4*epsilon*(pow((sigma/x),12)-pow((sigma/x),6)-pow((1/(float)2.5),12)+pow((1/(float)2.5),6));}	//4*eps*((sigma/r)**12-(sigma/r)**6)
  else {result=0;}
  return result;
}

/* GRADIENT OF THE POTENTIAL FUNCTION */
float LJprime(float epsilon, int sigma, float Vrange, float x) {
  float result;
  if (x<Vrange) {result=24*epsilon*(pow(x,6)-2)*pow(x,-13);}	//24*(x^6-2)/(x^13)
  else {result=0;}
  return result;
}

void main(int argc, char **argv){

  srand(time(NULL));					/* Initialize random number generator */

  /* PARAMETERS */
  int DEBUG; 						/* Debug variable: 0 = no info printed - 1 = some info printed - 2 = all info printed */
  int L;						/* Lattice size */
  float alpha, alpha_0, alpha_int, alpha_end;		/* Tumbling rate: initial value, interval, final value */
  float fi, fi_0, fi_int, fi_end;			/* Density: initial value, interval, final value */
  int N;						/* # of particles */
  int Tmax;						/* Total simulation time */
  int Tint;						/* Measuring interval */
  float Vrange;						/* Range of the LJ potential */
  float epsilon, epsilon_0, epsilon_int, epsilon_end;	/* LJ potential intensity parameter */
  int sigma;						/* LJ potential particle diameters parameter */
  float Fp, Fp_0, Fp_int, Fp_end;			/* Propulsion force magnitude */
  float beta, beta_0, beta_int, beta_end;		/* 1/KbT: initial value, interval, final value */
  float mu, mu_0, mu_int, mu_end;			/* Friction coefficient: initial value, interval, final value */
  float Dt;						/* Translational diffusion coefficient*/
  int CMOB, CMOB_0, CMOB_int, CMOB_end;			/* Cluster mobility selector: 1 = motile clusters - 0 = still clusters */
  int INIT_STATE;					/* Initial state of the system selector: 0=random 1=coarsened 2=gas */
  float cluster_cutoff;					/* Cutoff distance for considering particles belong to same cluster */
  int ierr;						/* Return value of fscanf */
  float dt;           /* Time step */ 

  /* READ PARAMETERS */
  FILE *pars;
  pars=fopen("input.dat","r");
  ierr=fscanf(pars, "%d\n", & DEBUG);
  ierr=fscanf(pars, "%d\n", & L);
  ierr=fscanf(pars, "%f - %f - %f\n", & alpha_0, & alpha_int, & alpha_end);
  ierr=fscanf(pars, "%f - %f - %f\n", & fi_0, & fi_int, & fi_end);
  ierr=fscanf(pars, "%d\n", & Tmax);
  ierr=fscanf(pars, "%d\n", & Tint);
  ierr=fscanf(pars, "%f\n", & Vrange);
  ierr=fscanf(pars, "%f - %f - %f\n", & epsilon_0, & epsilon_int, & epsilon_end);
  ierr=fscanf(pars, "%d\n", & sigma);
  ierr=fscanf(pars, "%f - %f - %f\n", & Fp_0, & Fp_int, & Fp_end);
  ierr=fscanf(pars, "%f - %f - %f\n", & mu_0, & mu_int, & mu_end);
  ierr=fscanf(pars, "%f - %f - %f\n", & beta_0, & beta_int, & beta_end);
  ierr=fscanf(pars, "%d - %d - %d\n", & CMOB_0, & CMOB_int, & CMOB_end);
  ierr=fscanf(pars, "%f\n", & cluster_cutoff);
  ierr=fscanf(pars, "%d\n", & INIT_STATE);
  ierr=fscanf(pars, "%d\n", & dt);
  fclose(pars);

alpha=alpha_0;
while (alpha<=alpha_end) {			/* START ALPHA LOOP */
fi=fi_0;
while (fi<=fi_end) {				/* START PHI LOOP */
Fp=Fp_0;
while (Fp<=Fp_end) {				/* START ACTIVITY LOOP */
beta=beta_0;
while (beta<=beta_end) {			/* START TEMPERATURE LOOP */
mu=mu_0;
while (mu<=mu_end) {				/* START MOBILITY LOOP */
CMOB=CMOB_0;
while (CMOB<=CMOB_end) {			/* START CLUSTER EXPLICIT COLLECTIVE DYNAMICS LOOP */
epsilon=epsilon_0;
while (epsilon<=epsilon_end) {			/* START POTENTIAL STRENGTH LOOP */
  epsilon=epsilon/24;

  /* VARIABLES */
  N=L*fi;
  Dt=mu/beta;
  float positions[N];				/* Particles' positions */
  int directions[N];				/* Particles' swimming directions (1=+x; 0=-x;) */
  float trueDirections[N];			/* Particles' resulting velocities */
  int i, j;					/* Counter */
  float V[N];					/* Potential on each particle*/
  float Vprime[N];				/* Gradient of the potential on each particle*/
  int step;					/* Time step counter */
  int clusters[N];				/* Particles' clusters */
  float LBclusters[N/2];			/* Vector with the left boundaries of the clusters */
  float RBclusters[N/2];			/* Vector with the right boundaries of the clusters */
  int Nsize[N];					/* Vector with the number of clusters of each size */
  int Ndir[3]={0,0,0};				/* Vector with the number of clusters in each direction : +, -, 0 */
  int Ndir_swim[2]={0,0};			/* Vector with the number of particles swimming in each direction : +, - */
  int Ndir_res[3]={0,0,0};			/* Vector with the number of particles moving in each direction : +, -, 0 */
  float C[L/2];					/* Correlations */
  float E=0;					/* Energy of the system (normalized by the number of particles) */
  clock_t tt;					/* Real time */
  clock_t t_0=clock();				/* Real starting time */

  /* FILE PARAMETERS */
  char output_dir[192];									/* Output directory */
  sprintf(output_dir, "../../output_1DsppBD");
  struct stat st = {0};
  if (stat(output_dir, &st) == -1) {
    mkdir(output_dir, 0777);
    //int status = system("mv ");
  }
  char newdir[192];									/* Directory for the new simulation */
  sprintf(newdir, "%s/sim_a%.3f_f%.3f_t%.10d_L%.5d_D%.3f_Fp%.2f_eps%.5f_CMOB%d_IS%d_tint%.5d_dt%.5f", output_dir, alpha, fi, Tmax, L, Dt, Fp, epsilon, CMOB, INIT_STATE,Tint,dt);	
  struct stat st_bis = {0};
  if (stat(newdir, &st_bis) == -1) {
    mkdir(newdir, 0777);
  }
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
  char bufferspars[192];								/* Saves the parameters of this simulation */
  sprintf(bufferspars, "%s/parameters.dat", newdir);
  pars=fopen(bufferspars, "wb");
  fprintf(pars, "Debug variable: %d (0 = no info printed - 1 = some info printed - 2 = all info printed)\n", DEBUG);
  fprintf(pars, "Lattice size: %d\n", L);
  fprintf(pars, "Tumbling rate: %f\n", alpha);
  fprintf(pars, "Particle density: %f\n", fi);
  fprintf(pars, "Number of particles: %d\n", (int)N);
  fprintf(pars, "Total simulation time: %d\n", Tmax);
  fprintf(pars, "Measuring interval: %d\n", Tint);
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
  fclose(pars);

  if (DEBUG>0) {
    printf("\nSize of the system: %d\n# of particles: %d\nTumbling rate: %.3f\nTotal simulation time: %d\nMeasuring interval: %d\n",L,N,alpha,(int)Tmax,(int)Tint);
  }

  /* INITIAL STATE OF THE SYSTEM */
  if (DEBUG==3) {
    getchar();
    printf("\n\n\nINITIALIZING THE SYSTEM");
    if (INIT_STATE==0) {printf(" - RANDOM\n");}
    else if (INIT_STATE==2) {printf(" - SINGLE CLUSTER\n");}
    else if (INIT_STATE==1) {printf(" - DILUTE GAS\n");}
  }

  if (INIT_STATE==0) {							/* RANDOM INITIAL STATE */
    for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
      float r1=floor(rand()/(float)RAND_MAX*L);				/* Random position (0 to L) */
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
      int r2=rand()%2;							/* Random direction (1=+x; 0=-x;) */
      int r3;
      if (r2==1) {r3=1;}
      else if (r2==0) {r2=-1; r3=-1;}
      directions[i]=r2;
      trueDirections[i]=r3;
      clusters[i]=0;
      Nsize[i]=0;
    }
  }
  else if (INIT_STATE==1) {						/* GAS INITIAL STATE */
    for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
      float interval=L/(float)N;
      positions[i]=i*interval;						/* Ordered particles in a dilute gas at constant distance */
      int r2=rand()%2;							/* Random direction (1=+x; 0=-x;) */
      int r3;
      if (r2==1) {r3=1;}
      else if (r2==0) {r2=-1; r3=-1;}
      directions[i]=r2;
      trueDirections[i]=r3;
      clusters[i]=0;
      Nsize[i]=0;
    }
  }
  else if (INIT_STATE==2) {						/* COARSENED INITIAL STATE */
    for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
      float start_pos=L/(float)2-N/(float)2;
      positions[i]=start_pos+i;						/* Ordered particles in a single cluster in the middle of the system */
      int r2=rand()%2;							/* Random direction (1=+x; 0=-x;) */
      int r3;
      if (r2==1) {r3=1;}
      else if (r2==0) {r2=-1; r3=-1;}
      directions[i]=r2;
      trueDirections[i]=r3;
      clusters[i]=0;
      Nsize[i]=0;
    }
  }

  if (DEBUG>0) {
    printf("\nIntial state of the system:\n");
    for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
    printf("\n");
    for (i=0; i<N; i++) {printf("%d ",directions[i]);}
    printf("\n");
  }

  /* TIME LOOP - DYNAMICS*/
  if (DEBUG==3) {getchar(); printf("\n\n\nSTARTING TIME LOOP\n\n\n");}
  for (step=0; step<Tmax; step++) {

    for (i=0; i<N; i++) {clusters[i]=0;}		/* Reset the particle-to-cluster tagger */

    /* INDIVIDUAL DYNAMICS */
    if (DEBUG==3) {getchar(); printf("\n\n\nINDIVIDUAL DYNAMICS - t=%d\n",step+1);}

    /* CHOOSE A RANDOM ORDER OF MOVEMENT */
    int orderofmovement[N];				/* Choice of an order for moving the particles */
    for (i=0; i<N; i++) {
      int par=rand()%N+1;				/* Random particle to move at turn i */
      orderofmovement[i]=par;				/* Save particle par in turn i (even if repetition happens) */
    }

    /* CHANGE EACH PARTICLE STATE IN THE ABOVE ORDER */
    int k;
    for (k=0; k<N; k++) {
      int ptm=orderofmovement[k]-1;			/* Particle to move now */
      if (DEBUG>=2) {printf("\nParticle: %d\n",ptm+1);}

      /* TUMBLE */
      int prob=rand();					/* Random number to compare to alpha */
      if (prob<=alpha*RAND_MAX) {
        int tmbl=rand();				/* Random number to redefine direction */
        if (tmbl<=0.5*RAND_MAX) {
          directions[ptm]=-directions[ptm];
        }
      }

      /* POTENTIAL */
      int v;
      float pos=positions[ptm];				/* Position on the array of the particle */
      V[ptm]=0;						/* Reset the potential on that particle */
      Vprime[ptm]=0;					/* Reset the gradient of the potential on that particle */
      if (epsilon!=0) {
        for (v=0; v<N; v++) {
          if (v==ptm) {continue;}
          float dist=fabsf(pos-positions[v]);   // Cambiado abs() por fabsf() ya que abs(-0.999)=0.
          int posneg=(pos-positions[v])/dist;
          if (dist>L/2) {dist=L-dist;}
          V[ptm]+=LJ(epsilon,sigma,Vrange,dist);
          Vprime[ptm]+=posneg*LJprime(epsilon,sigma,Vrange,dist);
        }
      }
      if (DEBUG==3) {printf("F: %.10f\n",-Vprime[ptm]); printf("V: %.10f\n",V[ptm]);}

      /* GAUSSIAN WHITE NOISE */
      float U1=rand()/(float)RAND_MAX;			/* Uniformly distributed RV */
      //if (DEBUG==3) {printf("U1: %.3f\n",U1);}
      float U2=rand()/(float)RAND_MAX;			/* Uniformly distributed RV */
      //if (DEBUG==3) {printf("U2: %.3f\n",U2);}
      float Z1=sqrt(-2*log(U1))*cos(2*M_PI*U2);		/* Normally distributed RV */
      //if (DEBUG==3) {printf("Z1: %.3f\n",Z1);}
      float Z2=sqrt(-2*log(U1))*sin(2*M_PI*U2);		/* Normally distributed RV */
      //if (DEBUG==3) {printf("Z2: %.3f\n",Z2);}
      float eta=Z1;
      if (DEBUG==3) {printf("eta: %.3f\n",eta);}

      /* NEW POSITION */
      float npos;					/* New position */
      int vel=directions[ptm];				/* Velocity of the particle: 1 site/time-step in the direction */
      /*npos=pos+beta*Dt*(Fp*vel-Vprime[ptm])+sqrt(2*Dt)*eta;*/
      npos=pos+beta*Dt*dt*(Fp*vel-Vprime[ptm])+sqrt(2*Dt*dt)*eta;
      if (DEBUG>=2) {
        printf("Old position: %.3f\n",pos);
        printf("Velocity: %d\n",vel);
        printf("Expected new position: %.3f\n",npos);
      }
      
      /* RESULTING DIRECTION */
      trueDirections[ptm]=npos-pos;
      if (DEBUG>=2) {printf("True velocity: %.3f\n",trueDirections[ptm]);}

      /* FIND IF MOVEMENT IS ALLOWED */
      if (npos<0) {npos=npos+L;}			/* PBCs */
      else if (npos>=L) {npos=npos-L;}
      if (DEBUG>=2) {printf("Expected new position (with PBCs): %.3f\n",npos);}
      int allowed=1;
      for (v=0; v<N; v++) {
        if (v==ptm) {continue;}
        float dist=fabsf(npos-positions[v]);
        // if (dist<0) {dist=-dist;}			/* Absolute value (lo he puesto en la linea de arriba). */
        if (dist>L/2) {dist=L-dist;}			/* PBCs */
        if (dist<sigma) {allowed=0; break;}
        if (DEBUG==3) {printf("  Distance to particle %d: %.3f - allowed: %d\n",v+1,dist,allowed);}
      }
      if (allowed==1) {positions[ptm]=npos;}
      else {trueDirections[ptm]=0;}
      if (DEBUG>=2) {printf("Real new position: %.3f - Real new velocity: %.3f\n",positions[ptm],trueDirections[ptm]);}

      /* Save the direction of the particle: swimming and resulting (if Tint) */
      if ((step+1)%Tint==0 && step>0) {
        if (directions[ptm]==-1) {Ndir_swim[0]++;}
        else {Ndir_swim[1]++;}
        if (trueDirections[ptm]<0) {Ndir_res[0]++;}
        else if (trueDirections[ptm]>0) {Ndir_res[2]++;}
        else {Ndir_res[1]++;}
      }

    }

    if (DEBUG>=2) {
      printf("\nState of the system after individual dynamics:\n");
      for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
      printf("\n");
      for (i=0; i<N; i++) {printf("%d ",directions[i]);}
      printf("\n");
    }

    /* COLLECTIVE DYNAMICS */

    /* EVALUATE CLUSTERS */
    if (DEBUG==3) {getchar(); printf("\n\n\nEVALUATING CLUSTERS - t=%d\n\n\n",step+1);}

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
    if (DEBUG==3) {printf("Left-most particle: %d @ position %.5f\nRight-most particle: %d @ position %.5f\n",R_part+1,R_x,L_part+1,L_x);}
    prv=R_x-cluster_cutoff;
    if (prv<0 && L_x>prv+L) {						/* If there most-right particle is less than cluster_cutoff away (w/ PBCs) from most-left particle... */
      broken=c;							  /* cluster is broken */
      clusters[R_part]=c;
      clusters[L_part]=c;					  /* both R_part and L_part are in the cluster (1st one) */
      size=2;							  /* size of this cluster is now 2 */
    }
    if (DEBUG==3) {printf("Broken? (1=yes - 0=no): %d\n",broken);}
    if (broken==1) {						/* If first cluster is broken... */ /* Compute its size and limits to the left */
      while (incluster_L==1) {					  /* as long as we're still in this cluster to the left... */
        prv=L_x-cluster_cutoff;					    /* maximum distance for a particle to be in the same cluster */
        if (DEBUG==3) {printf("    New limit for cluster inclusion: %.5f\n", prv);}
        for (i=0; i<N; i++) {
          incluster_L=0;					    /* by default assume it's the end of the cluster -> we exit this cluster */
          if (positions[i]>prv && positions[i]<L_x && i!=L_part) {  /* if there is a particle to the left of the current boundary less than cluster_cutoff away... */
            L_x=positions[i];					      /* update boundary to new particle */
            L_part=i;
            size++;						      /* increase size */
            clusters[i]=c;					      /* tag cluster */
            incluster_L=1;					      /* we're still inside the cluster! */
            if (DEBUG==3) {printf("    Particle %d @ position %.5f -> cluster #%d of size %d\n", i+1, positions[i], c, size);}
            break;						      /* exit the particle loop and start over with new boundary */
          }
        }
        if (DEBUG==3 && incluster==0) {printf("  End of cluster reached\n");}
      }
      if (DEBUG==3) {printf("Left-boundary of broken cluster is particle %d @ position %.5f\n", L_part+1, L_x);}
    }
    if (DEBUG==3) {printf("\nMeasure clusters towards the right...\n");}
    while (R_x<L_x) {						/* While new right-boundary is to the left of the left-boundary (we haven't probed the whole system yet) */
      while (incluster==1) {					  /* as long as we're still in this cluster */
        float nxt=R_x+cluster_cutoff;				    /* maximum distance for a particle to be in the same cluster */
        if (DEBUG==3) {printf("  New limit for cluster inclusion: %.5f\n", nxt);}
        for (i=0; i<N; i++) {
          incluster=0;						    /* by default assume it's the end of the cluster -> we exit this cluster */
          if (positions[i]<nxt && positions[i]>R_x && i!=R_part) {  /* if there is a particle to the right of the current boundary less than cluster_cutoff away... */
            R_x=positions[i];					      /* update boundary to new particle */
            size++;						      /* increase size */
            if (size==2) {clusters[R_part]=c;}			      /* tag cluster of previous boundary if this is the second particle in the cluster */
            R_part=i;
            clusters[i]=c;					      /* tag cluster */
            incluster=1;					      /* we're still inside the cluster! */
            if (DEBUG==3) {printf("    Particle %d @ position %.5f -> cluster #%d of size %d\n", i+1, positions[i], c, size);}
            break;						      /* exit the particle loop and start over with new boundary */
          }
        }
        if (DEBUG==3 && incluster==0) {printf("  End of cluster reached\n");}
      }
      if (c==1 && broken==1)  {limit_broken=R_x;}		  /* save the right-boundary of broken cluster */
      if (DEBUG==3 && broken==1) {printf("Right-limit of the broken cluster: %.5f\n",limit_broken);}
      if (size>1) {c++;}					  /* change cluster because we've reached the end */
      size=0;
      o_R_x=R_x;
      R_x=L_x;							  /* reset R_x to the left-most position that is still to the right of old R_x */
      for (i=0; i<N; i++) {					  /* Find most-left particle - new cluster right boundary */
        if (positions[i]<R_x && positions[i]>o_R_x && i!=R_part) {
          R_x=positions[i]; R_part=i; incluster=1; size=1;
        }
      }
      if (DEBUG==3) {printf("\nNew cluster #%d starting at position %.5f with size %d\n", c, R_x, size);}
    }
    if ((step+1)%(int)sqrt(Tint)==0 && step>0) {fprintf(nc, "%d	%d\n", step+1, c-1);}			/* save # of clusters */
    if (DEBUG>=2) {
      printf("\n");
      for (i=0; i<N; i++) {printf("%d ",clusters[i]);}
      printf("\n");
      printf("Number of clusters: %d\n",c-1);
    }

    if (CMOB==1 && Fp!=0) {						/* If moving clusters are chosen do collective dynamics */

      /* CHANGE POSITION OF CLUSTERS */
      if (DEBUG==3) {getchar(); printf("\n\n\nCOLLECTIVE DYNAMICS - t=%d\n\n\n",step+1);}
      int dir, size;						/* Cluster direction and size variables */
      float clusterprob;					/* Cluster jumping probability variable */
      float clRB, clLB;						/* Cluster boundaries */
      float d=limit_broken+1;					/* Displacement for broken cluster */
      int p;
      for (i=1; i<c; i++) {
        dir=0, size=0;						/* Reset direction and size for new cluster */
        if (DEBUG>=2) {printf("\nCluster: %d -- Zeroth boundaries: %.3f-%.3f",i,clLB,clRB);}
        for (p=0; p<N; p++) {					/* Compute the direction, size and boundaries of the cluster */
          if (clusters[p]==i) {
            float posit=positions[p];				/* Position */
            if (i==broken) {
              posit=posit-d;					/* Displacement of broken cluster */
              if (posit<0) {posit=posit+L;}			/* PBCs */
            }
            if (size==0) {clRB=posit; clLB=posit;}		/* Reset boundaries for new cluster */
            else {
              if (posit>clRB) {clRB=posit;}
              if (posit<clLB) {clLB=posit;}
            }
            size++;						/* Increment in size */
            dir=dir+directions[p];				/* Compute resulting direction */
            if (DEBUG==3) {
              if (i==broken) {printf("\n  Position of particle %d (displaced cluster by %.3f): %.3f -- New boundaries: %.3f-%.3f",size,d,posit,clLB,clRB);}
              else {printf("\n  Position of particle %d: %.3f -- New boundaries: %.3f-%.3f",size,posit,clLB,clRB);}
            }
          }
        }
        if (i==broken) {clRB+=d; clLB+=d;}			/* Undo displacement for broken cluster */
        if (clRB>=L) {clRB-=L;}					/* PBCs */
        LBclusters[i-1]=clLB;					/* Save the boundaries of the cluster */
        RBclusters[i-1]=clRB;
        clusterprob=fabsf(dir)/(float)size;			/* Define the jumping probability for the cluster (cambiado abs() por fabsf())*/
        int rclust=rand();
        if (DEBUG>=2) {
          printf("\nPosition: %.3f-%.3f\nSize: %d\nDirection: %d\nProbability: %.3f\nRandom number: %.3f\n",clLB,clRB,size,dir,clusterprob,rclust/(float)RAND_MAX);
        }
        if (rclust<=clusterprob*RAND_MAX) {			/* Jumping probability met */
          float nclRB=clRB+1;					/* Limit for excluded volume to the right */
          if (nclRB>=L) {nclRB-=L;}				/* PBCs */
          float pclLB=clLB-1;					/* Limit for excluded volume to the left */
          if (pclLB<0) {pclLB+=L;}				/* PBCs */
          int allowed=1;
          if (dir>0) {						/* If positive direction... */
            if (DEBUG==3) {printf("Right-moving cluster. New limit if jump: %.3f\n",nclRB);}
            for (p=0; p<N; p++) {				/* Check movement of cluster is allowed */
              if (clusters[p]==i) {continue;}
              float dd=fabsf(positions[p]-nclRB);			/* Distance to other particles */
              // if (dd<0) {dd=-dd;}				/* Absolute value (lo he puesto en la linea de arriba)*/
              if (dd>L/2) {dd=L-dd;}				/* PBCs */
              if (dd<sigma) {
                allowed=0;
                if (DEBUG==3) {printf("  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
                break;
              }							/* Movement isn't allowed if particles are close to it */
              if (DEBUG==3) {printf("  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
            }
            if (DEBUG==3) {printf("Allowed: %d\n",allowed);}
            if (allowed==1) {					/* ...and movement is allowed */
              for (p=0; p<N; p++) {
                if (clusters[p]==i) {				  /* Move cluster particles up */
                  /*positions[p]=positions[p]+1;*/
                  positions[p]=positions[p]+Fp*dt;
                  if (positions[p]>=L) {positions[p]=positions[p]-L;}/* PBCs */
                }
              }
            }
          }
          else if (dir<0) {					/* If negative direction... */
            if (DEBUG==3) {printf("Left-moving cluster. New limit if jump: %.3f\n",pclLB);}
            for (p=0; p<N; p++) {				/* Check movement of cluster is allowed */
              if (clusters[p]==i) {continue;}
              float dd=fabsf(positions[p]-pclLB);			/* Distance to other particles */
              // if (dd<0) {dd=-dd;}				/* Absolute value (lo he puesto en la linea de arriba)*/
              if (dd>L/2) {dd=L-dd;}				/* PBCs */
              if (dd<sigma) {
                allowed=0;
                if (DEBUG==3) {printf("  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
                break;
              }							/* Movement isn't allowed if particles are close to it */
              if (DEBUG==3) {printf("  Distance to particle %d: %.3f - allowed: %d\n",p+1,dd,allowed);}
            }
            if (DEBUG==3) {printf("Allowed: %d\n",allowed);}
            if (allowed==1) {					/* ...and movement is allowed */
              for (p=0; p<N; p++) {
                if (clusters[p]==i) {				  /* Move cluster particles down */
                  positions[p]=positions[p]-1;
                  if (positions[p]<0) {positions[p]=positions[p]+L;}/* PBCs */
                }
              }
            }
          }
          if (DEBUG>=2) {
            int ii;
            printf("State of the system:\n");
            for (ii=0; ii<N; ii++) {printf("%.3f ",positions[ii]);}
            printf("\n");
            for (ii=0; ii<N; ii++) {printf("%d ",directions[ii]);}	
            printf("\n");
            for (ii=0; ii<N; ii++) {printf("%d ",clusters[ii]);}
            printf("\n");
          }
        }
  
        /* Only at selected time steps - MEASURING THE SYSTEM */
        if ((step+1)%Tint==0 && step>0) {				
          Nsize[size-1]++;					/* Histogram of the cluster sizes at this time */
          if (DEBUG>=2) {printf("Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
          if (dir<0) {Ndir[0]++;}				/* Histogram of the cluster directions at this time */
          else if (dir==0) {Ndir[1]++;}
          else if (dir>0) {Ndir[2]++;}
        }
      }

    }								/* End of cluster dynamics */

    else {							/* ONLY MEASURE SIZES AND DIRECTIONS */
      if (DEBUG==3) {getchar(); printf("\n\n\nMEASURING CLUSTERS - t=%d\n\n\n",step+1);}
      int dir, size;						/* Cluster direction and size variables */
      int p;
      for (i=1; i<c; i++) {
        dir=0, size=0;						/* Reset direction and size for new cluster */
        for (p=0; p<N; p++) {					/* Compute the direction, size and boundaries of the cluster */
          if (clusters[p]==i) {
            size++;						/* Increment in size */
            dir=dir+directions[p];				/* Compute resulting direction */
          }
        }
        if (DEBUG>=2) {printf("\nCluster: %d -- Size: %d",i,size);}
  
        /* Only at selected time steps - MEASURING THE SYSTEM */
        if ((step+1)%Tint==0 && step>0) {				
          Nsize[size-1]++;					/* Histogram of the cluster sizes at this time */
          if (DEBUG>=2) {printf("Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
          if (dir<0) {Ndir[0]++;}				/* Histogram of the cluster directions at this time */
          else if (dir==0) {Ndir[1]++;}
          else if (dir>0) {Ndir[2]++;}
        }
      }
    }

    /* SPATIAL CORRELATIONS */  
    int tag1, tag2, idist;
    float dist;
    if ((step+1)%Tint==0 && step>0) {				/* Only at selected time steps */
      if (DEBUG==3) {
        printf("\n\n\nState of the system:\n");
        for (tag1=0; tag1<N; tag1++) {printf("%.3f ",positions[tag1]);}
        printf("\n");
        printf("\n\n\nCOMPUTING CORRELATIONS...\n\n");
      }
      for (tag1=0; tag1<L/2; tag1++) {C[tag1]=0;}		/* Reset correlations */
      if (DEBUG==3) {printf("Correlations reset\n\n");}
      for (tag1=0; tag1<N; tag1++) {				/* Sum over particles */
        for (tag2=tag1; tag2<N; tag2++) {			/* Sum over all other particles */
          dist=positions[tag1]-positions[tag2];
          if (dist<0) {dist=-dist;}
          if (dist>(L/2)) {dist=L-dist;}
          idist=(int)ceil(dist);
          C[idist]+=1/(float)N;
          if (DEBUG==3) {
            printf("Particles %d and %d\nDistance: %.3f\nInteger distance: %d\nC[%d]=%.3f\n\n",tag1+1,tag2+1,dist,idist,idist,C[idist]);
          }
        }
      }
      fprintf(corr, "%d	", step+1);				/* Save correlations */
      for (idist=0; idist<L/2; idist++) {fprintf(corr, "%.10f	", C[idist]);}
      fprintf(corr, "\n");
      if (DEBUG==3) {printf("\n\nFINISHED COMPUTING CORRELATIONS\n\n\n");}
    }

    /* Save total number of clusters and system snapshot and measure the Energy of the system */
    if ((step+1)%(int)sqrt(Tint)==0 && step>0) {
      int sn;
      if (DEBUG==3) {printf("\nComputing energy...\n");}
      E=0;							/* reset energy */
      for (sn=0; sn<N; sn++) {
        E+=0.5*pow(trueDirections[sn],2)+0.5*V[sn];
        if (DEBUG==3) {
          printf("  Particle %d - v=%.3f - V=%.3f - E=%.3f\n",sn,trueDirections[sn],V[sn],0.5*pow(trueDirections[sn],2)+V[sn]);
        }
      }
      E=E/(float)N;						/* compute and normalize E */
      if (DEBUG==3) {printf("\nEnergy = %.5f\n",E); getchar();}
      fprintf(nrj, "%d	%f\n", step+1, E);			/* save E */
      if (step+1<=Tint) {
        fprintf(snap, "%d ", step+1);				/* snap initial times */
        for (sn=0; sn<N; sn++) {
          fprintf(snap, "%.10f %d ", positions[sn], directions[sn]);
        }
        fprintf(snap, "\n");
      }
      else if (step+1>=Tmax-Tint) {
        fprintf(snapbis, "%d ", step+1);			/* snap final times */
        for (sn=0; sn<N; sn++) {
          fprintf(snapbis, "%.10f %d ", positions[sn], directions[sn]);
        }
        fprintf(snapbis, "\n");
      }
    }

    /* Print to terminal */
    if (DEBUG>0) {
      printf("\nState of the system @ t=%d:\n",step+1);
      for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
      printf("\n");
      for (i=0; i<N; i++) {printf("%d ",directions[i]);}
      printf("\n");
      for (i=0; i<N; i++) {printf("%d ",clusters[i]);}
      printf("\n");
      if ((step+1)%Tint==0 && step>0) {
        for (i=0; i<N; i++) {printf("%d ",Nsize[i]);}
        printf("\n");
      }
      getchar();
    }
    else if (DEBUG==0) {
      tt = clock() - t_0;
      double seconds_taken=((double)tt)/CLOCKS_PER_SEC;
      int minutes_taken=(int)floor(seconds_taken/(float)60);
      seconds_taken-=minutes_taken*60;
      int hours_taken=(int)floor(minutes_taken/(float)60);
      minutes_taken-=hours_taken*60;
      char message[96];
      sprintf(message,"%.2f %% elapsed after %d hours %d minutes and %.2f seconds",(step+1)/(double)Tmax*100,hours_taken,minutes_taken,seconds_taken);
      if (step==0) {
        printf("\nL=%d phi=%.3f alpha=%.3f epsilon=%.5f v=%.1f beta=%.3f D=%.3f CMOB=%d IS=%d\n",L,fi,alpha,epsilon,Fp,beta,Dt,CMOB,INIT_STATE);
        printf("%s",message);
      }
      else if (step<Tmax-1) {
        int mm;
        for (mm=0; mm<96; mm++) {printf("\b");}
        printf("%s",message);
      }
      else {printf("\n");}
    }

  }								/* END OF TIME LOOP */


  /*SAVE CLUSTER SIZE DISTRIBUTION (NSIZE)*/
  sc=fopen(buffersc, "wb");
  scn=fopen(bufferscn, "wb");
  int norm=0;
  for (i=0; i<N; i++) {norm+=Nsize[i];}
  for (i=0; i<N; i++) {if (Nsize[i]!=0) {
    fprintf(sc, "%d	%.10f\n", i+1, Nsize[i]/(float)norm);
    fprintf(scn, "%f	%.10f\n", (i+1)/(float)N, Nsize[i]/(float)norm);
  }}
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
epsilon+=epsilon_int;
}						/* END OF POTENTIAL STRENGTH LOOP */
CMOB+=CMOB_int;
}						/* END OF CLUSTER EXPLICIT COLLECTIVE DYNAMICS LOOP */
mu+=mu_int;
}						/* END OF MOBILITY LOOP */
beta+=beta_int;
}						/* END OF TEMPERATURE LOOP */
Fp+=Fp_int;
}						/* END OF ACTIVITY LOOP */
fi+=fi_int;
}						/* END OF PHI LOOP */
alpha+=alpha_int;
}						/* END OF ALPHA LOOP */

}								/* END OF MAIN */
