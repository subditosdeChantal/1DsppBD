/* BROWNIAN DYNAMICS OF SELF-PROPELLED PARTICLES IN A ONE-DIMENSIONAL OFF-LATTICE SYSTEM */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* GRADIENT OF THE POTENTIAL FUNCTION */
float LJprime(float epsilon, int sigma, float Vrange, float x) {
  float result;
  if (x<Vrange) {result=24*epsilon*(pow(x,6)-2)*pow(x,-13);}	//24*(x^6-2)/(x^13)
  else {result=0;}
  return result;
}

void main(int argc, char **argv){

  srand(time(NULL));				/* Initialize random number generator */

  /* PARAMETERS */
  int DEBUG; 					/* Debug variable: 0 = no info printed - 1 = some info printed - 2 = all info printed */
  int L;					/* Lattice size */
  float alpha;					/* Tumbling rate */
  float fi;					/* Density */
  int N;					/* # of particles */
  int Tmax;					/* Total simulation time */
  int Tint;					/* Measuring interval */
  float Vrange;					/* Range of the LJ potential */
  float epsilon;				/* LJ potential intensity parameter */
  int sigma;					/* LJ potential particle diameters parameter */
  float Fp;					/* Propulsion force magnitude */
  float beta;					/* 1/KbT */
  float Dt;					/* Translational diffusion coefficient */
  int CMOB;					/* Cluster mobility selector: 1 = motile clusters - 0 = still clusters */
  float cluster_cutoff;				/* Cutoff distance for considering particles belong to same cluster */
  int ierr;					/* Return value of fscanf */

  /* READ PARAMETERS */
  FILE *pars;
  pars=fopen("input.dat","r");
  ierr=fscanf(pars, "%d\n", & DEBUG);
  ierr=fscanf(pars, "%d\n", & L);
  ierr=fscanf(pars, "%f\n", & alpha);
  ierr=fscanf(pars, "%f\n", & fi);
  N=L*fi;
  ierr=fscanf(pars, "%d\n", & Tmax);
  ierr=fscanf(pars, "%d\n", & Tint);
  ierr=fscanf(pars, "%f\n", & Vrange);
  ierr=fscanf(pars, "%f\n", & epsilon);
  epsilon=epsilon/24;
  ierr=fscanf(pars, "%d\n", & sigma);
  ierr=fscanf(pars, "%f\n", & Fp);
  ierr=fscanf(pars, "%f\n", & beta);
  ierr=fscanf(pars, "%f\n", & Dt);
  ierr=fscanf(pars, "%d\n", & CMOB);
  ierr=fscanf(pars, "%f\n", & cluster_cutoff);
  fclose(pars);

  /* VARIABLES */
  float positions[N];				/* Particles' positions */
  int directions[N];				/* Particles' swimming directions (1=+x; 0=-x;) */
  float trueDirections[N];			/* Particles' resulting velocities */
  int i;					/* Counter */
  float Vprime[N];				/* Gradient of the potential on each particle*/
  int step;					/* Time step counter */
  int clusters[N];				/* Particles' clusters */
  float LBclusters[N/2];			/* Vector with the left boundaries of the clusters */
  float RBclusters[N/2];			/* Vector with the right boundaries of the clusters */
  int Nsize[N];					/* Vector with the number of clusters of each size */
  int Ndir[3];					/* Vector with the number of clusters in each direction : +, -, 0 */
  int Ndir_swim[2];				/* Vector with the number of particles swimming in each direction : +, - */
  int Ndir_res[3];				/* Vector with the number of particles moving in each direction : +, -, 0 */
  float C[L/2];					/* Correlations */
  clock_t tt;					/* Real time */
  clock_t t_0=clock();				/* Real starting time */

  /* FILE PARAMETERS */
  char output_dir[192];									/* Output directory */
  sprintf(output_dir, "../../output_1DSPPdynamics");
  struct stat st = {0};
  if (stat(output_dir, &st) == -1) {
    mkdir(output_dir, 0777);
    //int status = system("mv ");
  }
  char newdir[192];									/* Directory for the new simulation */
  sprintf(newdir, "%s/sim_a%.3f_f%.3f_t%d_L%d_D%.2f_Fp%.2f_CMOB%d", output_dir, alpha, fi, Tmax, L, Dt, Fp, CMOB);	
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
  sprintf(bufferdc, "%s/dirdistr", newdir);
  FILE *pdc;										/* Saves the number of particles swimming in each direction overall */
  char bufferpdc[192];
  sprintf(bufferpdc, "%s/swim_dirdistr.dat", newdir);
  FILE *rpdc;										/* Saves the number of particles moving (result of LJ) in each direction overall */
  char bufferrpdc[192];
  sprintf(bufferrpdc, "%s/result_dirdistr.dat", newdir);
  FILE *nc;										/* Saves the total number of clusters */
  char buffernc[192];
  sprintf(buffernc, "%s/nclusters.dat", newdir);
  nc=fopen(buffernc, "wb");
  FILE *snap;										/* Saves snapshots of the system */
  char buffersnap[192];
  sprintf(buffersnap, "%s/system.dat", newdir);
  snap=fopen(buffersnap, "wb");
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
  fprintf(pars, "Translational diffusivity: %f\n", Dt);
  fprintf(pars, "Mobility of the clusters: %d\n", CMOB);
  fprintf(pars, "Cluster cutoff distance: %f\n", cluster_cutoff);
  fclose(pars);

  if (DEBUG>0) {
    printf("\nSize of the system: %d\n# of particles: %d\nTotal simulation time: %d\nMeasuring interval: %d\n",L,N,(int)Tmax,(int)Tint);
  }

  /* INITIAL STATE OF THE SYSTEM - RANDOM */
  for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
    float r1=rand()/(float)RAND_MAX*L;				/* Random position (0 to L) */
    positions[i]=r1;
    int r2=rand()%2;						/* Random direction (1=+x; 0=-x;) */
    int r3;
    if (r2==1) {r3=1;}
    else if (r2==0) {r2=-1; r3=-1;}
    directions[i]=r2;
    trueDirections[i]=r3;
    clusters[i]=0;
  }
  if (DEBUG>0) {
    printf("\nIntial state of the system:\n");
    for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
    printf("\n");
    for (i=0; i<N; i++) {printf("%d ",directions[i]);}
    printf("\n");
  }

  /* TIME LOOP - DYNAMICS*/
  for (step=0; step<Tmax; step++) {

    for (i=0; i<N; i++) {clusters[i]=0;}		/* Reset the particel-to-cluster tagger */

    /* INDIVIDUAL DYNAMICS */

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
      if (DEBUG==2) {printf("\nParticle: %d\n",ptm+1);}

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
      Vprime[ptm]=0;					/* Reset the gradient of the potential on that particle */
      for (v=0; v<N; v++) {
        if (v==ptm) {continue;}
        float dist=abs(pos-positions[v]);
        int posneg=(pos-positions[v])/dist;
        if (dist>L/2) {dist=L-dist;}
        Vprime[ptm]+=posneg*LJprime(epsilon,sigma,Vrange,dist);
      }
      if (DEBUG==2) {printf("F: %.10f\n",-Vprime[ptm]);}

      /* GAUSSIAN WHITE NOISE */
      float U1=rand()/(float)RAND_MAX;			/* Uniformly distributed RV */
      if (DEBUG==2) {printf("U1: %.3f\n",U1);}
      float U2=rand()/(float)RAND_MAX;			/* Uniformly distributed RV */
      if (DEBUG==2) {printf("U2: %.3f\n",U2);}
      float Z1=sqrt(-2*log(U1))*cos(2*M_PI*U2);		/* Normally distributed RV */
      if (DEBUG==2) {printf("Z1: %.3f\n",Z1);}
      float Z2=sqrt(-2*log(U1))*sin(2*M_PI*U2);		/* Normally distributed RV */
      if (DEBUG==2) {printf("Z2: %.3f\n",Z2);}
      float eta=Z1;
      if (DEBUG==2) {printf("eta: %.3f\n",eta);}

      /* NEW POSITION */
      float npos;					/* New position */
      int vel=directions[ptm];				/* Velocity of the particle: 1 site/time-step in the direction */
      npos=pos+beta*Dt*(Fp*vel-Vprime[ptm])+sqrt(2*Dt)*eta;
      if (DEBUG==2) {
        printf("Old position: %.3f\n",pos);
        printf("Velocity: %d\n",vel);
        printf("Expected new position: %.3f\n",npos);
      }
      
      /* RESULTING DIRECTION */
      trueDirections[ptm]=npos-pos;
      if (DEBUG==2) {printf("True velocity: %.3f\n",trueDirections[ptm]);}

      /* FIND IF MOVEMENT IS ALLOWED */
      if (npos<0) {npos=npos+L;}				/* PBCs */
      else if (npos>=L) {npos=npos-L;}
      if (DEBUG==2) {printf("Expected new position (with PBCs): %.3f\n",npos);}
      int allowed=1;
      for (v=0; v<N; v++) {
        if (v==ptm) {continue;}
        float dist=abs(npos-positions[v]);
        if (dist<sigma) {allowed=0;}
      }
      if (allowed==1) {positions[ptm]=npos;}
      if (DEBUG==2) {printf("Real new position: %.3f\n",positions[ptm]);}

      /* Save the direction of the particle: swimming and resulting (if Tint) */
      if ((step+1)%Tint==0 && step>0) {
        if (directions[ptm]==-1) {Ndir_swim[0]++;}
        else {Ndir_swim[1]++;}
        if (trueDirections[ptm]<0) {Ndir_res[0]++;}
        else if (trueDirections[ptm]>0) {Ndir_res[2]++;}
        else {Ndir_res[1]++;}
      }

    }

    if (DEBUG==2) {
      printf("\nState of the system after individual dynamics:\n");
      for (i=0; i<N; i++) {printf("%.3f ",positions[i]);}
      printf("\n");
      for (i=0; i<N; i++) {printf("%d ",directions[i]);}
      printf("\n");
    }

    /* COLLECTIVE DYNAMICS */

    /* EVALUATE CLUSTERS */

    int c=1;							/* Cluster counter */
    int j=0;							/* Particle counter */
    int broken=0;						/* Tag for broken cluster @ the boundary */
    float limit_broken=0;					/* Limit of broken cluster @ the boundary */
    for (i=0; i<N-1; i++) {					/* First particle loop */
      for (j=i+1; j<N; j++) {					/* Second particle loop */
        float dist=abs(positions[j]-positions[i]);		/* Distance */
        if (dist>L/2) {
          dist=L-dist;						/* PBCs */
          if (dist<cluster_cutoff) {broken=c;}			/* Two particles at the boundaries so cluster is broken */
        }
        if (dist<cluster_cutoff) {				/* Inside same cluster */
          if (clusters[i]==0) {
            if (clusters[j]==0) {				/* If none of the two particles have been assigned to a cluster... */
              clusters[i]=c;
              clusters[j]=c;					  /* Both particles assigned to current cluster count */
              c++;						  /* Increment cluster count */
            }
            else {clusters[i]=clusters[j];}			/* If particle i has no cluster but particle j does: assign cluster of j to i */
          }
          else {clusters[j]=clusters[i];}			/* If particle i has a cluster assign cluster of i to j */
        }
      } 
    }
    for (i=0; i<N; i++) {					/* Measure mass and limit of the broken cluster @ the boundary */
      if (clusters[i]==broken) {
        if (positions[i]<(limit_broken+cluster_cutoff)) {limit_broken=positions[i];}
      }
    }
    if (DEBUG==2) {
      for (i=0; i<N; i++) {printf("%d ",clusters[i]);}
      printf("\n");
      printf("Number of clusters: %d\n",c);
    }

    if (CMOB==1) {						/* If moving clusters are chosen do collective dynamics */

      /* CHANGE POSITION OF CLUSTERS */
      int dir, size;						/* Cluster direction and size variables */
      float clusterprob;					/* Cluster jumping probability variable */
      float clRB, clLB;						/* Cluster boundaries */
      float d=limit_broken+1;					/* Displacement for broken cluster */
      int p;
      for (i=1; i<c; i++) {
        dir=0, size=0;						/* Reset direction and size for new cluster */
        if (DEBUG==2) {printf("\nCluster: %d -- Zeroth boundaries: %.3f-%.3f",i+1,clLB,clRB);}
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
            if (DEBUG==2) {
              if (i==broken) {printf("\n  Position of particle %d (displaced cluster by %.3f): %.3f -- New boundaries: %.3f-%.3f",size,d,posit,clLB,clRB);}
              else {printf("\n  Position of particle %d: %.3f -- New boundaries: %.3f-%.3f",size,posit,clLB,clRB);}
            }
          }
        }
        if (i==broken) {clRB+=d; clLB+=d;}			/* Undo displacement for broken cluster */
        if (clRB>=L) {clRB-=L;}					/* PBCs */
        LBclusters[i-1]=clLB;					/* Save the boundaries of the cluster */
        RBclusters[i-1]=clRB;
        clusterprob=abs(dir)/(float)size;			/* Define the jumping probability for the cluster */
        int rclust=rand();
        if (DEBUG==2) {
          printf("\nPosition: %.3f-%.3f\nSize: %d\nDirection: %d\nProbability: %.3f\nRandom number: %.3f\n",clLB,clRB,size,dir,clusterprob,rclust/(float)RAND_MAX);
        }
        if (rclust<=clusterprob*RAND_MAX) {			/* Jumping probability met */
          float nclRB=clRB+1;					/* Limit for excluded volume to the right */
          if (nclRB>=L) {nclRB-=L;}				/* PBCs */
          float pclLB=clLB-1;					/* Limit for excluded volume to the left */
          if (pclLB<0) {pclLB+=L;}				/* PBCs */
          int allowed=1;
          if (dir>0) {						/* If positive direction... */
            for (p=1; p<i; p++) {
              if (LBclusters[p-1]<=nclRB) {allowed==0;}		/* Movement isn't allowed if another cluster is next to it */
            }
            if (allowed==1) {					/* ...and movement is allowed */
              for (p=0; p<N; p++) {
                if (clusters[p]==i) {				  /* Move cluster particles up */
                  positions[p]=positions[p]+1;
                  if (positions[p]>=L) {positions[p]=positions[p]-L;}/* PBCs */
                }
              }
            }
          }
          else if (dir<0) {					/* If negative direction... */
            for (p=1; p<i; p++) {
              if (RBclusters[p-1]>=pclLB) {allowed==0;}		/* Movement isn't allowed if another cluster is next to it */
            }
            if (allowed==1) {					/* ...and movement is allowed */
              for (p=0; p<N; p++) {
                if (clusters[p]==i) {				  /* Move cluster particles down */
                  positions[p]=positions[p]-1;
                  if (positions[p]<0) {positions[p]=positions[p]+L;}/* PBCs */
                }
              }
            }
          }
          if (DEBUG==2) {
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
          if (DEBUG==2) {printf("Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
          if (dir<0) {Ndir[0]++;}				/* Histogram of the cluster directions at this time */
          else if (dir==0) {Ndir[1]++;}
          else if (dir>0) {Ndir[2]++;}
        }
      }

    }								/* End of cluster dynamics */

    /* SPATIAL CORRELATIONS */  
    int tag, dist, rest_i, rest_j;
    float norm=N*(N-1)/(float)2;
    if ((step+1)%Tint==0 && step>0) {				/* Only at selected time steps */
      for (dist=0; dist<L/2; dist++) {C[dist]=0;}		/* Reset correlations */
      for (tag=0; tag<N-1; tag++) {				/* Sum over particles */
        for (dist=tag+1; dist<N; dist++) {			/* Sum over all other particles */
          float distance=abs(positions[tag]-positions[dist]);
          C[(int)distance]+=1/norm;
        }
      }
      fprintf(corr, "%d	", step+1);				/* Save correlations */
      for (dist=0; dist<L/2; dist++) {fprintf(corr, "%.10f	", C[dist]);}
      fprintf(corr, "\n");
    }

    /* Save total number of clusters and system snapshot */
    if ((step)%(int)100==0) {
      fprintf(nc, "%d	%d\n", step+1, c);			/* # of clusters */
      int sn;
      //fprintf(snap, "%d ", step+1);				/* snap */
      for (sn=0; sn<L; sn++) {fprintf(snap, "%.10f ", positions[sn]);}
      fprintf(snap, "\n");
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
      getchar();
    }
    else {
      tt = clock() - t_0;
      double seconds_taken=((double)tt)/CLOCKS_PER_SEC;
      int minutes_taken=(int)floor(seconds_taken/(float)60);
      seconds_taken-=minutes_taken*60;
      int hours_taken=(int)floor(minutes_taken/(float)60);
      minutes_taken-=hours_taken*60;
      char message[96];
      sprintf(message,"%.2f %% elapsed after %d hours %d minutes and %.2f seconds",(step+1)/(double)Tmax*100,hours_taken,minutes_taken,seconds_taken);
      if (step==0) {printf("%s",message);}
      else if (step<Tmax-1) {
        if (minutes_taken==1 && seconds_taken<0.001) {printf("\n%s",message);}
        else {
          int mm;
          for (mm=0; mm<96; mm++) {printf("\b");}
          printf("%s",message);
        }
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

  fclose(nc);
  fclose(snap);
  fclose(corr);


}								/* END OF MAIN */
