#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

/* GRADIENT OF THE POTENTIAL FUNCTION */
float LJprime(float epsilon, int sigma, int x) {
  float result;
  result=24*epsilon*(pow(x,6)-2)*pow(x,-13);	//24*(x^6-2)/(x^13)
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
  int Vrange;					/* Range of the LJ potential */
  float epsilon;				/* LJ potential intensity parameter */
  int sigma;					/* LJ potential particle diameters parameter */
  float Fp;					/* Propulsion force magnitude */
  float beta;					/* 1/KbT */
  float Dt;					/* Translational diffusion coefficient */
  int CMOB;					/* Cluster mobility selector: 1 = motile clusters - 0 = still clusters */
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
  ierr=fscanf(pars, "%d\n", & Vrange);
  ierr=fscanf(pars, "%f\n", & epsilon);
  epsilon=epsilon/24;
  ierr=fscanf(pars, "%d\n", & sigma);
  ierr=fscanf(pars, "%f\n", & Fp);
  ierr=fscanf(pars, "%f\n", & beta);
  ierr=fscanf(pars, "%f\n", & Dt);
  ierr=fscanf(pars, "%d\n", & CMOB);
  fclose(pars);

  /* VARIABLES */
  int array[L];					/* System array: 0=free 1=occupied */
  int positions[N];				/* Particles' positions */
  int directions[N];				/* Particles' swimming directions (1=+x; 0=-x;) */
  int trueDirections[N];			/* Particles' resulting moving directions (1=+x; 0=Stopped; -1=-x;) */
  int i;					/* Counter */
  float Vprime[N];				/* Gradient of the potential on each particle*/
  int step;					/* Time step counter */
  int arraycluster[L];				/* Vector with the clusters on the above array */
  int clusters[N];				/* Particles' clusters */
  int Nsize[N];					/* Particles' clusters */
  int Ndir[3];					/* Vector with the number of clusters in each direction : +, -, 0 */
  int Ndir_swim[2];				/* Vector with the number of particles swimming in each direction : +, - */
  int Ndir_res[3];				/* Vector with the number of particles moving in each direction : +, -, 0 */
  int C[L/2];					/* Correlations */
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
  fprintf(pars, "Range of the LJ potential: %d\n", Vrange);
  fprintf(pars, "LJ potential intensity parameter (epsilon): %f\n", epsilon);
  fprintf(pars, "LJ potential particle diameters parameter (sigma): %d\n", sigma);
  fprintf(pars, "Propulsion force: %f\n", Fp);
  fprintf(pars, "Beta - Inverse of the temperature energy: %f\n", beta);
  fprintf(pars, "Translational diffusivity: %f\n", Dt);
  fprintf(pars, "Mobility of the clusters: %d\n", CMOB);
  fclose(pars);

  if (DEBUG>0) {
    printf("\nSize of the lattice: %d\n# of particles: %d\nTotal simulation time: %d\nMeasuring interval: %d\n",L,N,(int)Tmax,(int)Tint);
  }

  /* INITIAL STATE OF THE SYSTEM - RANDOM - ONLY ROOT PROCESS */
  for (i=0; i<L; i++) {array[i]=0;}				/* Initialize array */	
  for (i=0; i<N; i++) {						/* Particles are N, ordered from 1 to N */
    int r1=rand()%L+1;						/* Random position (1 to L) */
    if (array[r1-1]==0) {
      array[r1-1]=i+1;
      positions[i]=r1;
      int r2=rand()%2;						/* Random direction (1=+x; 0=-x;) */
      int r3;
      if (r2==1) {r3=1;}
      else if (r2==0) {r2=-1; r3=-1;}
      directions[i]=r2;
      trueDirections[i]=r3;
    }
    else {i=i-1;}
  }
  if (DEBUG>0) {
    printf("\nIntial state of the system:\n");
    for (i=0; i<L; i++) {printf("%d ",array[i]);}
    printf("\n");
    for (i=0; i<L; i++) {
      if (array[i]==0) {printf("- ");}
      else if (directions[array[i]-1]==-1) {printf("< ");}
      else if (directions[array[i]-1]==1) {printf("> ");}
    }
    printf("\n");
  }

  /* TIME LOOP - DYNAMICS*/
  for (step=0; step<Tmax; step++) {

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
      int pos=positions[ptm]-1;				/* Position on the array of the particle */
      Vprime[ptm]=0;					/* Reset the gradient of the potential on that particle */
      for (v=1; v<Vrange+1; v++) {
        int Vposplus=pos+v;
        int Vposminus=pos-v;
        if (Vposplus>=L) {Vposplus-=L;}			/* PBCs */
        if (Vposminus<0) {Vposminus+=L;}
        if (array[Vposplus]!=0) {Vprime[ptm]-=LJprime(epsilon, sigma, v);}
        if (array[Vposminus]!=0) {Vprime[ptm]+=LJprime(epsilon, sigma, v);}
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
      int npos;						/* New discretized position on the lattice */
      float npos_f;					/* New continuous position on the lattie */
      int vel=directions[ptm];				/* Velocity of the particle: 1 site/time-step in the direction */
      npos_f=pos+beta*Dt*(Fp*vel-Vprime[ptm])+sqrt(2*Dt)*eta;
      npos=round(npos_f);
      if (DEBUG==2) {
        printf("Old position: %d\n",pos+1);
        printf("Velocity: %d\n",vel);
        printf("Expected new position (continuous): %.3f\n",npos_f+1);
        printf("Expected new position: %d\n",npos+1);
      }
      
      /* RESULTING DIRECTION */ 
      if (npos_f<pos) {trueDirections[ptm]=-1;}
      else if (npos_f>pos) {trueDirections[ptm]=1;}
      else {trueDirections[ptm]=0;}
      if (DEBUG==2) {printf("True direction: %d\n",trueDirections[ptm]);}

     /* FIND IF MOVEMENT IS ALLOWED */
     if (npos<0) {npos=npos+L;}				/* PBCs */
     else if (npos>=L) {npos=npos-L;}
     if (DEBUG==2) {printf("Expected new position (with PBCs): %d\n",npos+1);}
     if (array[npos]==0) {array[pos]=0; array[npos]=ptm+1; positions[ptm]=npos+1;}
     if (DEBUG==2) {printf("Real new position: %d\n",positions[ptm]);}

     /* Save the direction of the particle: swimming and resulting (if Tint) */
     if ((step+1)%Tint==0 && step>0) {
       if (directions[ptm]==-1) {Ndir_swim[0]++;}
       else {Ndir_swim[1]++;}
       if (trueDirections[ptm]==-1) {Ndir_res[0]++;}
       else if (trueDirections[ptm]==0) {Ndir_res[1]++;}
       else if (trueDirections[ptm]==1) {Ndir_res[2]++;}
     }

    }

    if (DEBUG==2) {
      printf("\nState of the system after individual dynamics:\n");
      for (i=0; i<L; i++) {printf("%d ",array[i]);}
      printf("\n");
      for (i=0; i<L; i++) {
        if (array[i]==0) {printf("- ");}
        else if (directions[array[i]-1]==-1) {printf("< ");}
        else if (directions[array[i]-1]==1) {printf("> ");}
      }
      printf("\n");
    }

    /* COLLECTIVE DYNAMICS */

    if (CMOB==0) {continue;}					/* If still clusters are chosen skip collective dynamics */

    /* EVALUATE CLUSTERS */
    int d, p, pp, nn, j=0;
    for (i=0; i<L; i++) {if (array[i]==0) {d=i; break;}}	/* Displacement to avoid separating a cluster due to end of count */
    for (i=0; i<L; i++) {
      p=i+d;							/* Position */
      nn=p+1;							/* Next one */
      pp=p-1;							/* Previous one */
      if (p>=L) {p=p-L;}					/* PBCs */
      if (nn>=L) {nn=nn-L;}
      if (pp>=L) {pp=pp-L;}
      if (pp<0) {pp=pp+L;}
      if (array[p]!=0) {					/* Save cluster it is in for each particle */
        if (array[pp]!=0 || array[nn]!=0) {			/* Only if there is another particle before or after (cond. 4 cluster) */
          clusters[array[p]-1]=j+1;				/* Save cluster tag to the particle */
        }
      }
      else if (array[p]==0 && clusters[array[pp]-1]!=0) {j++;}	/* If not change cluster */
    }
    int Nclusters=j;
    if (array[p]!=0) {Nclusters++;}				/* If we finish on a cluster add 1 to the count (haven't j++ yet) */
    if (DEBUG==2) {
      for (i=0; i<L; i++) {
        if (array[i]==0) {printf("- ");}
        else {printf("%d ",clusters[array[i]-1]);}
      }
      printf("\n");
      printf("\nDisplacement to avoid separating a cluster due to end of count: %d\n",d);
      printf("Number of clusters: %d\n",Nclusters);
    }

    /* CHANGE POSITION OF CLUSTERS */
    int dir, size;						/* Cluster direction and size variables */
    float clusterprob;						/* Cluster jumping probability variable */
    int clRB, clLB;						/* Cluster boundaries */
    for (i=0; i<Nclusters; i++) {
      dir=0, size=0;						/* Reset direction and size for new cluster */
      clRB=0, clLB=L;						/* Reset boundaries for new cluster */
      if (DEBUG==2) {printf("\nCluster: %d -- Zeroth boundaries: %d-%d",i+1,clLB,clRB);}
      for (p=0; p<N; p++) {					/* Compute the direction, size and boundaries of the cluster */
        if (clusters[p]==i+1) {
          size++;
          dir=dir+directions[p];
          int posit=positions[p]-d;				/* Displacement to avoid cutting clusters due to PBCs */
          if (posit<=0) {posit+=L;}				/* PBCs */
          if (posit>clRB) {clRB=posit;}
          if (posit<clLB) {clLB=posit;}
          if (DEBUG==2) {printf("\n  Position of particle %d: %d+%d -- New boundaries: %d-%d",size,posit,d,clLB,clRB);}
        }
      }
      clRB+=d; clLB+=d;						/* Undo displacement */
      if (clRB>L) {clRB-=L;}					/* PBCs */
      clusterprob=abs(dir)/(float)size;				/* Define the jumping probability for the cluster */
      int rclust=rand();
      if (DEBUG==2) {
        printf("\nPosition: %d-%d\nSize: %d\nDirection: %d\nProbability: %.3f\nRandom number: %.3f\n",clLB,clRB,size,dir,clusterprob,rclust/(float)RAND_MAX);
      }
      if (rclust<=clusterprob*RAND_MAX) {			/* Jumping probability met */
        int nclRB=clRB+1;					/* Next site */
        if (nclRB>=L) {nclRB=nclRB-L;}				/* PBCs */
        int pclLB=clLB-1;					/* Previous site */
        if (pclLB<0) {pclLB=pclLB+L;}				/* PBCs */
        if (dir>0 && array[nclRB-1]==0) {				/* If next site is free and positive direction movement is allowed */
          array[clLB-1]=0;
          for (p=0; p<N; p++) {
            if (clusters[p]==i+1) {				/* Move cluster particles up */
              positions[p]=positions[p]+1;
              if (positions[p]>L) {positions[p]=positions[p]-L;}/* PBCs */
              array[positions[p]-1]=p+1;
            }
          }
        }
        else if (dir<0 && array[pclLB-1]==0) {			/* If previous site is free and negative direction movement is allowed */
          array[clRB-1]=0;
          for (p=0; p<N; p++) {
            if (clusters[p]==i+1) {				/* Move cluster particles down */
              positions[p]=positions[p]-1;
              if (positions[p]<=0) {positions[p]=positions[p]+L;}/* PBCs */
              array[positions[p]-1]=p+1;
            }
          }
        }
        if (DEBUG==2) {
          printf("State of the system:\n");
          for (int ii=0; ii<L; ii++) {printf("%d ",array[ii]);}
          printf("\n");
          for (i=0; i<L; i++) {
            if (array[i]==0) {printf("- ");}
            else if (directions[array[i]-1]==-1) {printf("< ");}
            else if (directions[array[i]-1]==1) {printf("> ");}
          }
          printf("\n");
        }
      }

      /* Only at selected time steps - MEASURING THE SYSTEM */
      if ((step+1)%Tint==0 && step>0) {				
        Nsize[size-1]++;					/* Histogram of the cluster sizes at this time */
        if (DEBUG==2) {printf("Cluster %d Size %d N[%d]=%d\n",i,size,size,Nsize[size-1]);}
        if (dir<0) {Ndir[0]++;}					/* Histogram of the cluster directions at this time */
        else if (dir==0) {Ndir[1]++;}
        else if (dir>0) {Ndir[2]++;}
      }
    }

    /* SPATIAL CORRELATIONS */  
    int tag, dist, rest_i, rest_j;
    if ((step+1)%Tint==0 && step>0) {				/* Only at selected time steps */
      for (dist=0; dist<L/2; dist++) {C[dist]=0;}		/* Reset correlations */
      for (tag=0; tag<N; tag++) {				/* Sum over particles */
        for (dist=1; dist<L/2; dist++) {			/* For all distances */
          rest_i=positions[tag]-1+dist;				/* Right position */
          rest_j=positions[tag]-1-dist;				/* Left position */
          if (rest_i>=L) {rest_i=rest_i-L;}			/* PBCs */
          if (rest_j<0) {rest_j=rest_j+L;}			/* PBCs */
          if (array[rest_i]!=0) {C[dist-1]++;}			/* Increment if occuppied */
          if (array[rest_j]!=0) {C[dist-1]++;}			/* Increment if occuppied */
        }
      }
      fprintf(corr, "%d	", step+1);				/* Save correlations */
      for (dist=0; dist<N/2; dist++) {fprintf(corr, "%d	", C[dist]);}
      fprintf(corr, "\n");
    }

    /* Save total number of clusters and system snapshot*/
    if ((step)%(int)100==0) {
      fprintf(nc, "%d	%d\n", step+1, Nclusters);		/* # of clusters */
      int sn;
      //fprintf(snap, "%d ", step+1);				/* snap */
      for (sn=0; sn<L; sn++) {
        if (array[sn]==0) {fprintf(snap, "%d ", 0);}
        else if (directions[array[sn]-1]==1) {fprintf(snap, "%d ", 2);}
        else if (directions[array[sn]-1]==-1) {fprintf(snap, "%d ", 1);}
      }
      fprintf(snap, "\n");
    }

    /* Print to terminal */
    if (DEBUG>0) {
      printf("\nState of the system @ t=%d:\n",step+1);
      for (i=0; i<L; i++) {printf("%d ",array[i]);}
      printf("\n");
      for (i=0; i<L; i++) {
        if (array[i]==0) {printf("- ");}
        else if (directions[array[i]-1]==-1) {printf("< ");}
        else if (directions[array[i]-1]==1) {printf("> ");}
      }
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
