/************************************************************************** 
*                         sequences.c
* 
* Computes constant neutron star sequences with constant baryonic mass
* 
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "surface.h"
#include "stableorbit.h"
#include "interpol.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
      
  int i, ierr;
  double
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
   spin_freq=100,                  /* Spin Frequency */
    Gamma_P;                      /* Gamma for polytropic EOS */  
                
  int j;

  int a = 0, numseq=2;
  int value = 0;
  float e_c[4], M_0[4];
  float M0, Mstat, Rstat, energy_value, temp_energy, ratio_r = 1.0, ej;
  float maxmass, maxradius;
  float T, W;
  //double Kfreq, Kfreq_j;

  FILE *fpointer;
  fpointer = fopen("NS_data.txt", "a");

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */



  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly" or "quark"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;  

      case 'b':
	sscanf(argv[i+1],"%lf",&B);
	B *= 1.602e33*KSCALE;
	break;       

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'd':
	/* CHOOSE THE NAME OF THE OUTPUT DIRECTORY */
	sscanf(argv[i+1],"%s",data_dir);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	break;

      case 'l':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_max);
	if(strcmp(eos_type,"poly")!=0)
	  e_max *= C*C*KSCALE;
	break;

      case 'n':  
  /* CHOOSE THE NUMBER OF SEQUENCES 
     PRODUCED */
	sscanf(argv[i+1],"%d",&numseq);
	break;

     case 'm':  
  /* CHOOSE MAXIMUM MASS FOR THE EOS */
  sscanf(argv[i+1],"%f",&maxmass);
  break;

     case 's':
	/* CHOOSE THE SPIN FREQUENCY (HZ) */
	sscanf(argv[i+1],"%lf",&spin_freq);
	printf("spin=%g\n",spin_freq);
	break;

     case 'r':
  /* CHOOSE THE RATIO OF r_p and r_e */
  sscanf(argv[i+1],"%f",&ratio_r);
  break;

     case 'x':
  /* CHOOSE THE MAXIMUM RADIUS FOR THE EOS*/
  sscanf(argv[i+1],"%f",&maxradius);
  break;

     case 't':
  /* CHOOSE TO COMPUTE SEQUENCE UP TO 800 HZ
   BY INPUT 1*/
  sscanf(argv[i+1],"%d",&value);
  break;

     case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        fprintf(stderr,"     quark : simple quark model \n"); 
	fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -d directory output goes to \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	fprintf(stderr,"  -l largest central energy density \n");
	fprintf(stderr,"  -h this menu\n\n");
	exit(1);
	break;  
      }
    }


  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("%s,  MDIVxSDIV=%dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */

  ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,
		    &eos, &star);

  //printf("The star infrastructure has been set up! \n");

  e_center = e_min;
  temp_energy = e_center;

  printf("e_c \t Mass \t Mass_0\t StatM \t Radius\tR-ratio\t StatR \t Spin\t K freq\n");
  printf("e15 \t Msun \t Msun\t Msun\t km\t --  \t km \t Hz \t Hz \n");

  if(ratio_r == 1.0 && value == 0){
    printf("Computing star with r_ratio = 1 up to Kepler limit\n");
  while ( a < numseq  ){
    ratio_r = 1.0;
    temp_energy = e_center;
    printf("Energy center = %g \n",e_center);
  
    ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    if((star.Mass/MSUN) > maxmass){
      printf("The maximum mass has been reached\n");
      break;
    }

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;

   T = (0.5*star.ang_mom*star.Omega)/(C*C);
   W = star.Mp + T - star.Mass;
  
    printf("%g \t%.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
	      star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);

    M0 = star.Mass_0/MSUN;
 while(1){   
 printf("---------------------------------------------------------------------------\n");
  ratio_r = ratio_r - 0.01;

  for(j=0;j<3;j++){
   ej = temp_energy - 0.01*j;

   ierr = MakeSphere(&eos, &star, ej);
   //rns(ratio_r, ej, &eos, &star); 
   if(ratio_r < 0.7)
      rns(0.7, ej, &eos, &star);
  
   rns(ratio_r, ej, &eos, &star); 

   printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
        ej, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI));

   e_c[j] = ej;
   M_0[j] = star.Mass_0/MSUN;
   }
 //printf("-----------------------------------------------\n");
   //printf("M_0 = %g\n",M0);
   //printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
   //printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

   energy_value = polyinter(M0, e_c, M_0);

 printf("---------------------------------------------------------------------------\n");
    temp_energy = energy_value;
    ierr = MakeSphere(&eos, &star, energy_value);
    //rns(ratio_r, energy_value, &eos, &star); 

    if(ratio_r < 0.7)
      rns(0.7, ej, &eos, &star);
    
    rns(ratio_r, ej, &eos, &star); 

    printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    if( (star.Omega/(2.0*PI)) > (star.Omega_K/(2.0*PI))) break;
    if(isnan(star.Mass/MSUN)){
      printf("Mass is NAN\n");
      break;
    }

    T = (0.5*star.ang_mom*star.Omega)/(C*C);
    W = star.Mp + T - star.Mass;
    //printf("Mp = %g T = %g Mass = %g  J = %g W = %g T/W = %g \n", star.Mp, T, star.Mass, star.ang_mom, W, T/W);
    //printf("M0 = %g \t Mass_0 = %g\n", M0, star.Mass_0/MSUN);
    if((round(M0*100.0)/100.0) == (round(star.Mass_0/MSUN * 100.0)/100.0))
    fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);

   }
   e_center = e_center + 0.1; //0.1
   a = a + 1;
  }
 }
// Computing sequences from 0 Hz to 800 Hz
else if(value == 1){
  printf("Computing star with small values of Omega\n");
  while ( a < numseq  ){
    ratio_r = 1.0;
    temp_energy = e_center;
    printf("Energy center = %g \n",e_center);
  
    ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    if((star.Mass/MSUN) > maxmass){
      printf("The maximum mass has been reached\n");
      break;
    }

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;

   T = (0.5*star.ang_mom*star.Omega)/(C*C);
   W = star.Mp + T - star.Mass;
  
    printf("%g \t%.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);

    M0 = star.Mass_0/MSUN;
 while(1){   
 printf("---------------------------------------------------------------------------\n");
  ratio_r = ratio_r - 0.0025;

  for(j=0;j<3;j++){
   ej = temp_energy - 0.01*j; 

   ierr = MakeSphere(&eos, &star, ej);
   //rns(ratio_r, ej, &eos, &star); 
   if(ratio_r < 0.7)
      rns(0.7, ej, &eos, &star);
  
   rns(ratio_r, ej, &eos, &star); 

   printf("%g %.5f  %.5f  %.5f %.5f %.3f %.4f %.4f \n",
        ej, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI));

   e_c[j] = ej;
   M_0[j] = star.Mass_0/MSUN;
   }
 //printf("-----------------------------------------------\n");
   //printf("M_0 = %g\n",M0);
   //printf("energies:  %g, %g, %g\n",e_c[0], e_c[1], e_c[2]);
   //printf("M_0's: %g, %g, %g\n",M_0[0], M_0[1], M_0[2]);

   energy_value = polyinter(M0, e_c, M_0);

 printf("---------------------------------------------------------------------------\n");
    temp_energy = energy_value;
    ierr = MakeSphere(&eos, &star, energy_value);
    //rns(ratio_r, energy_value, &eos, &star); 

    if(ratio_r < 0.7)
      rns(0.7, ej, &eos, &star);
    
    rns(ratio_r, ej, &eos, &star); 

    printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    if( (star.Omega/(2.0*PI)) > 800) break;
    if(isnan(star.Mass/MSUN)){
      printf("Mass is NAN\n");
      break;
    }

    T = (0.5*star.ang_mom*star.Omega)/(C*C);
    W = star.Mp + T - star.Mass;
    //printf("Mp = %g T = %g Mass = %g  J = %g W = %g T/W = %g \n", star.Mp, T, star.Mass, star.ang_mom, W, T/W);
    //printf("M0 = %g \t Mass_0 = %g\n", M0, star.Mass_0/MSUN);
    if((round(M0*100.0)/100.0) == (round(star.Mass_0/MSUN * 100.0)/100.0))
    fprintf(fpointer, "%7g %7g %7g %6g %4g %8g %4g %6g %7g %6g  %g  %g  %g  %g  %g\n", 
        energy_value, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, maxmass, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI), star.ang_mom, T, W, maxradius, maxmass/maxradius);

   }
   e_center = e_center + 0.1; //0.1
   a = a + 1;
  }

}

 // Computing one star
else{
    printf("Computing one neutron star\n");
    ratio_r = 0.70;
    ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;
  
    printf("%g %4.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));

    ratio_r = 0.644;
    //ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;
  
    printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));
        

    ratio_r = 0.643;
    //ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;
  
    printf("%g %.5f  %.5f  %.5f %.5f %.3f %.5f %.3f %.5f\n",
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));
  
    ratio_r = 0.640;
    //ierr = MakeSphere(&eos, &star, e_center);
    rns(ratio_r, e_center, &eos, &star); 

    Mstat = star.Mass/MSUN;
    Rstat = star.R_e*1e-5;
  
    printf("%g %4.5f  %4.5f  %4.5f %.5f %.3f %.5f %.3f %.5f\n",
        star.e_center, star.Mass/MSUN, star.Mass_0/MSUN, Mstat, star.R_e*1e-5, ratio_r, Rstat, star.Omega/(2.0*PI), star.Omega_K/(2.0*PI));


}

  fclose(fpointer);
  return 0;
}









