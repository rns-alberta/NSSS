/*****************************************************************************
                    Polynomial interpolation function.
It calculates a value of central energy density of a star that has a baryonic mass
that is equal to that of a fixed star. 

The function uses the value of M_0 of the fixed star and the baryonic mass and 
central energy density of three stars to do the interpolation.
******************************************************************************/

#include <stdio.h>
#include <string.h> 
#include <math.h>
#include "equil_util.h"
#include "consts.h"
#include "nrutil.h"
#include "equil.h"
#include "interpol.h"



/***********Interpolation**********/
float polyinter(float M0, float e_c[], float M_0[]){
  float a, b, c, termx, const_term, y, sol1, sol2;
  float dif[3];
  int i;

  for(i=0;i<3;i++)
   dif[i] = M0 - M_0[i];

  // 3 point Interpolation
  a = dif[0] / ((e_c[0]-e_c[1])*(e_c[0]-e_c[2]));
  b = dif[1] / ((e_c[1]-e_c[0])*(e_c[1]-e_c[2]));
  c = dif[2] / ((e_c[2]-e_c[0])*(e_c[2]-e_c[1]));

  termx = (a*e_c[2] + a*e_c[1] + b*e_c[2] + b*e_c[0] + c*e_c[1] + c*e_c[0]);
  const_term = a*e_c[1]*e_c[2] + b*e_c[0]*e_c[2] + c*e_c[0]*e_c[1];

  // Interpolation second order equation
  y = (a+b+c)*M0*M0 - termx*M0 + const_term;

  //Solution of the second order equation
  sol1 = (termx + sqrt(termx*termx - (4*(a+b+c)*const_term))) / (2*(a+b+c));
  sol2 = (termx - sqrt(termx*termx - (4*(a+b+c)*const_term))) / (2*(a+b+c));

  //printf("\ne_c = %f \n", sol1);
  //printf("e_c = %f \n", sol2);

  return sol2;
}




