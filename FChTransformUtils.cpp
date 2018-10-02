/*
  copyright (C) 2012 I. Irigaray, M. Rocamora

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "FChTransformUtils.h"
#include <math.h>

void cumtrapz(const double *x, const double *y, int N, double *accum)
/*Trapezoidal Integrator: 1/2(b-a)(F(a)+F(b))*/
{
    accum[0]=0.0;
    for (int i = 1; i < N; i++) {
        accum[i]=accum[i-1]+0.5*(x[i]-x[i-1])*(y[i]+y[i-1]);
    }
}


void interp1(const double *x1, const double *y1, int N1, const double *x2, double *y2, int N2){
/*1-D linear interpolation*/

    for (int i = 0; i < N2; i++) {
        /*Smaller or equal than the smallest, or larger or equal than the largest.*/
        if ( x2[i] <= x1[0] ) {
            y2[i] = y1[0];
        } else if ( x2[i] >= x1[N1-1] ) {
            y2[i] = y1[N1-1];
        } else {
            /*Search every value of x2 in x1*/
            int j = 1;
            int salir = 0;
            while ((j<N1)&&(!salir)) {			
                if ( x2[i] <= x1[j] ) {
                    y2[i] = y1[j-1] + ( ( y1[j] - y1[j-1] )*(x2[i] - x1[j-1] ) )/ ( x1[j] - x1[j-1] );
                    salir = 1;
                } // if
                j++;
            } // for
        } // else
    } // for

}

void interp1q(const double *y1, const int *x2_int, const double *x2_frac, double *y2, int N2){

    for(int i=0;i<N2;i++){
        y2[i] = y1[x2_int[i]]*(1.0-x2_frac[i])+y1[x2_int[i]+1]*x2_frac[i];
    } // for

}

void hanning_window(double *p_window, int n, bool normalize) {

    double accum=0;
    for (int i = 0; i < n; i++) {
        p_window[i] = 0.5*(1.0-cos(2*M_PI*(double)(i+1)/((double)n+1.0)));
        accum += p_window[i];
    }
    if (normalize) {
    	for (int i = 0; i < n; i++) { //window normalization
            p_window[i] = p_window[i]/accum;
    	}
    }

}
