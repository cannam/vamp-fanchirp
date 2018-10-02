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

#ifndef FCHTRANSFORMUTILS_H
#define FCHTRANSFORMUTILS_H

#include <string.h>

class Utils
{
public:
    static void interp1(const double *x1,const double *y1, int N1, const double *x2, double *y2, int N2);

    static void interp1q(const double *y1, const int *x2_int, const double *x2_frac, double *y2, int N2){
        for(int i=0;i<N2;i++){
            y2[i] = y1[x2_int[i]]*(1.0-x2_frac[i])+y1[x2_int[i]+1]*x2_frac[i];
        } // for
    }

    static void cumtrapz(const double *x, const double *y, int N, double *accum);

    static void hanning_window(double *p_window, int n, bool normalize);
};

#endif

