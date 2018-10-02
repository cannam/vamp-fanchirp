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

#include <string.h>

void interp1(const double *x1,const double *y1, int N1, const double *x2, double *y2, int N2);

void interp1q(const double *y1, const int *x2_int, const double *x2_frac, double *y2, int N2);

void cumtrapz(const double *x, const double *y, int N, double *accum);

void hanning_window(double *p_window, int n, bool normalize);
