#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include "skycalc.h"


// double jd = 2451545.;  // J2000 Julian Date

// Define globals
int update_on = 0;
double update_delta = 0.;
FILE *sclogfl = NULL;


double star_tzero, star_terr, star_period, star_perr;  /* for ephemeris calculations ... global */



double* calc_parallax(double ra, double dec, double jd){

  double coor_out[2];
  double dra, ddec;

  double ra_hour = ra/360.0 * 24.0;
  double epoch = 2015.0;

  double lat = 0.0;
  double longit = 0.0;

  double aberra, aberdec;

  parellipse(jd, ra_hour, dec, epoch, lat, longit, &dra, &ddec, &aberra, &aberdec);

  coor_out[0] = dra;
  coor_out[1] = ddec;

  // printf("HERE:%d %d\n", aberra, aberdec);

  return coor_out;
}

//
// static PyObject *
// parallax_get_shifts(PyObject *self, PyObject *args)
// {
//     const double ra, dec, jd;
//     double dra, ddec;
//
//     if (!PyArg_ParseTuple(args, "f", command))
//         return NULL;
//
//     jd_mod = jd + (double)i * 365.25 / 1000.0;
//
//     calc_parallax(ra, dec, jd_mod, &dra, &ddec);
//
//     //
//     // sts = system(command);
//     return Py_BuildValue("f", sts);
// }







//
// int main(){
//
//   double jd_mod;
//   double jd = 2451545.;  // J2000 Julian Date
//   // double ra = 210.0;
//   // double dec = 89.0;
//
//   double ra = 272.5126;
//   double dec = 66.540;
//
//   double dra, ddec;
//
//
//   for(int i=0; i<1000; i=i+1){
//
//     jd_mod = jd + (double)i * 365.25 / 1000.0;
//
//     calc_parallax(ra, dec, jd_mod, &dra, &ddec);
//     printf("dra, ddec = %f %f\n", dra, ddec);
//
//   }
//
//
//   return;
// }
