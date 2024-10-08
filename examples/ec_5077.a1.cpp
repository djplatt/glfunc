// Copyright Edgar Costa 2019
// See LICENSE file for license details.
/*
 * Make up a degree 2 L-function associated the Elliptic Curve with LMFDB label 5077.a1
 * https://www.lmfdb.org/EllipticCurve/Q/5077/a/1
 *
 * Change the following two lines, to modify the number of decimal digits printed, or to print raw format (not human friendly)
 */
#define DIGITS 20
#define RAW false
/*
 * with DIGITS = 20 and RAW = False, running this file should generate something
 * comparable to:
Order of vanishing = 3
Epsilon = (-1 + 4.1376012004041556323e-62j)  +/-  (2.53e-57, 3.23e-55j)
First non-zero Taylor coeff = 1.7318499001193006898 +/- 3.0255e-47
L(1.5) = (0.086781447504949672575 + 2.0382384424521477768e-58j)  +/-  (1.03e-39, 1.03e-39j)
L(2.5) = (0.51086394865482725426 + 1.9975501652723247176e-48j)  +/-  (6.79e-41, 6.79e-41j)
First 10 zeros
Zero 0 = 2.0524728584799397697 +/- 2.0881e-53
Zero 1 = 3.2624435559787574664 +/- 2.0881e-53
Zero 2 = 4.4705515133100979509 +/- 8.3524e-53
Zero 3 = 4.7544315159634058642 +/- 1.6705e-52
Zero 4 = 6.0119227529863951901 +/- 3.341e-52
Zero 5 = 6.6225046134077067814 +/- 5.3455e-51
Zero 6 = 7.3428149795396481469 +/- 2.1382e-50
Zero 7 = 7.7067946481132534446 +/- 2.1382e-50
Zero 8 = 8.4768019426235003774 +/- 2.1382e-50
Zero 9 = 9.3821789111719395491 +/- 3.4211e-49
Z-plot in [0, 10]:
0.00	0.00	                              o
0.50	0.22	                              |o
1.00	1.54	                              |----------o
1.50	2.98	                              |---------------------o
2.00	0.55	                              |---o
2.05	zero	                              Z
2.50	-4.96	------------------------------|
3.00	-4.10	o-----------------------------|
3.26	zero	                              Z
3.50	3.87	                              |----------------------------o
4.00	5.01	                              |------------------------------
4.47	zero	                              Z
4.50	-0.17	                             o|
4.75	zero	                              Z
5.00	2.41	                              |-----------------o
5.50	5.38	                              |------------------------------
6.00	0.12	                              o
6.01	zero	                              Z
6.50	-0.75	                         o----|
6.62	zero	                              Z
7.00	1.20	                              |-------o
7.34	zero	                              Z
7.50	-0.36	                            o-|
7.71	zero	                              Z
8.00	1.44	                              |---------o
8.48	zero	                              Z
8.50	-0.23	                             o|
9.00	-3.15	       o----------------------|
9.38	zero	                              Z
9.50	1.04	                              |------o
10.00	0.90	                              |-----o
10.20	zero	                              Z
 */
#define __STDC_FORMAT_MACROS
#include <chrono>
#include <cstdint>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
//#include <flint/fmpz.h>
//#include <flint/fmpzxx.h>
#include <flint/acb_poly.h>
#include "glfunc_internals.h"
//#include "examples_tools.h"
#include <cassert>

//using flint::fmpzxx;
using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;


// such dictionary can be obtained directly from the sidebar
// http://www.lmfdb.org/L/ModularForm/GL2/Q/holomorphic/23/1/b/a/
map<int64_t, vector<int64_t>> euler_factors  = {
  {2, {1, 2, 2}},
  {3, {1, 3, 3}},
  {5, {1, 4, 5}},
  {7, {1, 4, 7}},
  {521, {1, -18, 521}},
  {11, {1, 6, 11}},
  {13, {1, 4, 13}},
  {1039, {1, 5, 1039}},
  {17, {1, 4, 17}},
  {19, {1, 7, 19}},
  {23, {1, 6, 23}},
  {1049, {1, -31, 1049}},
  {1657, {1, -20, 1657}},
  {1051, {1, 43, 1051}},
  {29, {1, 6, 29}},
  {31, {1, 2, 31}},
  {1433, {1, 18, 1433}},
  {547, {1, -9, 547}},
  {37, {1, 0, 37}},
  {1063, {1, -8, 1063}},
  {41, {1, 0, 41}},
  {43, {1, 8, 43}},
  {557, {1, -17, 557}},
  {47, {1, 9, 47}},
  {563, {1, 33, 563}},
  {1553, {1, 51, 1553}},
  {53, {1, 9, 53}},
  {1033, {1, -8, 1033}},
  {569, {1, -2, 569}},
  {59, {1, 11, 59}},
  {769, {1, 50, 769}},
  {61, {1, 2, 61}},
  {1087, {1, -55, 1087}},
  {577, {1, 22, 577}},
  {67, {1, 12, 67}},
  {1093, {1, -46, 1093}},
  {991, {1, 13, 991}},
  {71, {1, 8, 71}},
  {73, {1, 14, 73}},
  {587, {1, 26, 587}},
  {1613, {1, -21, 1613}},
  {79, {1, -9, 79}},
  {593, {1, -5, 593}},
  {83, {1, 2, 83}},
  {1109, {1, -49, 1109}},
  {599, {1, -16, 599}},
  {89, {1, -11, 89}},
  {1187, {1, -12, 1187}},
  {1117, {1, -18, 1117}},
  {607, {1, -34, 607}},
  {97, {1, -6, 97}},
  {1123, {1, -45, 1123}},
  {101, {1, -4, 101}},
  {103, {1, 0, 103}},
  {617, {1, 4, 617}},
  {107, {1, -12, 107}},
  {109, {1, -7, 109}},
  {1213, {1, 19, 1213}},
  {113, {1, 5, 113}},
  {787, {1, 1, 787}},
  {1487, {1, 27, 1487}},
  {631, {1, -47, 631}},
  {1453, {1, -36, 1453}},
  {1061, {1, 30, 1061}},
  {127, {1, -7, 127}},
  {641, {1, 18, 641}},
  {131, {1, 7, 131}},
  {1669, {1, -59, 1669}},
  {601, {1, 10, 601}},
  {647, {1, 3, 647}},
  {137, {1, 3, 137}},
  {1607, {1, -7, 1607}},
  {139, {1, -4, 139}},
  {653, {1, -24, 653}},
  {1181, {1, 29, 1181}},
  {823, {1, 31, 823}},
  {659, {1, -14, 659}},
  {149, {1, -6, 149}},
  {151, {1, -9, 151}},
  {157, {1, -4, 157}},
  {709, {1, 8, 709}},
  {673, {1, -16, 673}},
  {163, {1, 4, 163}},
  {1301, {1, 30, 1301}},
  {677, {1, 48, 677}},
  {167, {1, 8, 167}},
  {1193, {1, 32, 1193}},
  {683, {1, 45, 683}},
  {173, {1, 9, 173}},
  {541, {1, -10, 541}},
  {1201, {1, 22, 1201}},
  {1489, {1, 15, 1489}},
  {179, {1, -9, 179}},
  {181, {1, -19, 181}},
  {1097, {1, -58, 1097}},
  {1481, {1, -45, 1481}},
  {1567, {1, 22, 1567}},
  {701, {1, -24, 701}},
  {191, {1, 16, 191}},
  {193, {1, 13, 193}},
  {1483, {1, 24, 1483}},
  {197, {1, 20, 197}},
  {1609, {1, 50, 1609}},
  {199, {1, 13, 199}},
  {859, {1, 41, 859}},
  {1399, {1, -20, 1399}},
  {1229, {1, -18, 1229}},
  {719, {1, 10, 719}},
  {211, {1, 12, 211}},
  {1237, {1, -46, 1237}},
  {727, {1, 33, 727}},
  {1031, {1, 4, 1031}},
  {1231, {1, 14, 1231}},
  {733, {1, 18, 733}},
  {223, {1, -21, 223}},
  {1249, {1, -51, 1249}},
  {1627, {1, -23, 1627}},
  {227, {1, 26, 227}},
  {229, {1, 10, 229}},
  {743, {1, -16, 743}},
  {1559, {1, 7, 1559}},
  {233, {1, -15, 233}},
  {1259, {1, -16, 1259}},
  {239, {1, 15, 239}},
  {241, {1, 29, 241}},
  {1423, {1, 64, 1423}},
  {757, {1, 33, 757}},
  {1321, {1, 26, 1321}},
  {1021, {1, 12, 1021}},
  {761, {1, 38, 761}},
  {251, {1, -8, 251}},
  {1277, {1, 18, 1277}},
  {1279, {1, -62, 1279}},
  {257, {1, 18, 257}},
  {1283, {1, -3, 1283}},
  {773, {1, -3, 773}},
  {263, {1, -2, 263}},
  {1543, {1, 48, 1543}},
  {1289, {1, -19, 1289}},
  {1151, {1, -24, 1151}},
  {269, {1, 14, 269}},
  {271, {1, -18, 271}},
  {1297, {1, 22, 1297}},
  {643, {1, 12, 643}},
  {277, {1, 13, 277}},
  {1373, {1, 45, 1373}},
  {1303, {1, -24, 1303}},
  {281, {1, -24, 281}},
  {283, {1, 16, 283}},
  {1583, {1, 17, 1583}},
  {797, {1, -39, 797}},
  {1523, {1, 23, 1523}},
  {1069, {1, 47, 1069}},
  {293, {1, -9, 293}},
  {1597, {1, 14, 1597}},
  {1319, {1, 27, 1319}},
  {809, {1, 54, 809}},
  {811, {1, -16, 811}},
  {1459, {1, -15, 1459}},
  {1327, {1, 8, 1327}},
  {619, {1, 3, 619}},
  {307, {1, 17, 307}},
  {821, {1, 24, 821}},
  {311, {1, 32, 311}},
  {313, {1, 14, 313}},
  {827, {1, -24, 827}},
  {317, {1, 22, 317}},
  {1579, {1, -31, 1579}},
  {1163, {1, 34, 1163}},
  {829, {1, 40, 829}},
  {839, {1, -15, 839}},
  {1153, {1, -17, 1153}},
  {331, {1, 4, 331}},
  {337, {1, -8, 337}},
  {739, {1, -44, 739}},
  {853, {1, -15, 853}},
  {1367, {1, -65, 1367}},
  {1307, {1, 51, 1307}},
  {857, {1, 22, 857}},
  {347, {1, -12, 347}},
  {349, {1, 35, 349}},
  {863, {1, -21, 863}},
  {523, {1, -6, 523}},
  {353, {1, 34, 353}},
  {571, {1, 31, 571}},
  {1171, {1, -64, 1171}},
  {1381, {1, 48, 1381}},
  {359, {1, 15, 359}},
  {691, {1, -2, 691}},
  {1511, {1, -24, 1511}},
  {1427, {1, 12, 1427}},
  {877, {1, 32, 877}},
  {367, {1, -32, 367}},
  {881, {1, 33, 881}},
  {883, {1, 22, 883}},
  {373, {1, 0, 373}},
  {887, {1, -1, 887}},
  {1291, {1, 42, 1291}},
  {379, {1, -15, 379}},
  {1223, {1, 45, 1223}},
  {383, {1, 27, 383}},
  {1409, {1, -11, 1409}},
  {1637, {1, 58, 1637}},
  {389, {1, 3, 389}},
  {1601, {1, 30, 1601}},
  {1663, {1, 23, 1663}},
  {907, {1, 28, 907}},
  {397, {1, -10, 397}},
  {911, {1, 24, 911}},
  {661, {1, -10, 661}},
  {401, {1, 18, 401}},
  {1091, {1, -17, 1091}},
  {1667, {1, -63, 1667}},
  {1429, {1, -10, 1429}},
  {919, {1, 25, 919}},
  {409, {1, 31, 409}},
  {751, {1, 4, 751}},
  {1439, {1, -45, 1439}},
  {929, {1, 28, 929}},
  {419, {1, 20, 419}},
  {421, {1, -10, 421}},
  {1447, {1, 41, 1447}},
  {937, {1, -37, 937}},
  {1129, {1, -43, 1129}},
  {1451, {1, 24, 1451}},
  {941, {1, -14, 941}},
  {1361, {1, -46, 1361}},
  {431, {1, 30, 431}},
  {433, {1, -12, 433}},
  {947, {1, 17, 947}},
  {439, {1, -10, 439}},
  {953, {1, -48, 953}},
  {443, {1, 32, 443}},
  {1549, {1, 38, 1549}},
  {1471, {1, 56, 1471}},
  {449, {1, -14, 449}},
  {1571, {1, 22, 1571}},
  {967, {1, 23, 967}},
  {457, {1, 21, 457}},
  {971, {1, 2, 971}},
  {461, {1, 6, 461}},
  {463, {1, -16, 463}},
  {977, {1, -10, 977}},
  {467, {1, 11, 467}},
  {1493, {1, 8, 1493}},
  {983, {1, -17, 983}},
  {1103, {1, -36, 1103}},
  {1499, {1, -55, 1499}},
  {479, {1, 26, 479}},
  {1531, {1, 31, 1531}},
  {997, {1, -49, 997}},
  {487, {1, 24, 487}},
  {491, {1, -40, 491}},
  {613, {1, -10, 613}},
  {1009, {1, 23, 1009}},
  {1619, {1, -57, 1619}},
  {499, {1, -8, 499}},
  {1217, {1, 28, 1217}},
  {1013, {1, 46, 1013}},
  {503, {1, 41, 503}},
  {1019, {1, -25, 1019}},
  {509, {1, -42, 509}},
  {1621, {1, 46, 1621}},
};


// compute the Euler poly for p
// just uses the map above
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec __attribute__((unused)), void *param __attribute__((unused)))
{
  acb_poly_zero(poly);
  auto it = euler_factors.find(p);
  if( it != euler_factors.end() ) {
  for(size_t i = 0; i < it->second.size(); ++i)
    acb_poly_set_coeff_si(poly, i, it->second[i]);
  }
}




int main ()
{
  Lfunc_t L;
  double mus[2] = {0, 1};
  Lerror_t ecode;

  // we have a degree 2 L-function of motivic weight 1, so normalisation = 0.5
  L = Lfunc_init(2, 5077, 0.5, mus, &ecode);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  // populate local factors
  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr, ecode);
    return 0;
  }

  // do the computation
  ecode|=Lfunc_compute(L);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  Lfunc *LL=(Lfunc *) L;
  printf("sqrt sign = ");acb_printd(LL->sqrt_sign,20);printf("\n");

  
  // now extract some information
  printf("Order of vanishing = %" PRIu64 "\n",Lfunc_rank(L));
  printf("Sign = ");
  acb_printd(Lfunc_sign(L),DIGITS);
  printf("\n");
  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L),DIGITS);
  printf("\n");

  acb_t ctmp;
  acb_init(ctmp);
  ecode|=Lfunc_special_value(ctmp, L, 1.5, 0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  // ~ 0.510863948654827254260483653705
  ecode|=Lfunc_special_value(ctmp, L, 2.5,0.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr, ecode);
    std::abort();
  }
  printf("L(2.5) = ");acb_printd(ctmp, DIGITS);printf("\n");
  acb_clear(ctmp);


  printf("First 10 zeros\n");
  // we could use Lfunc_zeros(L, 1) for the dual L-function
  arb_srcptr zeros=Lfunc_zeros(L, 0);
  for(int i  = 0; i < 10; ++i) {
    printf("Zero %d = ", i);
    arb_printd(zeros+i, DIGITS);
    printf("\n");
    if (RAW) cout<<"RAW: "<<zeros + i<< endl;
  }

  printf("Z-plot in [0, 10]:\n");
  Lplot_t *Lpp=Lfunc_plot_data(L, 0, 10.0, 20);
  int z = 0;
  double zero_double = arf_get_d(arb_midref(zeros + z), ARF_RND_NEAR);
  for(size_t k=0; k < Lpp->n_points; ++k) {
    printf("%.2f\t%.2f\t", k*Lpp->spacing , Lpp->points[k]);
    int y = 30 + int(7.5*Lpp->points[k]);
    int zero = 30;
    // assuming 60 columns
    for(int i = 0; i < 61; ++i) {
      if(i == y) {
        printf("o");
      } else if (i == zero) {
        printf("|");
      } else if ( (i > zero and i < y) or (i < zero and i > y) ) {
        printf("-");
      } else {
        printf(" ");
      }
    }
    printf("\n");
    if(k*Lpp->spacing < zero_double and (k+1)*Lpp->spacing >= zero_double){
      printf("%.2f\tzero\t", zero_double);
      for(int i = 0; i < 30; ++i)
        printf(" ");
      printf("Z\n");
      zero_double = arf_get_d(arb_midref(zeros + ++z), ARF_RND_NEAR);
    }
  }

  //free memory
  Lfunc_clear_plot(Lpp);
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}
