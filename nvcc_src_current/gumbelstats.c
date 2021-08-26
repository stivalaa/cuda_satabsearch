/*****************************************************************************
 * 
 * File:    gumbelstats.c
 * Author:  Alex Stivala
 * Created: July 2010
 *
 * $Id: gumbelstats.c 4709 2013-11-18 00:18:33Z astivala $
 *
 * Functions to compute z-score and p-value from tableau matching score,
 * according to Gumbel distribution parameters.
 * See thesis Ch. 6 (s6.2.2) and Ortiz et al (2002), Abagyan & Batalov (1997),
 * Versetr0m & Taylor (2006), Kolbeck et al (2006), Levitt & Gerstein (1998).
 *
 *****************************************************************************/

#include <math.h>

#include "gumbelstats.h"
/*****************************************************************************
 *
 * Constants
 *
 *****************************************************************************/

/* Euler-Mascheroni constant */
const double eulergamma = 0.5772156649015328606; 

/* pi/sqrt(6) is used in Gumbel distribution functions */
const double pi_over_sqrt6 = M_PI / sqrt(6.0);

/*****************************************************************************
 *
 * Functions
 *
 *****************************************************************************/



/* 
 * z_gumbel() -  compute Z-score from Gumbel distribution
 *
 * Parameters:
 *    x - score to compute Z-score for
 *    a - Gumbel distribution a (location) parameter
 *    b - Gumebel distribution b (scale) parameter
 *
 * Return value:
 *    Z-score computed for x according to Gumbel(a,b) distribution
 */
double z_gumbel(int x, double a, double b)
{

  double mu, sigma, z;
  mu = a + b * eulergamma;
  sigma = ( pi_over_sqrt6 * b );
  z =  (x - mu)/sigma;
  return z;
}

/*
 * pv_gumbel() - compute P-value for Z-score from Gumbel distribution
 *
 * Parameters: 
 *    z - z-score from z_gumbel()
 *
 * Return value:
 *    P-value for the Z-score
 */
double pv_gumbel(double z)
{
  return ( 1 - exp(-exp(-( (pi_over_sqrt6 * z + eulergamma) ) ) ) );
}

/*
 * norm2() - normalize tableau matching score by size
 *
 *  Normalization similar to norm2 in Pelta et al 2008 (from Xie & Sahinidis
 *   2006) for MAX-CMO.
 *
 *  norm2(struct1,struct2) = 2*tabmatch_score(struct1,struct2) / 
 *                               (#sses(struct1) + #sses(struct2))
 *
 *   Parameters:
 *       score  - tableau match score for the two structures
 *       size1  - number of SSEs in one structure
 *       size2  - number of SSEs in other structure
 *
 *   Return value:
 *       normalized score as described above.
 */
double norm2(int score, int size1, int size2)
{
  return 2.0 * score / ((double)(size1 + size2));
}
