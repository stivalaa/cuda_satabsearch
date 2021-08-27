#ifndef GUMBELSTATS_H
#define GUMBELSTATS_H
/*****************************************************************************
 * 
 * File:    gumbelstats.h
 * Author:  Alex Stivala
 * Created: July 2010
 *
 * $Id: gumbelstats.h 4709 2013-11-18 00:18:33Z astivala $
 *
 * Functions to compute z-score and p-value from tableau matching score,
 * according to Gumbel distribution parameters.
 * See thesis Ch. 6 (s6.2.2) and Ortiz et al (2002), Abagyan & Batalov (1997),
 * Versetr0m & Taylor (2006), Kolbeck et al (2006), Levitt & Gerstein (1998).
 *
 *****************************************************************************/


/* Gumbel distribution parameters estimated by MLE for structures in 
   different folds (see thesis Ch. 6 (s.6.2.2). */
/* These for 4096 restart on query200 */
const double gumbel_a = 0.3780327676087335;
const double gumbel_b = 0.3582596175507505;

/* compute z_score from  from tableau matching score */
double z_gumbel(int x, double a, double b);

/* compute p-value from z-score */
double pv_gumbel(double z);

/* normalize tableau match scores by size */
double norm2(int score, int size1, int size2);

#endif /* GUMBELSTATS_H */
