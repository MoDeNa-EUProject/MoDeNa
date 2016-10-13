#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "spdeclarations.h"

//Definitionen
#define PI   3.1415926535897932384626433832795
#define THRD 0.3333333333333333333333333333333

//externe Variablen
extern int     No_parts;
extern int     Distrfunct;
extern double  Alpha;
extern double  Lambda;
extern double  Mu;
extern double  GSigma;
extern double  Vv;
extern double  Exp;

extern double  Dmin;
extern double  Dmax;

extern int     seed;

extern double  *Percentage;
extern double  *Perc_values;
extern int     Perc_number;
extern double *Diam;

/*
****************************************************************
  
    Gleichverteilung

****************************************************************
*/
double uniform_distr(double small, double large)
{
   double x= 0.0;

   while(x< small)
      x= ((double)rand()/RAND_MAX)*large;

   return x;
}

/*
****************************************************************
  
    Zweipunktverteilung

****************************************************************
*/
double distr_2(double small, double large, double Vv)
{
   double k3= (small/large)*(small/large)*(small/large);
   double P=  (Vv / (Vv + k3));

   return ((double)rand()/RAND_MAX < P) ? small : large;
}

/*
****************************************************************
  
    Gaußverteilung nach Box-Muller-Verfahren

****************************************************************
*/
double gaussianBM(double g_sigma, double mu)
{
   long   seed = 0;
   double erg= 0.0, u1= 0.0, u2= 0.0;
   int i= 0;

   u1= (double)rand()/RAND_MAX;
   u2= (double)rand()/RAND_MAX;

   erg= sqrt(-2*log(u1))*cos(2*PI*u2);

   erg= fabs(g_sigma*erg+mu);

   return erg;
}

/*
****************************************************************
  
    Logarithmische Normalverteilung

****************************************************************
*/
double log_normal(double g_sigma, double mu)
{
   double x= 0.0;

   return exp(gaussianBM(g_sigma, mu));
}

/*
****************************************************************
  
    begrenzte Potenzverteilung

****************************************************************
*/
double power(double small, double large, double a)
{
   static long seed = 0;
   double x= 0.0, x0= pow(small, (a+1)), x1= pow(large, (a+1)), p= 0.0, y= 0.0;

   if (seed == 0)
   {
      time(&seed);
      seed |= 1L;
      srand((seed));
   }

   do
   {
      y=((double)rand()/RAND_MAX);

      x= pow(((x1-x0)*y + x0), 1/(a+1));

      p=((double)rand()/RAND_MAX);
   }
   while(p< 0.95);

   return x;
}

/*
****************************************************************
  
    diskrete Verteilung (Prozentwerte)

****************************************************************
*/
void percentage_dist(void)
{
   double limit= 0;
   int i, j;

   for(j= 0; j< Perc_number; j++)
   {
      *(Percentage+j+1)= *(Percentage+j+1) * (double)No_parts + *(Percentage+j);
      for(i= (int)*(Percentage+j); i< *(Percentage+j+1); i++)
         *(Diam+i)= *(Perc_values+j);
   }
}

/*
****************************************************************
  
    Auswahl der Verteilung

****************************************************************
*/
void setDistribution(int number, double *Diam)
{
   int j= 0;
printf("dist: %d\n", number);
   switch(number)
   {
   case 0://konstanter Durchmesser
      for(j= 0; j< No_parts; j++)
      {
		  *(Diam+j)= 1.0;
		  //printf("%.lf\n", *(Diam+j));
       }
      break;

   case 1://Two-Point
      for(j= 0; j< No_parts; j++)
		  *(Diam+j)= distr_2(Dmin, Dmax, Vv);
      break;

   case 2://uniform distribution - Gleichverteilung
      for(j= 0; j< No_parts; j++)
		   *(Diam+j)= uniform_distr(Dmin, Dmax);
      break;

   case 3://Normalverteilung
      for(j= 0; j< No_parts; j++)
         *(Diam+j)= gaussianBM(GSigma, Mu);
      break;

   case 4://Logarithmische Normalverteilung
      for(j= 0; j< No_parts; j++)
		  *(Diam+j)= log_normal(GSigma, Mu);
      break;

   case 5://Potenzvert.
      for(j= 0; j< No_parts; j++)
         *(Diam+j)= power(Dmin, Dmax, Exp);
      break;

   case 6://Prozentverteilung
      percentage_dist();
      break;

   default://gleiche Kugeln als Standard
      for(j= 0; j< No_parts; j++)
		  *(Diam+j)= 1.0;
      break;
   }
}
