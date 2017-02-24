 /*"
 License
     This file is part of Modena.

     Modena is free software; you can redistribute it and/or modify it under
     the terms of the GNU General Public License as published by the Free
     Software Foundation, either version 3 of the License, or (at your option)
     any later version.

     Modena is distributed in the hope that it will be useful, but WITHOUT ANY
     WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
     FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.

     You should have received a copy of the GNU General Public License along
     with Modena.  If not, see <http://www.gnu.org/licenses/>.

 Description
     Tool for Packing Spheres
****************************************************************

    Kugelpackprogramm

    packt Kugeln gleicher und ungleicher Durchmesser durch
    eine Rearrangementmethode nach Jodrey&Tory
    und Bargiel&Moscinski

    programmiert: Antje Elsner

    (Listenstruktur aus:
      Moscinski, J., and Bargiel, M., "C-language program
      for simulation of irregular close packing of hard spheres",
      Comp Phys Comm, 64 (1991) 183-192
      http://cpc.cs.qub.ac.uk/summaries/ABZF_v1_0.html)

    Version: 12.28

    Datum: 28.08.2008

****************************************************************
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <limits.h>

#include "spdeclarations.h"

//Definitionen
#define FANZ   7
#define TRUE   1
#define FALSE  0

#define PI        3.141592653589793238462643383279

//Flags
int    Flag_Random=  0;
int    Flag_End   =  FALSE;

short  Time_out[6];

//Dateierweiterungen
const char gen_par_psfx[]  = ".prj"; //.prj - generelle Parameterdatei
const char beg_coord_psfx[]= ".sco"; //.sco - Anfangszustand bei zuf�lliger Erzeugung
const char res_coord_psfx[]= ".rco"; //.rco - Ergebniskoordinaten
const char inp_coord_psfx[]= ".pco"; //.pco - vorgegebene Koordinatendatei
const char res_data_psfx[] = ".rst"; //.rst - Ergebnisdaten

//Dateinamen
char prj_name[80];            //Projektname
char gen_par_name[80];        //.prj - generelle Parameterdatei
char beg_coord_name[80];      //.sco - Anfangszustand bei zuf�lliger Erzeugung
char res_coord_name[80];      //.rco - Ergebniskoordinaten
char inp_coord_name[80];      //.pco - vorgegebene Koordinatendatei
char res_data_name[80];       //.rst - Ergebnisdaten

//Dateipointer
FILE* gen_par_file= NULL;        //.prj - generelle Parameterdatei
FILE* beg_coord_file= NULL;      //.sco - Anfangszustand bei zuf�lliger Erzeugung
FILE* res_coord_file= NULL;      //.rco - Ergebniskoordinaten
FILE* inp_coord_file= NULL;      //.pco - vorgegebene Koordinatendatei
FILE* res_data_file= NULL;       //.rst - Ergebnisdaten

//Globale Variablen
int     No_parts;
int     No_cells_x, No_cells_y, No_cells_z;
int     No_cells, Ncell_min;
long    Max_steps= 100000, Ntau= 1638400;
double  Epsilon= 0.6, Epsilon_scl;
double  Din, Dout0, Dout, Dout2;
double  Pactual, Pnom0= 0.6, Pnomin, Pactual_old;
double  Diam_dens;
double  Box_x, Box_y, Box_z, Half_x, Half_y, Half_z;
double  Relax;
double  Diam_min, Diam_max;
double  Dmin, Dmax;

//Speicherbereiche
int     *Link_head, *Link_list;
int     *Ncell_bound_x, *Ncell_bound_y, *Ncell_bound_z;
double  *Pbc_x, *Pbc_y, *Pbc_z;
double  *X, *Y, *Z;
double  *Diam;
double  *Force_x, *Force_y, *Force_z;
double  *k_freq, *l_freq, *g_r;

//Parameter f�r Verteilungsfunktionen
int     Distrfunct= 0;
double  Alpha  = 2.0;
double  Lambda = 0.04;
double  Mu     = 1.0;
double  GSigma = 1.0;
double  Vv     = 0.5;
double  Exp    = -0.3;
double  *Percentage;
double  *Perc_values;
int     Perc_number= 0;

int     seed= 0;

long    Step;

/*
****************************************************************

    Hauptfunktion

****************************************************************
*/
int main(int argc, char * argv[])
{
   strcpy(prj_name, "Project01");

   if (argc == 2)
   {
		strcpy(prj_name, argv[1]);
   }

   make_project();

   read_param();

   set_coordinates();

   list_structure();

   packing_loop();

   write_sphere_data();

   end_run();

   return 0;
}

/*
****************************************************************

    Anlegen und �ffnen der Projektdateien

****************************************************************
*/
void make_project(void)
{
   strcpy(gen_par_name, prj_name); strcat(gen_par_name, gen_par_psfx);

   strcpy(beg_coord_name, prj_name); strcat(beg_coord_name, beg_coord_psfx);
   strcpy(res_coord_name, prj_name); strcat(res_coord_name, res_coord_psfx);
   strcpy(inp_coord_name, prj_name); strcat(inp_coord_name, inp_coord_psfx);
   strcpy(res_data_name, prj_name); strcat(res_data_name, res_data_psfx);

   if (( gen_par_file= fopen(gen_par_name, "r")) == NULL)
   {
      printf("Keine Steuerparameterdatei %s.prj gefunden.\n", prj_name);
      fehlerAusgabe("read_param", -3);
   }

   if (( beg_coord_file= fopen(beg_coord_name, "w") ) == NULL)
   {
      printf("Konnte %s.sco nicht erzeugen.\n", prj_name);
      fehlerAusgabe("read_param", -3);
   }

   if (( inp_coord_file= fopen(inp_coord_name, "r") ) == NULL)
   {
      printf("Keine Datendatei %s.pco gefunden.\n", prj_name);
  //    fehlerAusgabe("read_param", -3);
   }

   if (( res_data_file= fopen(res_data_name, "w") ) == NULL)
   {
      printf("Konnte %s.rst nicht erzeugen.\n", prj_name);
      fehlerAusgabe("read_param", -3);
   }
}

/*
****************************************************************

    Einlesen der Parameter

****************************************************************
*/
int read_param(void)
{
   int k= 0, m= 0, l= 0;
   double sum= 0;
   char line[80];

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //zweite Zeile
   fgets(line, 80, gen_par_file);
   l++;

   if ((Flag_Random = atoi(line))== EOF)
   {
      printf("Fehler in Zeile %d beim Einlesen von Flag_Random.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //vierte Zeile: Anzahl Kugeln
   fgets(line, 80, gen_par_file);
   l++;

   if (!(No_parts = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //6. 7. 8. Zeile Boxl�ngen
   fgets(line, 80, gen_par_file);
   l++;

   if (!(No_cells_x = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(No_cells_y = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(No_cells_z = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //Epsilon
   fgets(line, 80, gen_par_file);
   l++;

   if (!(Epsilon = atof(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Epsilon <= 0)
      {
         printf("Parameter Epsilon muss groesser Null sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //Ntau
   fgets(line, 80, gen_par_file);
   l++;

   if (!(Ntau = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Ntau <= 0)
      {
         printf("Parameter Ntau muss groesser Null sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //nom. Dichte
   fgets(line, 80, gen_par_file);
   l++;

   if (!(Pnom0 = atof(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Pnom0 <= 0 || Pnom0 > 1.0)
      {
         printf("Parameter Pnom muss zwischen Null  und Eins liegen.\n");
         fehlerAusgabe("read_param", -1);
      }

   }

   //Zeile ist ein Kommentar
   fgets(line, 80, gen_par_file);
   l++;

   //min. und max. Durchmesser
   fgets(line, 80, gen_par_file);
   l++;

   if (!(Diam_max = atof(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Diam_max <= 0)
      {
         printf("Maximaler Durchmesser muss groesser Null sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(Diam_min = atof(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Diam_min <= 0 || Diam_min > Diam_max)
      {
         printf("Minimaler Durchmesser muss groesser Null und kleiner als maximaler Durchmesser sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   //Max. Anzahl Schritte
   fgets(line, 80, gen_par_file);
   l++;

   fgets(line, 80, gen_par_file);
   l++;

   if (!(Max_steps = atoi(line)))
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Max_steps < 0)
      {
         printf("Maximale Anzahl Verdichtungsschritte muss positiv sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   //Verteilungsfunktion
   fgets(line, 80, gen_par_file);
   l++;

   fgets(line, 80, gen_par_file);
   l++;

   if ((Distrfunct = atoi(line))== EOF)
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Distrfunct < 0)
      {
         printf("Nummer der Verteilungsfunktion muss positiv oder 0 sein.\nEs werden gleich gro�e Kugeln erzeugt.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   if (Flag_Random == 1)
      Distrfunct= 11;

   //weitere Parameter f�r die Verteilungen
   fgets(line, 80, gen_par_file);
   l++;

   fgets(line, 80, gen_par_file);
   l++;

   if (!(Mu = atof(line)))//
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Mu < 0)
      {
         printf("Mittelwert fuer Normal- oder Log-Normalvert. muss positiv sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(GSigma = atof(line)))//
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(GSigma < 0)
      {
         printf("Standardabweichung muss positiv sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(Vv = atof(line)))//
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Vv < 0 || Vv > 1)
      {
         printf("Volumenverh�ltnis muss positiv sein und zwischen 0 und 1 liegen.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   fgets(line, 80, gen_par_file);
   l++;

   if (!(Exp = atof(line)))//
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }

   //Prozentbins
   fgets(line, 80, gen_par_file);
   l++;

   if (!(Perc_number = atoi(line)))//
   {
      printf("Fehler in Zeile %d.\n", l);
      fehlerAusgabe("read_param", -1);
   }
   else
   {
      if(Perc_number < 1)
      {
         printf("Anzahl Prozentbins muss groesser 1 sein.\n");
         fehlerAusgabe("read_param", -1);
      }
   }

   Percentage   = (double *) calloc(Perc_number+2, sizeof(double));
   Perc_values  = (double *) calloc(Perc_number+1, sizeof(double));

   if (Distrfunct == 6)
   {
      *Percentage= 0;

      for (k= 1; k<= Perc_number; k++)
      {
         fgets(line, 80, gen_par_file);
         l++;

         if (!(*(Percentage+k) = atof(line)))//
         {
            printf("Fehler in Zeile %d.\n", l);
            fehlerAusgabe("read_param", -1);
         }

         sum+= *(Percentage+k);
      }

      if (sum != 1.0)
      {
         printf("\nProgrammausfuehrung fehlerhaft beendet.\nBitte geben Sie eine gueltige Prozentverteilung an!\n");
         exit(1);
      }

      for (m= 0; m < Perc_number; m++)
      {
         fgets(line, 80, gen_par_file);
         l++;

         if (!(*(Perc_values+m) = atof(line)))//
         {
            printf("Fehler in Zeile %d.\n", l);
            fehlerAusgabe("read_param", -1);
         }
      }
   }

   fclose(gen_par_file);
   gen_par_file = NULL;

   return 0;
}

/*
****************************************************************

    Belegen der x,y,z-Koordinaten und Durchmesser
    durch Einlesen oder zuf�lliges Erzeugen

****************************************************************
*/
int set_coordinates(void)
{
   int  ipart= 0, cnt= 0, r;
   double x, y, z, d, max_x= 0, max_y= 0, max_z= 0;
   char line[128];

   static long seed = 0;
   time(&seed);
   srand(seed);

   if (Flag_Random != 0)
   {
      while(!feof(inp_coord_file))
      {
         if(fscanf(inp_coord_file,"%lf %lf %lf %lf", &x, &y, &z, &d) > 1)
            cnt++;
      }

      No_parts = cnt;

      rewind(inp_coord_file);
   }

   Diam    = (double *) calloc(No_parts, sizeof(double));
   X       = (double *) calloc(No_parts, sizeof(double));
   Y       = (double *) calloc(No_parts, sizeof(double));
   Z       = (double *) calloc(No_parts, sizeof(double));

   Force_x = (double *) calloc(No_parts, sizeof(double));
   Force_y = (double *) calloc(No_parts, sizeof(double));
   Force_z = (double *) calloc(No_parts, sizeof(double));

   if (Flag_Random != 0)
   {
      r= fscanf(inp_coord_file,"%lf %lf %lf %lf", &x, &y, &z, &d);

      if(r != 4)
      {
         printf("Fehler in Koordinatendatei in Zeile 0.\n");
         fehlerAusgabe("set_coordinates", -1);
      }

      ipart= 0;

      while(ipart < No_parts)
      {
         if(max_x< x) max_x= x;
         if(max_y< y) max_y= y;
         if(max_z< z) max_z= z;

         *(Diam+ipart)= d;
         *(X+ipart) = x;
         *(Y+ipart) = y;
         *(Z+ipart) = z;

         if (*(X+ipart) < 0.0  ||  *(X+ipart) > No_cells_x  ||
             *(Y+ipart) < 0.0  ||  *(Y+ipart) > No_cells_y  ||
             *(Z+ipart) < 0.0  ||  *(Z+ipart) > No_cells_z   )
         {
            fehlerAusgabe("set_coordinates", -1);
         }

         r= fscanf(inp_coord_file,"%lf %lf %lf %lf", &x, &y, &z, &d);

         if(r == 4)
            ipart++;
         else
            break;
      }

      if (max_x <= 1.0 && max_y <= 1.0 && max_z <= 1.0)
      {
         for (ipart = 0;  ipart < No_parts;  ipart++)
         {
            *(X+ipart) *= No_cells_x;
            *(Y+ipart) *= No_cells_y;
            *(Z+ipart) *= No_cells_z;

            if (*(X+ipart) < 0.0  ||  *(X+ipart) > No_cells_x  ||
                *(Y+ipart) < 0.0  ||  *(Y+ipart) > No_cells_y  ||
                *(Z+ipart) < 0.0  ||  *(Z+ipart) > No_cells_z   )
            {
               fehlerAusgabe("set_coordinates", -1);
            }
         }
      }
      sort();

      fclose(inp_coord_file);
      inp_coord_file = NULL;
   }
   else
   {
      for (ipart = 0;  ipart < No_parts; ipart++)
      {
         *(X+ipart) = ((double)rand()/(double)RAND_MAX) * (double)No_cells_x;
         *(Y+ipart) = ((double)rand()/(double)RAND_MAX) * (double)No_cells_y;
         *(Z+ipart) = ((double)rand()/(double)RAND_MAX) * (double)No_cells_z;

         if (*(X+ipart) < 0.0  ||  *(X+ipart) > No_cells_x  ||
             *(Y+ipart) < 0.0  ||  *(Y+ipart) > No_cells_y  ||
             *(Z+ipart) < 0.0  ||  *(Z+ipart) > No_cells_z   )
            {
               fehlerAusgabe("set_coordinates", -1);
            }
      }
   }

   return 0;
}

/*
****************************************************************

    Erzeugen der Durchmesserliste und Anlegen der Hilfslisten
    f�r die Verwaltung der Zellzerlegung und
    der periodischen Randbedingungen

    -> Verwaltungslisten progr. nach:
      Moscinski, J., and Bargiel, M., "C-language program
      for simulation of irregular close packing of hard spheres",
      Comp Phys Comm, 64 (1991) 183-192
      http://cpc.cs.qub.ac.uk/summaries/ABZF_v1_0.html

****************************************************************
*/
int list_structure(void)
{
   int    ipart, i, imult;
   int    ndim_x, ndim_y, ndim_z;
   double corr, Summe_Kubikdurchm;
   double Power      = 0.3333333333333333333333333333333;
   double Sphere_vol = 1.9098593171027440292266051604702;

   No_cells   = No_cells_x * No_cells_y * No_cells_z;
   Ncell_min  = (No_cells_x < No_cells_y) ? No_cells_x : No_cells_y;
   Ncell_min  = (No_cells_z < Ncell_min) ? No_cells_z : Ncell_min;
   Box_x      = (double) No_cells_x;
   Box_y      = (double) No_cells_y;
   Box_z      = (double) No_cells_z;
   Half_x     = 0.5 * Box_x;
   Half_y     = 0.5 * Box_y;
   Half_z     = 0.5 * Box_z;
   Diam_dens  = No_cells * Sphere_vol;

   Diam_max /= Diam_min;
   Diam_min = 1.0;

   Dmin= Diam_min;
   Dmax= Diam_max;

   if (Flag_Random == 0)
   {
      setDistribution(Distrfunct, Diam);
      qsort ((void *) Diam, No_parts, sizeof(double), dblcmp);
   }

   for (ipart = 0; ipart < No_parts; ipart++)
      fprintf(beg_coord_file, "%.6f\t%.6f\t%.6f\t%.6f\n", *(X+ipart), *(Y+ipart), *(Z+ipart), *(Diam+ipart));

   fflush(beg_coord_file);

   fclose(beg_coord_file);
   beg_coord_file = NULL;

   Summe_Kubikdurchm =  0.0;

   for (ipart = 0; ipart < No_parts; ipart++)
      Summe_Kubikdurchm += *(Diam+ipart) * *(Diam+ipart) * *(Diam+ipart);

   Diam_dens /= Summe_Kubikdurchm;

   Pnomin       = Pnom0;
   Dout0        = pow (Diam_dens * Pnom0 , Power);
   Dout         = Dout0;
   Relax        = 0.5 * Dout / Ntau;
   Dout2        = Dout * Dout;
   Epsilon_scl  = Epsilon / Dout0;

   Link_head =  (int *) calloc(No_cells, sizeof(int));
   Link_list =  (int *) calloc(No_parts, sizeof(int));

   ndim_x  =  3 * No_cells_x - 2;
   ndim_y  =  3 * No_cells_y - 2;
   ndim_z  =  3 * No_cells_z - 2;

   Ncell_bound_x  =  (int *) calloc(ndim_x+1, sizeof(int));
   if (Ncell_bound_x == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   Ncell_bound_y  =  (int *) calloc(ndim_y+1, sizeof(int));
   if (Ncell_bound_y == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   Ncell_bound_z  =  (int *) calloc(ndim_z+1, sizeof(int));
   if (Ncell_bound_z == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   Pbc_x =  (double *) calloc(ndim_x+1, sizeof(double));
   if (Pbc_x == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   Pbc_y =  (double *) calloc(ndim_y+1, sizeof(double));
   if (Pbc_y == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   Pbc_z =  (double *) calloc(ndim_z+1, sizeof(double));
   if (Pbc_z == NULL)
   {
      fehlerAusgabe("list_structure", -2);
   }

   imult = 1;
   corr  = -No_cells_x;

   for (i= 1; i <= ndim_x; i++)
   {
      *(Ncell_bound_x+i) = imult * No_cells_y * No_cells_z;
      *(Pbc_x+i)         = corr;

      ++imult;

      if (imult >= No_cells_x)
      {
         imult = 0;
         corr += No_cells_x;
      }
   }

   imult = 1;
   corr  = -No_cells_y;

   for (i = 1; i <= ndim_y; i++)
   {
      *(Ncell_bound_y+i) = imult * No_cells_z;
      *(Pbc_y+i)         = corr;

      ++imult;

      if (imult >= No_cells_y)
      {
         imult = 0;
         corr += No_cells_y;
      }
   }

   imult = 1;
   corr  = -No_cells_z;

   for (i= 1; i <= ndim_z; i++)
   {
      *(Ncell_bound_z+i) = imult;
      *(Pbc_z+i)         = corr;

      ++imult;

      if (imult >= No_cells_z)
      {
         imult = 0;
         corr += No_cells_z;
      }
   }

   return 0;
}

/*
****************************************************************

    Hauptschleife

****************************************************************
*/
int packing_loop(void)
{
   date_time (res_data_file, (long *) Time_out);

   if (Max_steps > 0)
   {
      interaction();

      Pactual = Din*Din*Din / Diam_dens;

      for (Step= 1; Step <= Max_steps && Flag_End == FALSE; Step++)
      {
         Dout       -=  Relax;

         Epsilon_scl = Epsilon / Dout;

         Dout2       =  Dout * Dout;

         motion();

         interaction();

         stepon();
      }
      result(Step-1);
   }

   return 0;
}

/*
****************************************************************

    Errechnen der aktuellen Dichten, Erneuern der Relaxationsrate

****************************************************************
*/
int stepon(void)
{
   int  iprec, Nsf= 10;

   Pnomin    =  Dout * Dout * Dout / Diam_dens;
   Pactual   =  Din  * Din  * Din  / Diam_dens;

   if (Dout <= Din)
   {
      double dsave = Dout;

      Dout  = 1.1 * Din;

      Dout2 = Dout * Dout;

      interaction();

      Dout  = dsave;

      Dout2 = Dout * Dout;

      if (Din * *(Diam+No_parts-1) > Ncell_min)
         Din = Ncell_min / *(Diam+No_parts - 1);

      Pactual   =  Din  * Din  * Din  / Diam_dens;

      Flag_End= TRUE;

      return 0;
   }

   iprec= (int) -log10 (Pnomin - Pactual);

   if (iprec >= Nsf)
   {
      Relax  *=  0.5;

      Nsf++;
   }
   return 0;
}

/*
****************************************************************

    Interaktionen zwischen den Kugeln

****************************************************************
*/
int interaction(void)
{
   int    ipart, icell;
   double Rcut;

   Din  =  Dout * Dout;

   for(icell= 0; icell < No_cells; icell++)
   {
      *(Link_head+icell)= -1;
   }

   for(ipart= 0; ipart < No_parts; ipart++)
   {
      *(Force_x+ipart)= 0.0;
      *(Force_y+ipart)= 0.0;
      *(Force_z+ipart)= 0.0;
   }

   for (ipart= 0; ipart < No_parts; ipart++)
   {
      Rcut= ((ipart == 0) ? 0.0 : 0.5) * Dout * (*(Diam + ipart - 1) + *(Diam + ipart));

      if ((int) (2.0 * Rcut) < Ncell_min - 1)
         force_part(ipart, Rcut);
      else
         force_all(ipart);
   }

   Din  =  sqrt(Din);

   return 0;
}

/*
****************************************************************

    Neuanordnen der Zentren

****************************************************************
*/
void motion(void)
{
   int    i;
   double x_buff, y_buff, z_buff;

   for (i = 0; i < No_parts; i++)
   {
      x_buff = *(X+i)+(*(Force_x+i) * Epsilon_scl) / *(Diam+i);
      y_buff = *(Y+i)+(*(Force_y+i) * Epsilon_scl) / *(Diam+i);
      z_buff = *(Z+i)+(*(Force_z+i) * Epsilon_scl) / *(Diam+i);

      if (x_buff < 0.0)
         x_buff+= (ceil(-(x_buff/Box_x))*Box_x);

      if (x_buff >= Box_x)
         x_buff-= (floor(x_buff/Box_x)*Box_x);

      if (y_buff < 0.0)
         y_buff+= (ceil(-(y_buff/Box_y))*Box_y);

      if (y_buff >= Box_y)
         y_buff-= (floor(y_buff/Box_y)*Box_y);

      if (z_buff < 0.0)
         z_buff+= (ceil(-(z_buff/Box_z))*Box_z);

      if (z_buff >= Box_z)
         z_buff-= (floor(z_buff/Box_z)*Box_z);

      *(X+i) = x_buff;
      *(Y+i) = y_buff;
      *(Z+i) = z_buff;
   }
}

/*
****************************************************************

    Verschiebungen bei Nutzung der Nachbarschaftslisten

****************************************************************
*/
void force_part(int ipart_p, double Rcut)
{
   int    jpart= 0, ifirst= 0;
   int    leap_x= 0, leap_y= 0, leap_z= 0;
   int    icell_x= 0, icell_xy= 0, icell_y= 0, icell_z= 0, icell= 0;
   int    low_cell_x, low_cell_y, low_cell_z;
   int    lim_cell_x, lim_cell_y, lim_cell_z;
   int    idif_x= 0, idif_y= 0, idif_z= 0, idif= 0;

   double pos_x, pos_y, pos_z, dist;
   double diam_i, dsum, sigma;
   double dif_x, dif_y, dif_z, difpot;
   double f_x, f_y, f_z;

   pos_x  =  *(X+ipart_p);
   pos_y  =  *(Y+ipart_p);
   pos_z  =  *(Z+ipart_p);
   diam_i =  *(Diam+ipart_p);

   low_cell_x = (int) (pos_x + No_cells_x - Rcut);
   low_cell_y = (int) (pos_y + No_cells_y - Rcut);
   low_cell_z = (int) (pos_z + No_cells_z - Rcut);
   lim_cell_x = (int) (pos_x + No_cells_x + Rcut);
   lim_cell_y = (int) (pos_y + No_cells_y + Rcut);
   lim_cell_z = (int) (pos_z + No_cells_z + Rcut);
   idif_x  =  (*(Ncell_bound_x+lim_cell_x) > *(Ncell_bound_x+low_cell_x))  ? 0 : 4;
   idif_y  =  (*(Ncell_bound_y+lim_cell_y) > *(Ncell_bound_y+low_cell_y))  ? 0 : 2;
   idif_z  =  (*(Ncell_bound_z+lim_cell_z) > *(Ncell_bound_z+low_cell_z))  ? 0 : 1;
   idif    =  idif_x + idif_y + idif_z;

   switch (idif)
   {
      case 0:
      for (leap_x= low_cell_x; leap_x<= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y<= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z<= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x;
                  dif_y  =  *(Y+jpart) - pos_y;
                  dif_z  =  *(Z+jpart) - pos_z;

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 1:
      for (leap_x=low_cell_x; leap_x<=lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y=low_cell_y; leap_y<=lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z=low_cell_z; leap_z<=lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x;
                  dif_y  =  *(Y+jpart) - pos_y;
                  dif_z  =  *(Z+jpart) - pos_z + *(Pbc_z+leap_z);

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 2:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x;
                  dif_y  =  *(Y+jpart) - pos_y + *(Pbc_y+leap_y);
                  dif_z  =  *(Z+jpart) - pos_z;

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 3:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x;
                  dif_y  =  *(Y+jpart) - pos_y + *(Pbc_y+leap_y);
                  dif_z  =  *(Z+jpart) - pos_z + *(Pbc_z+leap_z);

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
              }
           }
        }
      }
      break;

      case 4:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x + *(Pbc_x+leap_x);
                  dif_y  =  *(Y+jpart) - pos_y;
                  dif_z  =  *(Z+jpart) - pos_z;

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 5:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x + *(Pbc_x+leap_x);
                  dif_y  =  *(Y+jpart) - pos_y;
                  dif_z  =  *(Z+jpart) - pos_z + *(Pbc_z+leap_z);

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart  = *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 6:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x + *(Pbc_x+leap_x);
                  dif_y  =  *(Y+jpart) - pos_y + *(Pbc_y+leap_y);
                  dif_z  =  *(Z+jpart) - pos_z;

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart= *(Link_list+jpart);
               }
            }
         }
      }
      break;

      case 7:
      for (leap_x= low_cell_x; leap_x <= lim_cell_x; leap_x++)
      {
         icell_x = *(Ncell_bound_x+leap_x);

         for (leap_y= low_cell_y; leap_y <= lim_cell_y; leap_y++)
         {
            icell_xy = *(Ncell_bound_y+leap_y) + icell_x;

            for (leap_z= low_cell_z; leap_z <= lim_cell_z; leap_z++)
            {
               icell = *(Ncell_bound_z+leap_z) + icell_xy;
               jpart = *(Link_head+icell);

               while (jpart != -1)
               {
                  dif_x  =  *(X+jpart) - pos_x + *(Pbc_x+leap_x);
                  dif_y  =  *(Y+jpart) - pos_y + *(Pbc_y+leap_y);
                  dif_z  =  *(Z+jpart) - pos_z + *(Pbc_z+leap_z);

                  dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
                  dsum   =  0.5 * (diam_i + *(Diam+jpart));
                  sigma  =  Dout * dsum;

                  if (dist < sigma)
                  {
                     if (dist < sqrt(Din)*dsum)
                        Din = (dist*dist)/(dsum*dsum);

                     difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

                     f_x    = difpot * dif_x;
                     f_y    = difpot * dif_y;
                     f_z    = difpot * dif_z;

                     *(Force_x+ipart_p) -= f_x;
                     *(Force_y+ipart_p) -= f_y;
                     *(Force_z+ipart_p) -= f_z;
                     *(Force_x+jpart)   += f_x;
                     *(Force_y+jpart)   += f_y;
                     *(Force_z+jpart)   += f_z;
                  }
                  jpart= *(Link_list+jpart);
               }
            }
         }
      }
      break;
   }

   icell_x  =  (int) pos_x;
   icell_y  =  (int) pos_y;
   icell_z  =  (int) pos_z;
   icell    =  No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
   ifirst   =  *(Link_head+icell);
   *(Link_list+ipart_p)  =  ifirst;
   *(Link_head+icell)    =  ipart_p;
}

/*
****************************************************************

    Berechnen der Verschiebung f�r alle Zellen

****************************************************************
*/
void force_all(int ipart_p)
{
   int    jpart= 0, ifirst= 0;
   int    icell_x= 0, icell_y= 0, icell_z= 0, icell= 0;

   double pos_x= 0.0, pos_y= 0.0, pos_z= 0.0, dist= 0.0, sigma= 0.0;
   double diam_i= 0.0, dsum= 0.0;
   double dif_x= 0.0, dif_y= 0.0, dif_z= 0.0, difpot= 0.0;
   double f_x= 0, f_y= 0, f_z= 0;

   pos_x  =  *(X+ipart_p);
   pos_y  =  *(Y+ipart_p);
   pos_z  =  *(Z+ipart_p);
   diam_i =  *(Diam+ipart_p);

   for (icell= 0; icell < No_cells; icell++)
   {
      jpart = *(Link_head+icell);

      while (jpart != -1)
      {
         dif_x  =  *(X+jpart) - pos_x;
         dif_y  =  *(Y+jpart) - pos_y;
         dif_z  =  *(Z+jpart) - pos_z;

         if (dif_x < -Half_x)
            dif_x += Box_x;

         if (dif_y < -Half_y)
            dif_y += Box_y;

         if (dif_z < -Half_z)
            dif_z += Box_z;

         if (dif_x >  Half_x)
            dif_x -= Box_x;

         if (dif_y >  Half_y)
            dif_y -= Box_y;

         if (dif_z >  Half_z)
            dif_z -= Box_z;

         dist   =  sqrt(dif_x * dif_x + dif_y * dif_y + dif_z * dif_z);
         dsum   =  0.5 * (diam_i + *(Diam+jpart));
         sigma  =  Dout * dsum;

         if (dist < sigma)
         {
            if (dist < sqrt(Din)*dsum)
               Din = (dist*dist)/(dsum*dsum);

            difpot= (0.5 * Dout2 * diam_i * (*(Diam+jpart))) * (1.0 - (dist * dist) / (sigma * sigma))/dist;

            f_x    = difpot * dif_x;
            f_y    = difpot * dif_y;
            f_z    = difpot * dif_z;

            *(Force_x+ipart_p) -= f_x;
            *(Force_y+ipart_p) -= f_y;
            *(Force_z+ipart_p) -= f_z;
            *(Force_x+jpart)   += f_x;
            *(Force_y+jpart)   += f_y;
            *(Force_z+jpart)   += f_z;
         }
         jpart  = *(Link_list+jpart);
      }
   }

   icell_x  =  (int) pos_x;
   icell_y  =  (int) pos_y;
   icell_z  =  (int) pos_z;
   icell    =  No_cells_z * (No_cells_y * icell_x + icell_y) + icell_z;
   ifirst   =  *(Link_head+icell);
   *(Link_list+ipart_p)  =  ifirst;
   *(Link_head+icell)    =  ipart_p;
}


/*
****************************************************************

    Hilfsfunktion Tauschen

****************************************************************
*/
void swap(int i, int j)
{
   double temp;
   temp = *(Diam+j),   *(Diam+j)  = *(Diam+i),  *(Diam+i)  = temp;
   temp = *(X+j),      *(X+j)     = *(X+i),     *(X+i)     = temp;
   temp = *(Y+j),      *(Y+j)     = *(Y+i),     *(Y+i)     = temp;
   temp = *(Z+j),      *(Z+j)     = *(Z+i),     *(Z+i)     = temp;
}

/*
****************************************************************

    Hilfsfunktion Sortieren

****************************************************************
*/
void sort(void)
{
	int gap, i, j;
	for(gap = No_parts/2; gap > 0; gap /= 2)
		for(i = gap; i < No_parts; i++)
			for(j = i - gap; j >= 0 && *(Diam+j) > *(Diam+j+gap); j -= gap)
				swap(j, j+gap);
}

/*
****************************************************************

    Ausgabe des Packungsergebnisses

****************************************************************
*/
void result(int nstep)
{
   int noDist;
   char* DistFunct[]= {"constant","two-point","uniform","gaussian","log-normal","power","percentage","predefined (file)"};

   if(Distrfunct == 11)
      noDist= 7;
   else
      noDist= Distrfunct;

   fprintf(res_data_file, "----------------------------------------\n");
   if(Distrfunct == 11)
      fprintf(res_data_file,"\nDistribution     : %s", "predefined (file)");
   else
      fprintf(res_data_file,"\nDistribution     : %s", DistFunct[Distrfunct]);

   fprintf(res_data_file,"\nCells in x-dir.  : %d", No_cells_x);
   fprintf(res_data_file,"\nCells in y-dir.  : %d", No_cells_y);
   fprintf(res_data_file,"\nCells in z-dir.  : %d", No_cells_z);
   fprintf(res_data_file,"\nNumber of Spheres: %d", No_parts);
   fprintf(res_data_file,"\nDensity          : %5.13f",Pactual);
   fprintf(res_data_file,"\nCalculation steps: %d\n", nstep);
   fprintf(res_data_file, "----------------------------------------\n");

   date_time (res_data_file, (long *) Time_out);

   fflush(res_data_file);

   return;
}

/*
****************************************************************

    Freigeben von Speicherbereichen und Schlie�en von Dateien

****************************************************************
*/
void end_run(void)
{
   if (Link_head != NULL)
      free (Link_head);

   if (Link_list != NULL)
      free (Link_list);

   if (Force_x != NULL)
      free (Force_x);
   if (Force_y != NULL)
      free (Force_y);
   if (Force_z != NULL)
      free (Force_z);

   if (Ncell_bound_x != NULL)
      free (Ncell_bound_x);
   if (Ncell_bound_y != NULL)
      free (Ncell_bound_y);
   if (Ncell_bound_z != NULL)
      free (Ncell_bound_z);

   if (Pbc_x != NULL)
      free (Pbc_x);
   if (Pbc_y != NULL)
      free (Pbc_y);
   if (Pbc_z != NULL)
      free (Pbc_z);

   if (Diam != NULL)
      free (Diam);
   if (X != NULL)
      free (X);
   if (Y != NULL)
      free (Y);
   if (Z != NULL)
      free (Z);

   if (Percentage != NULL)
      free (Percentage);
   if (Perc_values != NULL)
      free (Perc_values);

   if (gen_par_file != NULL)   fclose(gen_par_file);
   if (beg_coord_file != NULL) fclose(beg_coord_file);
   if (res_coord_file != NULL) fclose(res_coord_file);
   if (inp_coord_file != NULL) fclose(inp_coord_file);
   if (res_data_file != NULL)  fclose(res_data_file);
}

/*
****************************************************************

    Ausgeben der x,y,z-Koordinaten und Durchmesser

****************************************************************
*/
void write_sphere_data()
{
   int count= 1000, i;
   char filename[80]= "";

   strcpy(filename, res_coord_name);

   res_coord_file= fopen(filename, "w");

   for (i= 0; i < No_parts; i++)
   {
      fprintf(res_coord_file, "%.14lf\t%.14lf\t%.14lf\t%.14lf\n", *(X+i), *(Y+i), *(Z+i), *(Diam+i)*Din);
      fflush(res_coord_file);
   }

   return;
}

/*
****************************************************************

    Hilfsfunktion Vergleich

****************************************************************
*/
int dblcmp(const void *v1, const void *v2)
{
   double t = *(double *)v1 - *(double *)v2;

   return ((t > 0) ? 1 : -1);
}

/*
****************************************************************

    Hilfsfunktion Zeitausgabe

****************************************************************
*/
void date_time(FILE * pfile_p, long * time_p)
{
   long clock;

   *time_p = time(&clock);
   fprintf(pfile_p, "%s", ctime(&clock));
}

/*
****************************************************************

    Fehlerbehandlung

****************************************************************
*/
int fehlerAusgabe(char *funktion, int fehlerNummer)
{
   int fehlerIndex= (-fehlerNummer) - 1;
	int fehlerart= 1; /* 0: behebbarer Fehler; 1: harter Fehler */

   const char fehlerNachricht[FANZ][60+1]=
      {
      " -1 Fehler in den Argumenten",
      " -2 Speicherplatz konnte nicht bereit gestellt werden",
      " -3 Datei konnte nicht zum Lesen geoeffnet werden",
      " -4 Datensatz konnte nicht gelesen werden",
      " -5 Datei konnte nicht zum Schreiben geoeffnet werden",
      " -6 Datensatz konnte nicht geschrieben werden",
      " -7 Datei konnte nicht geschlossen werden",
      };

   if (fehlerIndex >= 0 && fehlerIndex < FANZ)
   {
      printf("\nFehlermeldung: %s\n", fehlerNachricht[fehlerIndex]);
   }

   if (fehlerart)
   {
      end_run();

      exit(0);
   }
	return fehlerart;
}
