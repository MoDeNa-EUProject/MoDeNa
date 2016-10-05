/*
****************************************************************
  
    Funktionsdeklarationen

****************************************************************
*/

//Anlegen Listen und Dateien
void   make_project(void);
int    read_param(void);
int    set_coordinates(void);
int    list_structure(void);

//Algorithmus
int    packing_loop(void);
int    stepon(void);
int    interaction(void);
void   force_part(int ipart_p, double Rcut);
void   force_all(int ipart_p);
void   motion(void);

//Ausgabe
void   result(int nr);
void   write_sphere_data(void);

//Aufräumen
void   end_run(void);

//Hilfsfunktionen
void   sort(void);
int    dblcmp(const void*, const void*);
void   date_time(FILE * pfile_p, long * time_p);
int    fehlerAusgabe(char *funktion, int fehlerNummer);

//Verteilungen
void   setDistribution(int,double*);






