using namespace std;
#ifndef INOUT_H
#define INOUT_H

void readParameters(string, string &, string &, string &, string &, string &);
int ***allocateFromVTK(string, int ***);
void importFromVTK(string, int ***);
void saveToVTK(const char*, int ***);
void saveToDX(const char*, int ***);
void saveToGnuplot(string, int, int, float **, int **);
void saveDescriptors(string, double, double);

#endif
