using namespace std;
#ifndef INOUT_H
#define INOUT_H

void readParameters(string, string &, string &, string &, string &, string &, \
    string &);
int ***allocateFromVTK(string, int ***);
void importFromVTK(string, int ***);
void saveToVTK(const char*, int ***);
void saveToDX(const char*, int ***);
void saveToGnuplot(string, int, int, double **, int **);
void saveDescriptors(string, double, double);
void saveParameters(string, double);

#endif
