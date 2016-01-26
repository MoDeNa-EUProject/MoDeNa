using namespace std;
#ifndef INOUT_H
#define INOUT_H

void readParameters(string, string &, string &, string &, string &, string &, \
    string &);
int ***allocateFromVTK(string, int ***, bool);
void importFromVTK(string, int ***, bool);
void saveToVTK(const char*, int ***, bool);
void saveToDX(const char*, int ***, bool);
void saveToGnuplot(string, int, int, double **, int **, bool);
void saveDescriptors(string, double, double, bool);
void saveParameters(string, double, bool);

#endif
