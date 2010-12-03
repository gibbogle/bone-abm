#ifndef TRANSFER_H
#define TRANSFER_H

#include <QMutex.h>

extern int showingVTK;

#define MAX_MONO 10000
#define MAX_CAP 1000
#define MAX_PIT 10000
#define MAX_CLAST 100

extern int VTKbuffer[100];
extern int mono_list[5*MAX_MONO];
extern int nmono_list;
extern float cap_list[7*MAX_CAP];
extern int ncap_list;
extern float pit_list[4*MAX_PIT];
extern int npit_list;
extern float clast_list[MAX_CLAST];
extern int nclast_list;
extern QMutex mutex1, mutex2;

extern int summaryData[100];
extern int NX, NY, NZ, NBY;
extern int nt_vtk;
extern bool leftb;

#endif // TRANSFER_H
