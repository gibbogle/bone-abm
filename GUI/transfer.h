#ifndef TRANSFER_H
#define TRANSFER_H

#include <QMutex.h>

extern int showingVTK;

extern int VTKbuffer[100];
extern int mono_list[5*10000];
extern int nmono_list;
extern float cap_list[7*100];
extern int ncap_list;
extern float pit_list[4*1000];
extern int npit_list;
extern QMutex mutex1, mutex2;

extern int summaryData[100];
extern int NX, NY, NZ, NBY;
extern int nt_vtk;

#endif // TRANSFER_H
