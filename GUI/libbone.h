#ifndef LIBBONE_H
#define LIBBONE_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
#ifdef __GFORTRAN_DLL__
//void __omp_main_mod_MOD_execute(int *, char *, char *, char *, char *, int, int, int, int);
//void __bone_mod_MOD_execute(int *);
void __bone_mod_MOD_execute(char *, int);
#else
//void EXECUTE(int *, char *, int, char *, int, char *, int, char *, int);
//void EXECUTE(int *);
void EXECUTE(char *, int);
#endif
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBBONE_H
