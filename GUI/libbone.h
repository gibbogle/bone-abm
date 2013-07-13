#ifndef LIBBONE_H
#define LIBBONE_H

#ifdef __cplusplus
extern "C" {
#endif
//
//
void execute(char *, int *);
void simulate_step(int *);
void terminate_run(int *);
void get_dimensions(int *,int *,int *,int *);
void get_scene(int *, float *, int *, int *, int *, float *, int *, float *, int *, int *);
void get_summary(int *);
//
//
#ifdef __cplusplus
}
#endif

#endif // LIBBONE_H
