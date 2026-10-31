/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef DRSOLVE_PROFILE_H
#define DRSOLVE_PROFILE_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/* Inclusive monotonic wall times; nested entries must not be added together.
 * Configure the environment before solving. No shared mutable counters. */
static inline int dr_profile_enabled(void) {
    const char *s=getenv("DRSOLVE_PROFILE_DETAIL");
    return s && *s && strcmp(s,"0");
}
static inline double dr_profile_now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t);
    return (double)t.tv_sec+1e-9*t.tv_nsec;
}
static inline double dr_profile_start(void) {
    return dr_profile_enabled()?dr_profile_now():0;
}
static inline void dr_profile_end(const char *name,double start) {
    if(start)fprintf(stderr,"Profile wall %s %.9f\n",name,dr_profile_now()-start);
}
#endif
