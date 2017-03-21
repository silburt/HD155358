//This is testing whether the pal->particle and particle->pal functions are working correctly in python. To run this, you need to place it in a rebound subdirectory with a Makefile. It's not gonna work here. 

#include "rebound.h"
#include <stdio.h>
#include <stdlib.h>

void heartbeat(struct reb_simulation* r){
    printf("%f\n",r->t);
}

int main(int argc, char* argv[]) {
    struct reb_simulation* r = reb_create_simulation();
    r->dt = 0.1;
    r->heartbeat = heartbeat;
    r->exact_finish_time = 1; // Finish exactly at tmax in reb_integrate(). Default is already 1.
    
    struct reb_particle p1 = {0};
    p1.m = 1.;
    reb_add(r, p1);
    
    {
        struct reb_particle p2 = {0};
        double m=0.00271308265758,a=0.64141678,e=0.17145134,w=2.92629096,M=2.03300709;
        double f=reb_tools_M_to_f(e, M);
        p2 = reb_tools_orbit_to_particle(r->G, p1, m, a, e, 0, 0, w, f);
        reb_add(r, p2);
        
        double h=0, k=0, lambda=0, ix=0, iy=0, a2=0;
        reb_tools_particle_to_pal(r->G,p2,p1,&a2,&lambda,&k,&h,&ix,&iy);
        struct reb_particle pal = reb_tools_pal_to_particle(r->G,p1,m,a2,lambda,k,h,ix,iy);
        printf("\n p_inner: %e,%e,%e,%e,%e,%e",a2,lambda,k,h,ix,iy);
        printf("\n pal: %e,%e,%e,%e,%e,%e\n",pal.x-p2.x,pal.y-p2.y,pal.z-p2.z,pal.vx-p2.vx,pal.vy-p2.vy,pal.vz-p2.vz);
    }
    
    {
        struct reb_particle p3 = {0};
        double m=0.00245564203971,a=1.02299252,e=0.103073,w=5.5671592,M=1.86941606;
        double f=reb_tools_M_to_f(e, M);
        p3 = reb_tools_orbit_to_particle(r->G, p1, m, a, e, 0, 0, w, f);
        reb_add(r, p3);
        
        double h=0, k=0, lambda=0, ix=0, iy=0, a2=0;
        reb_tools_particle_to_pal(r->G,p3,p1,&a2,&lambda,&k,&h,&ix,&iy);
        struct reb_particle pal = reb_tools_pal_to_particle(r->G,p1,m,a2,lambda,k,h,ix,iy);
        printf("\n p_outer: %e,%e,%e,%e,%e,%e",a2,lambda,k,h,ix,iy);
        printf("\n pal: %e,%e,%e,%e,%e,%e\n",pal.x-p3.x,pal.y-p3.y,pal.z-p3.z,pal.vx-p3.vx,pal.vy-p3.vy,pal.vz-p3.vz);
        
    }
    
}
