/**
 * HD155358 - Part 1 - Migrate into resonance.
 *
 * This example integrates the HD155358 system + N planetesimal disk, with the hopes 
 * of reproducing the libration amplitudes of the system. In this part of the evolution
 * I am simply migrating the planets into resonance (planetesimal disk is added later).
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include <time.h>
#define MAX(a, b) ((a) < (b) ? (b) : (a))       ///< Returns the maximum of a and b

void heartbeat(struct reb_simulation* r);
void heartbeat_getRV(struct reb_simulation* r);
void migration_forces(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);
double calc_a(struct reb_simulation* r, int index);
double calc_P(struct reb_simulation* r, int index);

double E0;
char output_name[100] = {0};
time_t t_ini;
double* tau_a; 	/**< Migration timescale in years for all particles */
double* tau_e; 	/**< Eccentricity damping timescale in years for all particles */
double* omega;
double* lambda;
double* a;
double* e;
double mig_time, dispersal_time, dispersal_rate, mig_rate, K1, K2, dispersal_fac;
int mig_index, print_disk_dispersal=0, output_rate;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    
    double m1 = atof(argv[1]);      //units of jupiter mass
    double m2 = atof(argv[2]);      //units of jupiter mass
    double sini = atof(argv[3]);
    mig_rate = atof(argv[4]);
    K1 = atof(argv[5]);              //Lee & Peale (2002) K.
    K2 = atof(argv[6]);
    srand(atoi(argv[7]));
    strcat(output_name,argv[8]);
    
    r->integrator	= REB_INTEGRATOR_WHFAST;
    
	// Simulation Setup
    r->heartbeat	= heartbeat;
    r->additional_forces = migration_forces;
    r->force_is_velocity_dependent = 1;
    
    // Boundaries
    //r->boundary	= REB_BOUNDARY_OPEN;
    //const double boxsize = 30;
    //reb_configure_box(r,boxsize,2,2,1);
    
	// Star
	struct reb_particle star = {0};
	star.m 		= 0.92;
    star.r		= 0.005;        // Radius of particle is in AU!
    star.hash = 0;
	reb_add(r, star);
    
    double mJ = 0.00095424836;  //mass of Jupiter
    //double a1 = 0.6, a2=1.02;
    double a1 = 1, a2=1.8;
    r->dt = 2*M_PI*sqrt(a1*a1*a1/star.m)/50;
    // Planet 1
    {
        double a=a1, m=m1*mJ/sini, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, 0, inc, 0, 0, reb_random_uniform(0,2.*M_PI));
        p.r = 0.000467;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    //Planet 2
    {
        double a=a2, m=m2*mJ/sini, inc=reb_random_normal(0.00001);
        struct reb_particle p = {0};
        p = reb_tools_orbit_to_particle(r->G, star, m, a, 0, inc, 0, 0, reb_random_uniform(0,2.*M_PI));
        p.r = 0.000467;
        p.hash = r->N;
        reb_add(r, p);
    }
    
    r->N_active = r->N;
    
    //Migration times and rates
    mig_index = 2;                                                          //index of migrating planet
    mig_time = MAX(5*mig_rate,3e3*calc_P(r,1));                             //migration time
    dispersal_time = mig_time;                                              //1e4 orbital periods of inner planet
    dispersal_rate = pow(5e7/mig_rate + 1, 1./(dispersal_time/r->dt - 1));  //rate of disk dispersal
    dispersal_fac = 1;
    double tmax = 1.25*mig_time + dispersal_time;
    
    //Migraiton arrays
    tau_a = calloc(sizeof(double),r->N);
    tau_e = calloc(sizeof(double),r->N);
    tau_a[mig_index] = 2.*M_PI*mig_rate/sqrt(calc_a(r,mig_index));  //tau_a = a/dot(a), scale-free
    tau_e[1] = tau_a[mig_index]/K1;                                 //tau_e = e/dot(e)
    tau_e[2] = tau_a[mig_index]/K2;
    printf("\nmig_time=%e, dispersal_time=%e, dispersal_rate=%.10e\n",mig_time,dispersal_time,dispersal_rate);
    printf("\ntau_a=%e, tau_e1=%e, tau_e2=%e\n",tau_a[mig_index],tau_e[1],tau_e[2]);
    
    //Other arrays
    a = calloc(sizeof(double),r->N);
    e = calloc(sizeof(double),r->N);
    omega = calloc(sizeof(double),r->N);
    lambda = calloc(sizeof(double),r->N);
    output_rate = (int)tmax/5000;
    
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    
    //naming
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*"); system(syss);
    char binary[120] = {0}; strcat(binary,output_name); strcat(binary,".bin");
    //char RVout[120] = {0}; strcat(RVout, output_name); strcat(RVout, "_RV.txt");
    
    strcat(output_name,".txt");
    
    // Integrate!
    reb_integrate(r, tmax);
    
    if(r->N == 3) reb_output_binary(r, binary);
    
}

double tout = 1;
void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        //tout *= 1.01;
        tout += output_rate;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        
        FILE* f = fopen(output_name, "a");
        fprintf(f,"%e,%e,%d,%e,%e,%e,%e,%e,",r->t,relE,r->N,mig_rate,K1,K2,mig_time,dispersal_time+mig_time);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    if (reb_output_check(r, 500.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
    
    //double check timestep
    if(r->dt > calc_P(r,1)/30){
        r->dt /= 4;
        printf("\nreduced timestep\n");
    }
}

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<r->N;i++){
        struct reb_particle p = r->particles[i];
        const double mu = r->G*(r->particles[0].m + p.m);
        const double dvx = p.vx-com.vx;
        const double dvy = p.vy-com.vy;
        const double dvz = p.vz-com.vz;
        const double dx = p.x-com.x;
        const double dy = p.y-com.y;
        const double dz = p.z-com.z;
        
        const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
        const double d = sqrt ( dx*dx + dy*dy + dz*dz );
        const double vr = (dx*dvx + dy*dvy + dz*dvz)/d;
        const double ex = 1./mu*( (v*v-mu/d)*dx - d*vr*dvx );
        const double ey = 1./mu*( (v*v-mu/d)*dy - d*vr*dvy );
        const double ez = 1./mu*( (v*v-mu/d)*dz - d*vr*dvz );
        e[i] = sqrt( ex*ex + ey*ey + ez*ez );   // eccentricity
        a[i] = -mu/(v*v - 2.*mu/d);
        const double rdote = dx*ex + dy*ey + dz*ez;
        const double cosf = rdote/(e[i]*d);
        
        omega[i] = atan2(ey,ex);
        if(ey < 0.) omega[i] += 2*M_PI;
        double cosE = (a[i] - d)/(a[i]*e[i]);
        double E;
        if(cosf > 1. || cosf < -1.){
            E = M_PI - M_PI*cosE;
        } else {
            E = acos(cosE);
        }
        if(vr < 0.) E = 2.*M_PI - E;
        double MA = E - e[i]*sin(E);
        lambda[i] = MA + omega[i];
    }
    double phi = 0, phi2 = 0, phi3 = 0;
    phi = 2.*lambda[2] - lambda[1] - omega[1];
    phi2 = 2.*lambda[2] - lambda[1] - omega[2];
    phi3 = omega[1] - omega[2];
    while(phi >= 2*M_PI) phi -= 2*M_PI;
    while(phi < 0.) phi += 2*M_PI;
    while(phi2 >= 2*M_PI) phi2 -= 2*M_PI;
    while(phi2 < 0.) phi2 += 2*M_PI;
    while(phi3 >= 2*M_PI) phi3 -= 2*M_PI;
    while(phi3 < 0.) phi3 += 2*M_PI;
    
    fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",a[1],e[1],a[2],e[2],phi,phi2,phi3,r->particles[1].m,r->particles[2].m,tau_a[1]*a[mig_index],tau_e[1]*a[mig_index],tau_a[2]*a[mig_index],tau_e[2]*a[mig_index]);
    
}

double calc_a(struct reb_simulation* r, int index){
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = reb_get_com(r);
    struct reb_particle p = particles[index];
    const double mu = r->G*(com.m + p.m);
    const double dvx = p.vx-com.vx;
    const double dvy = p.vy-com.vy;
    const double dvz = p.vz-com.vz;
    const double dx = p.x-com.x;
    const double dy = p.y-com.y;
    const double dz = p.z-com.z;
    
    const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
    const double d = sqrt(dx*dx + dy*dy + dz*dz);    //distance
    const double dinv = 1./d;
    const double a = -mu/(v2 - 2.*mu*dinv);
    
    return a;
}

double calc_P(struct reb_simulation* r, int index){
    double a = calc_a(r, index);
    double Ms = r->particles[0].m;
    return 2*M_PI*sqrt(a*a*a/Ms);  //units of yr/2pi
}

void migration_forces(struct reb_simulation* r){
    if(r->t < (mig_time + dispersal_time)){
        //const double G = r->G;
        const int N = r->N;
        struct reb_particle* const particles = r->particles;
        struct reb_particle com = reb_get_com(r);
        for(int i=1;i<N;i++){
            if (tau_e[i]!=0||tau_a[i]!=0){
                struct reb_particle* p = &(particles[i]);
                const double dvx = p->vx-com.vx;
                const double dvy = p->vy-com.vy;
                const double dvz = p->vz-com.vz;
                
                if (tau_a[i]!=0){ 	// Migration
                    p->ax -=  dvx/(2.*tau_a[i]);
                    p->ay -=  dvy/(2.*tau_a[i]);
                    p->az -=  dvz/(2.*tau_a[i]);
                }
                
                if (tau_e[i]!=0){ 	// Eccentricity damping
                    const double dx = p->x-com.x;
                    const double dy = p->y-com.y;
                    const double dz = p->z-com.z;
                    
                    //Papalouzou & Larwood (2000)
                    const double rinv = 1/sqrt( dx*dx + dy*dy + dz*dz );
                    const double vr = (dx*dvx + dy*dvy + dz*dvz)*rinv;
                    const double term = -2.*vr*rinv/tau_e[i];
                    
                    p->ax += term*dx;
                    p->ay += term*dy;
                    p->az += term*dz;
                    
                    /*
                    //Lee & Peale (2002)
                    const double mu = r->G*(com.m + p->m);
                    const double hx = dy*dvz - dz*dvy;
                    const double hy = dz*dvx - dx*dvz;
                    const double hz = dx*dvy - dy*dvx;
                    const double h = sqrt ( hx*hx + hy*hy + hz*hz );
                    const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
                    const double r = sqrt ( dx*dx + dy*dy + dz*dz );
                    const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
                    const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
                    const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
                    const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
                    const double e = sqrt( ex*ex + ey*ey + ez*ez );		// eccentricity
                    const double a = -mu/( v*v - 2.*mu/r );			// semi major axis
                    const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
                    const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
                    p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
                    p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
                    p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
                     */
                }
            }
            com = reb_get_com_of_pair(com,particles[i]);
        }
        
        //update migration rate to keep scale-free
        tau_a[mig_index] = 2.*M_PI*mig_rate/sqrt(calc_a(r,mig_index));
        tau_e[1] = tau_a[mig_index]/K1;
        tau_e[2] = tau_a[mig_index]/K2;
        
        //disk_dispersal
        if(r->t > mig_time){
            if(print_disk_dispersal == 0){
                printf("\n **Disk Dispersal Started**\n");
                print_disk_dispersal = 1;
            }
            dispersal_fac *= pow(5e7/mig_rate + 1, 1./(dispersal_time/r->dt - 1));
            tau_a[mig_index] *= dispersal_fac;
            tau_e[1] *= dispersal_fac;
            tau_e[2] *= dispersal_fac;
        }
    } else if (print_disk_dispersal==1){
        printf("\n **Disk Dispersal Finished**\n");
        print_disk_dispersal = 2;
    }
}

//old Get RV signal
/*
 if(get_RV && (r->N == 3)){
 double day_to_yr2pi = 2*M_PI/365;
 double AUyr2ms = 29682.77;
 
 double tmaxRV = 60.0;     //yr/2pi into the future we want to simulate
 int nsteps = 600;
 r->dt = (calc_a(r, 1),1.5)/100;
 
 reb_reset_function_pointers(r);
 E0 = reb_tools_energy(r);
 printf("\nGetting RV signal.\n\n");
 for(int i=0;i<nsteps;i++){
 reb_integrate(r, tmax + i*tmaxRV/nsteps);
 double relE = fabs((reb_tools_energy(r)-E0)/E0);
 FILE* f = fopen(RVout, "a");
 fprintf(f,"%e,%e,%e,%e,%d\n",i*tmaxRV/nsteps,AUyr2ms*r->particles[0].vx*sini, AUyr2ms*r->particles[0].vy*sini,relE,r->N);
 fclose(f);
 
 reb_output_timing(r, 0);
 printf("%e",relE);
 }
 }*/

//int binary_out = 0;
//if(binary_out) reb_output_binary(r, binary);
//printf("\nSimulation complete.\n\n");
