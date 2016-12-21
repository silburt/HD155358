/**
 * HD155358 - Part 2 - Add Planetesimal disk to resonant planets.
 *
 * This example integrates the HD155358 system + N planetesimal disk, with the hopes
 * of reproducing the libration amplitudes of the system. In this part of the evolution
 * the planets are now in resonance, and a planetesimal disk is added. Subsequent 
 * evolution is analyzed.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include "rebound.h"
#include <time.h>

void heartbeat(struct reb_simulation* r);
double calc_a(struct reb_simulation* r, int index);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);

double E0;
int N_prev;
char output_name[100] = {0};
time_t t_ini;
double tout = 0;

//binary save
char binary_output_name[100] = {0};
double binary_output_time;

int main(int argc, char* argv[]){
    char binary[100] = {0}; strcat(binary, argv[1]); strcat(binary,".bin");
    struct reb_simulation* r = reb_create_simulation_from_binary(binary);
    int N_planetesimals = atoi(argv[2]);
    int seed = atoi(argv[3]);
    srand(seed);
    strcat(output_name,argv[1]); strcat(output_name,"planetesimals_Np"); strcat(output_name,argv[2]); strcat(output_name,"sd"); strcat(output_name,argv[3]);
    
	// Simulation Setup
	r->integrator	= REB_INTEGRATOR_HERMES;
    r->heartbeat	= heartbeat;
    r->ri_hermes.hill_switch_factor = 3;
    r->ri_hermes.solar_switch_factor = 20.;
    r->testparticle_type = 1;
    //r->gravity_ignore_10 = 0; //Use if created binary with WHFAST but now using !WHFAST.
    double tmax = 1e6;
    tout = r->t;
    
    // Collisions
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    r->track_energy_offset = 1;
    r->collision_resolve_keep_sorted = 1;
    
    // Boundaries
    r->boundary	= REB_BOUNDARY_OPEN;
    const double boxsize = 5;
    reb_configure_box(r,boxsize,2,2,1);
    
    // Planetesimal disk parameters (Planets already added)
    double total_disk_mass = r->particles[1].m/10.;
    double planetesimal_mass = total_disk_mass/N_planetesimals;
    printf("%e,%e\n",total_disk_mass,planetesimal_mass);
    double amin = calc_a(r, 1) - 0.5, amax = calc_a(r, 2) + 0.5;
    double powerlaw = 0;
    
    // Generate Planetesimal Disk
    struct reb_particle star = r->particles[0];
    while(r->N<(N_planetesimals + r->N_active)){
		struct reb_particle pt = {0};
		double a    = reb_random_powerlaw(amin,amax,powerlaw);
        double e    = reb_random_rayleigh(0.001);
        double inc  = reb_random_rayleigh(0.001);
        double Omega = reb_random_uniform(0,2.*M_PI);
        double apsis = reb_random_uniform(0,2.*M_PI);
        double phi 	= reb_random_uniform(0,2.*M_PI);
        pt = reb_tools_orbit_to_particle(r->G, star, r->testparticle_type?planetesimal_mass:0., a, e, inc, Omega, apsis, phi);
		pt.r 		= 0.00000934532;
        pt.hash = r->N;
		reb_add(r, pt);
    }

    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    N_prev = r->N;
    
    //binary (temp)
    strcat(binary_output_name, output_name); strcat(binary_output_name, "_t=");
    binary_output_time = r->t;
    
    //naming
    char timeout[200] = {0};
    char info[200] = {0};
    strcat(timeout,output_name); strcat(info,output_name);
    char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,output_name); strcat(syss,"*");
    system(syss);
    strcat(output_name,".txt");
    time_t t_ini = time(NULL);
    struct tm *tmp = gmtime(&t_ini);
    
    //info output
    strcat(info,"_info.txt");
    FILE* out1 = fopen(info,"w");
    fprintf(out1, "Simulation Details:\n");
    int coll_on = 0; if(r->collision_resolve == reb_collision_resolve_merge) coll_on =1;
    fprintf(out1, "\nSetup Parmaeters:\nHSF=%.2f, RSF=%.1f, dt=%e, tmax=%e, collisions_on=%d\n",r->ri_hermes.hill_switch_factor,r->ri_hermes.solar_switch_factor,r->dt,tmax,coll_on);
    fprintf(out1, "\nPlanet(s):\n");
    for(int i=1;i<r->N_active;i++){
        struct reb_particle p = r->particles[i];
        fprintf(out1,"Planet %d: m=%e, r=%e, a=%e\n",i,p.m,p.r,calc_a(r,i));
    }
    fprintf(out1, "\nPlanetesimal Disk:\nNumber of planetesimals=%d\ntotal mass of planetesimal disk=%e\nmass of each planetesimal=%e\nsemi-major axis limits of planetesimal disk: a_min=%f, amax_pl=%f\npowerlaw of planetesimal disk=%.2f\n",N_planetesimals,total_disk_mass,planetesimal_mass,amin,amax,powerlaw);
    fclose(out1);
    
    // Integrate!
    reb_integrate(r, tmax);
    
    //time output
    time_t t_fini = time(NULL);
    struct tm *tmp2 = gmtime(&t_fini);
    double time = t_fini - t_ini;
    strcat(timeout,"_elapsedtime.txt");
    FILE* outt = fopen(timeout,"w");
    fprintf(outt,"\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
    fclose(outt);
    printf("\nSimulation complete. Elapsed simulation time is %.2f s. \n\n",time);
}

void heartbeat(struct reb_simulation* r){
    if (tout <r->t){
        //tout += 0.01;
        tout += 25;
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        int N_mini = 0;
        if (r->ri_hermes.mini_active){
            N_mini = r->ri_hermes.mini->N;
        }
        
        FILE* f = fopen(output_name, "a");
        fprintf(f,"%e,%e,%d,%d,",r->t,relE,r->N,N_mini);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    if (reb_output_check(r, 100.*r->dt)){
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
    }
    
    if(r->N < N_prev){
        N_prev = r->N;
        double E = reb_tools_energy(r);
        reb_move_to_com(r);
        r->energy_offset += E - reb_tools_energy(r);
    }
    
    //output binary, temp
    if(binary_output_time < r->t){
        char out_time[10] = {0}; sprintf(out_time,"%.0f",r->t);
        char out[200] = {0}; strcat(out, binary_output_name); strcat(out, out_time); strcat(out, ".bin");
        reb_output_binary(r, out);
        //if(r->t >= 149999){
        //    binary_output_time += 1e2;
        //} else {
        binary_output_time += 1e5;
        //}
    }
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

void calc_resonant_angles(struct reb_simulation* r, FILE* f){
    double e[3] = {0}; //Hardcoding, probably should change in the future.
    double a[3] = {0};
    double omega[3] = {0};
    double lambda[3] = {0};
    struct reb_particle com = reb_get_com(r);
    for(int i=1;i<r->N_active;i++){
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
    
    fprintf(f,"%f,%f,%f,%f,%f,%f,%f\n",a[1],e[1],a[2],e[2],phi,phi2,phi3);
    
}
