/**
 * HD155358 - Stability
 *
 * Using simulation archive, simulate a given set of parameters for X years into the future. 
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
#define MIN(a, b) ((a) > (b) ? (b) : (a))       ///< Returns the maximum of a and b

void heartbeat(struct reb_simulation* r);
void calc_resonant_angles(struct reb_simulation* r, FILE* f);

double E0;
char filename[200] = {0};
char binary_out[100] = {0};
double* omega;
double* lambda;
double* a;
double* e;
double* MA;
double log_constant, tlog_output;
clock_t start;

int main(int argc, char* argv[]){
    strcat(filename,argv[1]);
    char simname[300] = {0}; strcat(simname,argv[1]); strcat(simname,".bin");
    double tmax = atof(argv[2]);
    struct reb_simulation* r = reb_create_simulation_from_simulationarchive(simname);
    
    if (r==NULL){
        printf("No simulation archive found. Creating new simulation.\n");
        r= reb_create_simulation();
        
        double mJ = 0.00095424836;  //mass of Jupiter
        double m0 = 0.92;
        
        double m1sini = atof(argv[3]);      //units of jupiter mass
        double m2sini = atof(argv[4]);      //units of jupiter mass
        double a1 = atof(argv[5]);
        double a2 = atof(argv[6]);
        double h1 = atof(argv[7]);
        double h2 = atof(argv[8]);
        double k1 = atof(argv[9]);
        double k2 = atof(argv[10]);
        double lambda1 = atof(argv[11]);
        double lambda2 = atof(argv[12]);
        double sini = atof(argv[13]);
        
        // Simulation Setup
        r->integrator	= REB_INTEGRATOR_WHFAST;
        r->dt = 2*M_PI*sqrt(a1*a1*a1/m0)/200;
        r->exit_min_distance = MIN(a1*pow(m1sini*mJ/sini/m0/3.,1./3.) , a2*pow(m2sini*mJ/sini/m0/3.,1./3.));
        r->exit_max_distance = 10;
        r->simulationarchive_interval = 2.*M_PI*tmax/1000; // 1000 outputs
        r->ri_whfast.safe_mode = 0;
        
        // Star
        struct reb_particle star = {0};
        star.m 		= m0;
        star.r		= 0.005;        // Radius of particle is in AU!
        star.hash = 0;
        reb_add(r, star);
        
        // Planet 1
        {
            struct reb_particle p = {0};
            p = reb_tools_pal_to_particle(r->G, star, m1sini*mJ/sini, a1, lambda1, k1, h1, 0, 0);
            p.hash = r->N;
            reb_add(r, p);
        }
        
        //Planet 2
        {
            struct reb_particle p = {0};
            p = reb_tools_pal_to_particle(r->G, star, m2sini*mJ/sini, a2, lambda2, k2, h2, 0, 0);
            p.hash = r->N;
            reb_add(r, p);
        }
        
        r->N_active = r->N;
        
        //naming
        char syss[100] = {0}; strcat(syss,"rm -v "); strcat(syss,filename); strcat(syss,"*"); system(syss);
    } else {
        printf("Found simulation archive. Loaded snapshot at t=%.16f.\n",r->t);
    }
    
    //ini
    reb_move_to_com(r);
    E0 = reb_tools_energy(r);
    r->heartbeat	= heartbeat;
    
    //arrays
    a = calloc(sizeof(double),r->N);
    e = calloc(sizeof(double),r->N);
    omega = calloc(sizeof(double),r->N);
    lambda = calloc(sizeof(double),r->N);
    MA = calloc(sizeof(double),r->N);
    
    //naming
    r->simulationarchive_filename = simname;
    strcat(filename,".csv");
    
    //output frequency
    int n_output = 5000;
    log_constant = pow(tmax + 1, 1./(n_output - 1));
    tlog_output = r->t + r->dt;
    
    // Integrate!
    start = clock();
    reb_integrate(r, tmax);
    
    free(a); free(e); free(omega); free(lambda); free(MA);
    
}

void heartbeat(struct reb_simulation* r){
    //outputs
    if(r->t > tlog_output){
        tlog_output = r->t*log_constant;
        
        reb_integrator_synchronize(r);
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        
        FILE* f = fopen(filename, "a");
        fprintf(f,"%e,%e,",r->t,relE);
        calc_resonant_angles(r,f);
        fclose(f);
    }
    
    if (reb_output_check(r, 10000.*r->dt)){
        reb_integrator_synchronize(r);
        double E = reb_tools_energy(r);
        double relE = fabs((E-E0)/E0);
        reb_output_timing(r, 0);
        printf("%e",relE);
        
        //check that not going overtime
        float hours = (clock() - start)*1.0 / CLOCKS_PER_SEC / 3600;
        if(hours > 47.75){
            printf("\n***Max time elapsed for %s. Exiting cleanly.***\n",filename);
            exit(0);
        }
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
        MA[i] = E - e[i]*sin(E);
        lambda[i] = MA[i] + omega[i];
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
    
    fprintf(f,"%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",a[1],e[1],lambda[1],omega[1],MA[1],a[2],e[2],lambda[2],omega[2],MA[2],phi,phi2,phi3);
    
}


