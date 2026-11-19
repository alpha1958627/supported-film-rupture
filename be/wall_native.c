// v4: y-index-corrected wall MG with consistent multiplicative noise flux.
#include "poisson.h"
#include "run.h"

#include "grid/tree.h"



#include <sys/stat.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <math.h>

face vector prescribed_noise[];
double noise_phi=.001;
static int mode=0, wn=128, started=0;
static double prefactor=1.;
static double attraction;
static double current_dt;
#include <omp.h>
void _init_solver(void);
int outer_count, mg_cycles;
double outer_error, mg_res;

static double Lambda(double x) {
  double z=x*x;
  if(fabs(x)<.2) return z*z*(1.5+z*(-1./6.+z*(1./40.+z*(-1./280.+z*(1./2160.+z*(-1./18480.+z/174720.))))));
  return 1.-3*z-exp(-z)*(1.-2*z)+2*sqrt(pi)*z*x*erf(x);
}
static double Upsilon(double x){return expm1(-x*x)+sqrt(pi)*x*erf(x);}
static double Theta(double x){return fabs(x)<1e-3 ? x*x*(.5+x*(-1./6.+x*(1./24.+x*(-1./120.+x/720.)))) : x+expm1(-x);}
static void stochastic_flux(double hf,double oldhf,double input,double dx,double dt,double *value,double *deriv) {
  if(mode==0) {
    *value=sqrt(2*noise_phi*prefactor)*pow(hf,1.5)*input;
    *deriv=1.5*sqrt(2*noise_phi*prefactor)*sqrt(hf)*input;
  } else {
    double sigma=.05,tau=.01,h=hf;
    double ratio=sqrt(pow(sigma,4)*tau*Lambda(h/sigma)*Upsilon(dx/sigma)*Theta(dt/tau)/(2*pi*h*h*h*dx*dt));
    ratio*=sqrt(sigma*Upsilon(dx/sigma)/(sqrt(pi)*dx));
    *value=sqrt(2*noise_phi*prefactor)*pow(h,1.5)*input/(dx*sqrt(dt))*ratio;
    double x=h/sigma;
    *deriv=(*value)*3*x*Upsilon(x)/(sigma*Lambda(x)); // Lambda'(x)=6*x*Upsilon(x).
  }
}
struct Layer {
    scalar f;
    double dt;
    double m;
    double n;
    double dd;
    face vector D;
    face vector Noise1;
    face vector beta;
    face vector gamma;
    scalar lambda;
    double tolerance;
  };
  
  static double residual_layer (scalar * al, scalar * bl,
                    scalar * resl, void * data)
  {
    scalar h = al[0], b = bl[0], res = resl[0];
    struct Layer * p = (struct Layer *) data;
    (const) face vector alpha = p->D;
    (const) face vector Noise1 = p->Noise1;
    (const) face vector beta = p->beta;
    (const) face vector gamma = p->gamma;
    (const) scalar lambda = p->lambda;
    double mytolerance = p->tolerance;
    double maxres = 0.;
    
    face vector g[], hf[] ,No[];
  
    foreach_face()
      hf.x[] = face_gradient_x (h, 0);
    boundary_flux ({hf});
  
    scalar pres[];
    foreach() {
      
      pres[] = 0.;
      foreach_dimension()
        pres[] -= (hf.x[1]- hf.x[0])/Delta;
    }
    boundary({pres});
    
    foreach_face() {
      
      g.x[] = alpha.x[]*face_gradient_x (pres, 0);
      hf.x[] = beta.x[]*face_gradient_x (h, 0) +
        gamma.x[]*(h[]+h[-1])/2.;
      No.x[]=Noise1.x[]*(h[]+h[-1])/2.; // v4: same face-height factor as relax_layer
    }
    boundary_flux ({g, hf,No});
    foreach (reduction(max:maxres)) {
      res[] = b[] - lambda[]*h[];
      foreach_dimension()
          res[] -= (g.x[1] - g.x[])/Delta + (hf.x[1] - hf.x[])/Delta +(No.x[1] - No.x[])/Delta;
      if (fabs (res[]) > maxres)
        maxres = fabs (res[]);
    }
    boundary (resl);
    return maxres;
  }
  static void relax_layer (scalar * al, scalar * bl, int l, void * data)
  {
    scalar h = al[0], b = bl[0];
    struct Layer * p = (struct Layer *) data;
    (const) face vector alpha = p->D;
    (const) face vector Noise1 = p->Noise1;
    (const) face vector beta = p->beta;
    (const) face vector gamma = p->gamma;
    (const) scalar lambda = p->lambda;
  
    
  #if JACOBI
    scalar c[];
  #else
    scalar c = h;
  #endif
  
    foreach_level_or_leaf (l) {
      double n =  sq(sq(Delta))*b[], d =  lambda[]*sq(sq(Delta));
      d -= 5*(alpha.x[] + alpha.x[1] + alpha.y[] + alpha.y[0,1]);
      d += 0.5*Delta*sq(Delta)*(gamma.x[1]-gamma.x[] + gamma.y[0,1]-gamma.y[]);
      d += 0.5*Delta*sq(Delta)*(Noise1.x[1]-Noise1.x[] + Noise1.y[0,1]-Noise1.y[]);
      foreach_dimension() {
        n -= (beta.x[1]*h[1] + beta.x[]*h[-1])*sq(Delta);
        d -= (beta.x[1] + beta.x[])*sq(Delta);
      }
  
      n -= alpha.x[]*(-h[-2] - h[-1,-1] + 5*h[-1] - h[-1,1] + h[0,-1] + h[0,1] + h[1]) +
        alpha.x[1]*(h[-1] + h[0,-1] + h[0,1] - h[1,-1] + 5*h[1] - h[1,1] - h[2]) +
        alpha.y[]*(-h[-1,-1] + h[-1] - h[0,-2] + 5*h[0,-1] + h[0,1] - h[1,-1] + h[1]) +
        alpha.y[0,1]*(h[-1] - h[-1,1] + h[0,-1] + 5*h[0,1] - h[0,2] + h[1] - h[1,1]);
  
      n -= 0.5*Delta*sq(Delta)*(h[1]*gamma.x[1] - h[-1]*gamma.x[] +
                    h[0,1]*gamma.y[0,1] - h[0,-1]*gamma.y[]); 
      n -= 0.5*Delta*sq(Delta)*(h[1]*Noise1.x[1] - h[-1]*Noise1.x[] +
                    h[0,1]*Noise1.y[0,1] - h[0,-1]*Noise1.y[]); 
      c[] = n/d;
    }
    
  #if JACOBI
    foreach_level_or_leaf (l)
      h[] = (h[] + 2.*c[])/3.;
  #endif
    
  #if TRASH
    scalar a1[];
    foreach_level_or_leaf (l)
      a1[] = h[];
    trash ({h});
    foreach_level_or_leaf (l)
      h[] = a1[];
  #endif
  }
  trace
  mgstats layer (scalar myf,double mydt,double mym,double myn,double mydd, int timestep){
  struct Layer p = {0};
  p.tolerance = 1e-3*805306368.; // equation residual scales, not physical BE tolerance
  mg_cycles = 0;
  p.dt=mydt;
  p.f=myf;
  p.m=mym;
  p.n=myn;
  p.dd=mydd;
  
    if (p.dt == 0.) {
      mgstats s = {0};
      return s;
    }
  
    scalar f = p.f;
    scalar b[], fold[], fa[];
    const scalar lambda[] = - 1./p.dt;
    double m = p.m, n = p.n ,dd=p.dd;
    face vector D[], beta[], gamma[] ,Noise1[];
      foreach() {
      fold[] = f[];
      }
  
    double maxres;
    mgstats s;
    int niter = 0, nitermax = 10;
  
    double defaultol = TOLERANCE;
    if (p.tolerance>1e-10)
      TOLERANCE = p.tolerance;
  
    // 在循环外部计算Noise场
  
//   int NN =256;
//   double Noise8[NN][NN];
//   double dl1=100.0/NN;
//   for (int iy = 0; iy < NN; iy++) {
//       // double yy = (iy + 0.5) * dl1;
//       for (int ix = 0; ix < NN; ix++) {
//           // double xx = (ix + 0.5) * dl1;
//           double del=100.0/NN;
//           Noise8[ix][iy] = noise()/(sqrt(del*p.dt*del));
//       }
//   }
  
    do {
        foreach_face()
            beta.x[] = face_gradient_x(f, 0);
        boundary_flux({beta});

        scalar pres[];
        foreach() {
            pres[] = 0.;
            foreach_dimension()
                pres[] -= (beta.x[1] - beta.x[0]) / Delta;
        }
        boundary({pres});
        
        face vector g[];
        face vector NRight[];
        
        foreach_face() {
            double aa = attraction;
            double varphi = noise_phi;
            // Repulsive van der Waals term removed, not merely multiplied by zero.
            double ff = (f[] + f[-1]) / 2.;
          
            double Noiseff = prescribed_noise.x[];

            D.x[] = prefactor*pow(ff, n);
            double noise_value,noise_derivative;
            stochastic_flux(ff,(fold[]+fold[-1])/2.,Noiseff,Delta,p.dt,&noise_value,&noise_derivative);
            Noise1.x[] = noise_derivative;
            beta.x[] = prefactor*(-3 * aa * pow(ff, n - m - 1.));
            gamma.x[] = prefactor*(n * pow(ff, n - 1.) * (pres[] - pres[-1]) / Delta -
                        3 * aa * (n - m - 1.) * pow(ff, n - m - 2.) * (f[] - f[-1]) / Delta);
            g.x[] = gamma.x[] * ff;
            NRight.x[] = noise_derivative*ff-noise_value;
        }
        boundary_flux({g});

        foreach() {
            b[] = -fold[] / p.dt;
            fa[] = f[];
            foreach_dimension()
                b[] += (g.x[1] - g.x[]) / Delta + 1 * (NRight.x[1] - NRight.x[]) / Delta;
        }

        restriction({beta, D, lambda, gamma, Noise1});
        p.lambda = lambda;
        p.D = D;
        p.Noise1 = Noise1;
        p.beta = beta;
        p.gamma = gamma;
        s = mg_solve({f}, {b}, residual_layer, relax_layer, &p);
        mg_cycles += s.i;

        maxres = 0;
        foreach(reduction(max:maxres)) {
            double res = fabs(f[] - fa[]);
            if(res > maxres)
                maxres = res;
        }
        boundary({fa});
        niter++;
    } while(maxres > 1e-8 && niter < nitermax);

    if(p.tolerance > 1e-10)
        TOLERANCE = defaultol;

    outer_count=niter; outer_error=maxres; mg_res=s.resa;
    return s;
}




scalar height[], previous_height[];
static int wrap(int k){return (k%wn+wn)%wn;}
int wall_init(int level,int noise_type) {
  if(level<6 || level>8 || noise_type!=0)return 80;
  if(!started){init_solver();_init_solver();started=1;}
  wn=1<<level;mode=noise_type;omp_set_num_threads(4);
  if(grid)free_grid();
  periodic(right);periodic(top);origin(0,0);size(1);init_grid(wn);
  attraction=sq(2*pi)/3.;return 0;
}
void wall_set_phi(double value){noise_phi=value;}
void wall_noise_probe(double h,double input,double dx,double dt,double *result){stochastic_flux(h,h,input,dx,dt,result,result+1);}
int wall_step(const double *input,const double *noise,double dt,double *output,double *diagnostic) {
  if(!(dt>0))return 80;current_dt=dt;
  foreach(){int ix=wrap((int)floor(x/Delta)),iy=wrap((int)floor(y/Delta));height[]=previous_height[]=input[ix*wn+iy];}
  boundary({height,previous_height});
  foreach_face(x){int ix=wrap((int)floor(x/Delta+1e-8)),iy=wrap((int)floor(y/Delta));prescribed_noise.x[]=noise[ix*wn+iy];}
  foreach_face(y){int ix=wrap((int)floor(x/Delta)),iy=wrap((int)floor(y/Delta+1e-8));prescribed_noise.y[]=noise[wn*wn+ix*wn+iy];}
  boundary((scalar*){prescribed_noise});
  layer(height,dt,3.,3.,5.,0);
  double minimum=HUGE,maximum=-HUGE,mass=0.;int invalid=0;
  foreach(reduction(min:minimum) reduction(max:maximum) reduction(+:mass) reduction(+:invalid)){
    int ix=wrap((int)floor(x/Delta)),iy=wrap((int)floor(y/Delta));output[ix*wn+iy]=height[];
    if(!isfinite(height[]) || height[]<=0)invalid++;
    if(height[]<minimum)minimum=height[];if(height[]>maximum)maximum=height[];mass+=height[]/(wn*wn);
  }
  diagnostic[0]=minimum;diagnostic[1]=maximum;diagnostic[2]=mass;diagnostic[4]=outer_error;diagnostic[5]=mg_res;diagnostic[6]=outer_count;diagnostic[7]=mg_cycles;
  if(invalid){diagnostic[3]=HUGE;return 82;}
  face vector slope[],flux[];scalar pressure[];
  foreach_face()slope.x[]=face_gradient_x(height,0);
  boundary_flux({slope});
  foreach(){pressure[]=0.;foreach_dimension()pressure[]-=(slope.x[1]-slope.x[])/Delta;}
  boundary({pressure});
  foreach_face(){double hf=(height[]+height[-1])/2.,oldhf=(previous_height[]+previous_height[-1])/2.,nv,nd;
    stochastic_flux(hf,oldhf,prescribed_noise.x[],Delta,dt,&nv,&nd);
    flux.x[]=prefactor*(pow(hf,3)*face_gradient_x(pressure,0)-3*attraction/hf*face_gradient_x(height,0))+nv;
  }
  boundary_flux({flux});double error=0.;
  foreach(reduction(max:error)){double r=height[]-previous_height[];foreach_dimension()r-=dt*(flux.x[1]-flux.x[])/Delta;if(fabs(r)>error)error=fabs(r);}
  diagnostic[3]=error;return isfinite(error) ? 0 : 83;
}
int main(void){return 0;}
