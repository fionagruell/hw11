#include <iostream>
#include <fstream>
#include <complex>
#include <sstream>
//-----------------------------------
using namespace std;
//-----------------------------------
typedef complex<double> cmplx;
//-----------------------------------
void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin);


void writeToFile(const cmplx* const v, const string s, const double dx,
         const int Nx, const double xmin, const double alpha,
         const double lambda, const double omega, const double t);

void step(cmplx* const psi0, cmplx* const psi1,
          const double dt, const double dx,
          const double omega,const int Nx, const double xmin, const double k,double*v);
//-----------------------------------
int main(){
  
	const int Nx =300 ;
	const double xmin =-40. ;
  const double xmax = 40.;
	const double Tend = 10.0*M_PI;
	const double dx = (xmax-xmin)/(Nx-1);
	const double dt =  1.*dx ;
  double t = 0;
	const int Na = 10;
	int Nk = int(Tend / Na / dt + 0.5);
	double v[Nx];
	
	const double lambda = 10;
  const double omega = 0.2;
  const double k =omega*omega;
  const double alpha = sqrt(omega);
 

  stringstream strm;
  
  
	cmplx* psi0 	= new cmplx[Nx];
	cmplx* psi1 	= new cmplx[Nx];
	cmplx* h 	= new cmplx[Nx]; //hilfsvektor

	init(psi0, alpha, lambda, dx, dt, Nx, xmin);

	writeToFile(psi0,"psi_0", dx,Nx,xmin, alpha, lambda, omega,t);



	for (int i = 1; i <= Na; i++) {
		for (int j = 1; j <= Nk-1; j++) {
		  

	step(psi0,psi1,dt,dx,omega,Nx,xmin,k,v);
	h=psi0;
	psi0=psi1;
	psi1=h;
		  
         t+=dt;
		}
		strm.str("");
		strm << "psi_" << i;
		writeToFile(psi0,strm.str(), dx,Nx,xmin, alpha, lambda, omega,t);
	}
  cout << "t = " << t << endl;
  delete[] psi0;
  delete[] psi1;
  delete[] h;
  

	return 0;
}
//-----------------------------------
void step(cmplx* const psi0, cmplx* const psi1,
          const double dt, const double dx,
          const double omega, const int Nx,const double xmin, const double k, double* v)

{
  cmplx* d= new cmplx[Nx];
  cmplx* a= new cmplx[Nx];
  cmplx* acon= new cmplx[Nx];
  cmplx* dcon= new cmplx[Nx];
  cmplx* aconu= new cmplx[Nx];
  cmplx* acono= new cmplx[Nx];
  

double x;
  for(int i=0; i<Nx; i++){
		  x = xmin + i * dx;
		 v[i]=k*x*x/2.;
  }

  

  for(int i=0;i<Nx;i++) d[i] = 1.+cmplx(0.0,(dt/(2.*dx*dx)+(dt/2.0)*v[i]));
  for(int i=0;i<Nx;i++) dcon[i] = conj(1.+cmplx(0.0,(dt/(2.*dx*dx)+(dt/2.0)*v[i])));
  for(int i=0;i<Nx;i++) a[i] = cmplx(0.0,-dt/(4.0*dx*dx));
  for(int i=0;i<Nx;i++) aconu[i] = conj(cmplx(0.0,-dt/(4.0*dx*dx)));
  for(int i=0;i<Nx;i++) acono[i] = conj(cmplx(0.0,-dt/(4.0*dx*dx)));
  
    
  for(int i=1;i<Nx;i++){
     d[i]-=a[i]/d[i-1]*a[i-1];
     dcon[i]-=a[i]/d[i-1]*acon[i-1];
     psi0[i]-=a[i]/d[i-1]*psi0[i-1];
     aconu[i]-=a[i]/d[i-1]*aconu[i-1];
     acono[i]-=a[i]/d[i-1]*acono[i-1];
  }
  
   
   psi1[Nx-1]=(dcon[Nx-1]*psi0[Nx-1]+aconu[Nx-1]*psi0[Nx-2])/(d[Nx-1]);
   for(int i=Nx-2; i>=1;i--){
   psi1[i]=(aconu[i-1]*psi0[i-1]+dcon[i]*psi0[i]+acono[i+1]*psi0[i+1]-a[i+1]*psi1[i+1])/(d[i]);
   }
   psi1[0]=(dcon[0]*psi0[0]+acono[1]*psi0[1]-a[1]*psi1[1])/(d[0]);
   
   delete[] a;
   delete[] acon;
   delete[] d;
   delete[] dcon;
   delete[] aconu;
   delete[] acono;

}

//-----------------------------------
void writeToFile(const cmplx* const v, const string s, const double dx,
                 const int Nx, const double xmin, const double alpha,
                 const double lambda, const double omega, const double t)
{
	ofstream out(s.c_str());
  double x, xi, xil;
  double h1, h2, h3;
  cmplx ana;
	for(int i=0; i<Nx; i++){
		x = xmin + i * dx;
    xi = alpha * x;
    xil = alpha * lambda;
    h1 = -0.5 * pow(xi - xil*cos(omega*t),2 );
    h2 = omega*t/2 + xi * xil * sin(omega*t);
    h3 =  - 0.25 * xil*xil* sin(2*omega*t);
    ana = cmplx( h1 , h2 + h3  );
    ana = sqrt(alpha / sqrt(M_PI) ) * exp(ana);
		out << x << "\t" << norm(v[i]) << "\t" << v[i].real() << "\t" << v[i].imag()
         << "\t" << norm(ana) << "\t" << ana.real() << "\t" << ana.imag() <<  endl;
	}
	out.close();
}
//-----------------------------------

void init( cmplx* const psi0, const double alpha, const double lambda,
           const double dx, const double dt,
           const int Nx, const double xmin)
{
	const double x0 = dx*Nx * 0.5;
	for(int i=0;i<Nx; i++){
		double x = xmin + i*dx ;
		psi0[i] = sqrt(alpha/sqrt(M_PI)) * exp(- pow(alpha*(x-lambda),2)/2 );
	}
}
