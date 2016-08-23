#include "geometry.h"
#include <cstdlib>
#include "utils.h"
#include <math.h>
#include <complex>
namespace Cartesian{

  /**
   * A 1D initializer.
   * Simply allocate the arrays and instanciate the number of points and the spatial step.
   * Set dx1=spatial_dx1
   * Set n1=N_points_dx1 
   * Allocate the wavefunction (wave) as an array of size (n1+2)
   * Allocate the potential (v) as an array of size (n1+2)
   * Allocate the x1 (x1) as an array of (n1+2)
   * Write x[i] such that the 0 is not included. wf.x1[i]=(-(wf.n1+1)/2.+i)*wf.dx1;
   */

  void initialize1D_cartesian( wavefunction &wf, int N_points_x1, double spatial_x1)
  {

    wf.n1=N_points_x1;
    wf.dx1=spatial_x1;
    
    for(int i=0;i<(wf.n1+2);i++)
      {
	wf.wave.push_back(0.);
	wf.v.push_back(0.);
      }
    
    for(int i=0;i<(wf.n1+2);i++)
      wf.x1.push_back(0.);
    
    for(int i=1;i<=wf.n1;i++)
      wf.x1[i]=(-(wf.n1+1)/2.+i)*wf.dx1;
  }

  void softcore_potential_1D(wavefunction &wf, double Charge, double soft_parameter)
  {
    for(int i=1;i<=wf.n1;i++)
      wf.v[i]=Charge/sqrt(soft_parameter+wf.x1[i]*wf.x1[i]);
     
  }

  void guess_function_gauss(wavefunction &wf,double sigma, double x0,double initial_momentum)
  {
    complex<double>  I=complex<double> (0.,1);
    for(int i=1;i<=wf.n1;i++)
      {
	double arg=(wf.x1[i]-x0)/sigma;
	wf.wave[i]=exp(-arg*arg)*exp(-I*initial_momentum*wf.x1[i]);
      }
  }

  void normalize1D(wavefunction &wf)
  {
    double norm=0.;
    
    for(int i=1;i<=wf.n1;i++)
      norm+=real(wf.dx1*conj(wf.wave[i])*wf.wave[i]);
    
    //printf("\nNorm: In=%e",norm);
    for(int i=1;i<=wf.n1;i++)
      wf.wave[i]/=sqrt(norm); 
    
    norm=0.;  
    for(int i=1;i<=wf.n1;i++)
      norm+=real(wf.dx1*conj(wf.wave[i])*wf.wave[i]);
    
    //printf(" Out=%e",norm);
    
  }

  void ImaghamiltonianX1( double time_step, wavefunction &wf)
  {
  
    double lightc=137.;
    int N1=wf.n1;
    double dx1=wf.dx1;
    int k,j,i;
    complex<double>  ax1,cx1;
  
    double kinetic_x1=2.;
  
    vector<complex<double> > p1;
    vector<complex<double> > sl1;
    vector<complex<double> > vl1;
  
    for (i=0; i<(wf.n1+2); i++)
      {
	p1.push_back(0.);
	sl1.push_back(0.);
	vl1.push_back(0.);
      }
  
    ax1=complex<double> (-time_step/kinetic_x1/dx1/dx1,0.);
    cx1=complex<double> (-time_step/kinetic_x1/dx1/dx1,0.);
  
    for(i=1;i<=wf.n1;i++)
      {
	vl1[i]=complex<double> (1.+2.*time_step/dx1/dx1/kinetic_x1,0.)+complex<double> (wf.v[i],0.)*time_step;
	sl1[i]=
	  -ax1*wf.wave[i-1]
	  +(complex<double> (1.-2.*time_step/dx1/dx1/kinetic_x1,0.)-complex<double> (wf.v[i],0.)*time_step)*wf.wave[i]
	  -cx1*wf.wave[i+1];
      }
  
    tridag(ax1,vl1,cx1,sl1,p1,wf.n1);
  
    for(i=1;i<=wf.n1;i++)
      wf.wave[i]=p1[i];
  }

  
  void HamiltonianX1(double time_step, wavefunction &wf)
  {
  
    double lightc=137.;
    int N1=wf.n1;
    double dx1=wf.dx1;
    int k,j,i;
    complex<double>  ax1,cx1;
  
    double kinetic_x1=2.;
  
    vector<complex<double> > p1;
    vector<complex<double> > sl1;
    vector<complex<double> > vl1;
  
    for (i=0; i<(wf.n1+2); i++)
      {
	p1.push_back(0.);
	sl1.push_back(0.);
	vl1.push_back(0.);
      }
  
    ax1=complex<double> (0.,-time_step/kinetic_x1/dx1/dx1);
    cx1=complex<double> (0.,-time_step/kinetic_x1/dx1/dx1);
  
  
    for(i=1;i<=wf.n1;i++)
      {
	vl1[i]=complex<double> (1.,2.*time_step/dx1/dx1/kinetic_x1)+complex<double> (0.,wf.v[i])*time_step;
	sl1[i]=
	  -ax1*wf.wave[i-1]
	  +conj(vl1[i])*wf.wave[i]
	  -cx1*wf.wave[i+1];
      }
  
    tridag(ax1,vl1,cx1,sl1,p1,wf.n1);
  
    for(i=1;i<=wf.n1;i++)
      wf.wave[i]=p1[i];
  }

  double energy1D(wavefunction &wf)
  {
    double kinetic_x1=2.;
    complex<double>  energy=complex<double> (0.,0.);
    complex<double>  kinetic;
   
    for(int i=1;i<=wf.n1;i++)
      {
	kinetic=-(wf.wave[i-1]-2.*wf.wave[i]+wf.wave[i+1])/wf.dx1/wf.dx1/kinetic_x1;
	energy+=wf.dx1*conj(wf.wave[i])*(kinetic+wf.v[i]*wf.wave[i]);
      }

    return real(energy);
  }

  double norm1D(wavefunction &wf)
  {
    double norm=0.;
    for(int i=1;i<=wf.n1;i++)
      norm+=real(wf.dx1*conj(wf.wave[i])*wf.wave[i]);    
    return norm;
  }



  /***********************************************************************************************/
  //2D cartesian
  /***********************************************************************************************

  void initialize2D( wavefunction &wf, int Number_points_x1 ,double spatial_step_x1, int Number_points_x2 ,double spatial_step_x2)
  {
    
    wf.n1=Number_points_x1;
    wf.n2=Number_points_x2;

    for(int i=0;i<(Number_points_x1+2)*(Number_points_x2+2);i++)
      wf.wave.push_back(0.);
    
    for(int i=0;i<(Number_points_x1+2);i++)
      wf.x1.push_back(0.);

    for(int i=0;i<(Number_points_x2+2);i++)
      wf.x2.push_back(0.);
    
       
    wf.dx1=spatial_step_x1;
    wf.dx2=spatial_step_x2;
    
    for(int i=1;i<=(Number_points_x1);i++)
      wf.x1[i]=(-(Number_points_x1+1)/2.+i)*wf.dx1;

    for(int i=1;i<=(Number_points_x2);i++)
      wf.x2[i]=(-(Number_points_x2+1)/2.+i)*wf.dx2;
    
  }

  void normalize2D(wavefunction &wf)
  {
  //printf("\nNormalizing..");
    double norm=0.;//normalize(wf);
    
    for(int i=0;i<wf.wave.size();i++)
      norm+=real(wf.dx1*wf.dx2*conj(wf.wave[i])*wf.wave[i]);
    
    
    printf("\nIn=%e",norm);
    for(int i=0;i<=wf.wave.size();i++)
      wf.wave[i]/=sqrt(norm); 

    norm=0.;   
    for(int i=0;i<wf.wave.size();i++)
      norm+=real(wf.dx1*wf.dx2*conj(wf.wave[i])*wf.wave[i]);
    
    printf(" Out=%e",norm);
 }

  double norm2D(wavefunction &wf)
  { 

    double norm=0.;//normalize(wf);
    
    for(int i=0;i<wf.wave.size();i++)
      norm+=real(wf.dx1*wf.dx2*conj(wf.wave[i])*wf.wave[i]);

    return norm;

  }

 double energy2D(wavefunction &wf, vector<double> &v)
 {  

   double kinetic_x1=2.;
   double kinetic_x2=2.;
   complex<double>  energy=complex<double> (0.,0.);
   complex<double>  kin1,kin2;

   //int ee=wf.in2(10,10);
   for(int j=1;j<=wf.n2;j++)
     for(int i=1;i<=wf.n1;i++)
       {

	 kin1=-(wf.wave[wf.in2(j,i-1)]-2.*wf.wave[wf.in2(j,i)]+wf.wave[wf.in2(j,i+1)])/wf.dx1/wf.dx1/kinetic_x1;
	 kin2=-(wf.wave[wf.in2(j-1,i)]-2.*wf.wave[wf.in2(j,i)]+wf.wave[wf.in2(j+1,i)])/wf.dx2/wf.dx2/kinetic_x2;

	 energy+=wf.dx1*wf.dx2*conj(wf.wave[j*(wf.n1+2)+i])*(kin1+kin2+v[j*(wf.n1+2)+i]*wf.wave[j*(wf.n1+2)+i]);
       }
   
   return real(energy);
 }
  
  double obs_r_2D(wavefunction &wf)
  { 
    
    double obs_r=0.;//normalize(wf);
    
    for(int j=1;j<=wf.n2;j++)
      for(int i=1;i<wf.n1;i++)
      obs_r+=real(wf.dx1*wf.dx2*conj(wf.wave[wf.in2(j,i)])*
		  (sqrt(wf.x1[i]*wf.x1[i]+wf.x2[j]*wf.x2[j]))*
		  wf.wave[wf.in2(j,i)]);

    return obs_r;

  }


  void ImaghamiltonianX1_2D( double time_step, wavefunction &wf, vector<double> &v)
  {
    
    double lightc=137.;
   
    int k,j,i;
    complex<double>  ax1,cx1;
  
    double kinetic_x1=2.;
    double potential_split=2.;

    vector<complex<double> > p1;
    vector<complex<double> > sl1;
    vector<complex<double> > vl1;
  
    for (i=0; i<(wf.n1+2); i++)
      {
	p1.push_back(0.);
	sl1.push_back(0.);
	vl1.push_back(0.);
      }
  
    ax1=complex<double> (-time_step/kinetic_x1/wf.dx1/wf.dx1,0.);
    cx1=complex<double> (-time_step/kinetic_x1/wf.dx1/wf.dx1,0.);
  
    for(j=1;j<=wf.n2;j++)
      {
	for(i=1;i<=wf.n1;i++)
	  {
	    vl1[i]=complex<double> (1.+2.*time_step/wf.dx1/wf.dx1/kinetic_x1,0.)+complex<double> (v[j*(wf.n1+2)+i],0.)*time_step/potential_split;
	    sl1[i]=
	      -ax1*wf.wave[wf.in2(j,i-1)]
	      +(complex<double> (1.-2.*time_step/wf.dx1/wf.dx1/kinetic_x1,0.)-complex<double> (v[j*(wf.n1+2)+i],0.)*time_step/potential_split)*wf.wave[wf.in2(j,i)]
	      -cx1*wf.wave[wf.in2(j,i+1)];
	  }
	
	
	tridag(ax1,vl1,cx1,sl1,p1,wf.n1);

	
	for(i=1;i<=wf.n1;i++)
	  wf.wave[wf.in2(j,i)]=p1[i];
      }  
  }

   void ImaghamiltonianX2_2D( double time_step, wavefunction &wf, vector<double> &v)
  {
    
    double lightc=137.;
   
    int k,j,i;
    complex<double>  ax2,cx2;
  
    double kinetic_x2=2.;
    double potential_split=2.;

    vector<complex<double> > p2;
    vector<complex<double> > sl2;
    vector<complex<double> > vl2;
  
    for (i=0; i<(wf.n2+2); i++)
      {
	p2.push_back(0.);
	sl2.push_back(0.);
	vl2.push_back(0.);
      }
  
    ax2=complex<double> (-time_step/kinetic_x2/wf.dx2/wf.dx2,0.);
    cx2=complex<double> (-time_step/kinetic_x2/wf.dx2/wf.dx2,0.);
  
    for(i=1;i<=wf.n1;i++)
      {
	for(j=1;j<=wf.n2;j++)
	  {
	    vl2[j]=complex<double> (1.+2.*time_step/wf.dx2/wf.dx2/kinetic_x2,0.)+complex<double> (v[j*(wf.n1+2)+i],0.)*time_step/potential_split;
	    sl2[j]=
	      -ax2*wf.wave[wf.in2(j-1,i)]
	      +(complex<double> (1.-2.*time_step/wf.dx2/wf.dx2/kinetic_x2,0.)-complex<double> (v[j*(wf.n1+2)+i],0.)*time_step/potential_split)*wf.wave[wf.in2(j,i)]
	      -cx2*wf.wave[wf.in2(j+1,i)];
	  }
	
	tridag(ax2,vl2,cx2,sl2,p2,wf.n2);
	
	for(j=1;j<=wf.n2;j++)
	  wf.wave[wf.in2(j,i)]=p2[j];
      }  
  }

  void HamiltonianX1_2D( double time_step, wavefunction &wf, vector<double> &v)
  {
    
    double lightc=137.;
   
    int k,j,i;
    complex<double>  ax1,cx1;
  
    double kinetic_x1=2.;
    double potential_split=2.;

    vector<complex<double> > p1;
    vector<complex<double> > sl1;
    vector<complex<double> > vl1;
  
    for (i=0; i<(wf.n1+2); i++)
      {
	p1.push_back(0.);
	sl1.push_back(0.);
	vl1.push_back(0.);
      }
  
    ax1=complex<double> (0.,-time_step/kinetic_x1/wf.dx1/wf.dx1);
    cx1=complex<double> (0.,-time_step/kinetic_x1/wf.dx1/wf.dx1);
  
    for(j=1;j<=wf.n2;j++)
      {
	for(i=1;i<=wf.n1;i++)
	  {
	    vl1[i]=complex<double> (1.,+2.*time_step/wf.dx1/wf.dx1/kinetic_x1)+complex<double> (0.,v[j*(wf.n1+2)+i])*time_step/potential_split;
	    sl1[i]=
	      -ax1*wf.wave[wf.in2(j,i-1)]
	      +conj(vl1[i])*wf.wave[wf.in2(j,i)]
	      -cx1*wf.wave[wf.in2(j,i+1)];
	  }
		
	tridag(ax1,vl1,cx1,sl1,p1,wf.n1);
	
	for(i=1;i<=wf.n1;i++)
	  wf.wave[wf.in2(j,i)]=p1[i];
      }
      
  }

  
  void HamiltonianX2_2D( double time_step, wavefunction &wf, vector<double> &v)
  {
    
    double lightc=137.;
   
    int k,j,i;
    complex<double>  ax2,cx2;
  
    double kinetic_x2=2.;
    double potential_split=2.;

    vector<complex<double> > p2;
    vector<complex<double> > sl2;
    vector<complex<double> > vl2;
  
    for (i=0; i<(wf.n2+2); i++)
      {
	p2.push_back(0.);
	sl2.push_back(0.);
	vl2.push_back(0.);
      }
  
    ax2=complex<double> (0.,-time_step/kinetic_x2/wf.dx2/wf.dx2);
    cx2=complex<double> (0.,-time_step/kinetic_x2/wf.dx2/wf.dx2);
  
    for(i=1;i<=wf.n1;i++)
      {
	for(j=1;j<=wf.n2;j++)
	  {
	    vl2[j]=complex<double> (1.,+2.*time_step/wf.dx2/wf.dx2/kinetic_x2)+complex<double> (0.,v[j*(wf.n1+2)+i])*time_step/potential_split;
	    sl2[j]=
	      -ax2*wf.wave[wf.in2(j-1,i)]
	      +conj(vl2[j])*wf.wave[wf.in2(j,i)]
	      -cx2*wf.wave[wf.in2(j+1,i)];
	  }
	
	
	tridag(ax2,vl2,cx2,sl2,p2,wf.n2);

	
	for(j=1;j<=wf.n2;j++)
	  wf.wave[wf.in2(j,i)]=p2[j];
      }  
  }

  */

}

namespace Cylindrical{


}





