//created by Jing on 02/16/2010
#include <cstdlib>
#ifndef CARTESIAN1D_H
#define CARTESIAN1D_H
#include "constant.h"
#include "wavefunction.h"
#include "Projection.h"
#include "potentials.h"

using namespace std; // to save the std:: in front of vector and cout

namespace Cartesian_1D{

/**
 * Declare everything what is needed for the Crank-Nicholson-Scheme.  
 * independent of geometry namespace
 */

struct temp_data_t {
/**
 * temporary storage for Tridag
 */
  vector<complex<double> > gam;

/**
 * the 3 diagonals
 */
  vector<complex<double> > tridag_upp;
  vector<complex<double> > tridag_mid;
  vector<complex<double> > tridag_low;

/**
 * the potential with respect to the coordinate
 */
  vector<complex<double> > v_1D;
  vector<complex<double> > wf_1D;
  vector<complex<double> > wf_1D_rightside;
  vector<complex<double> > wf_1D_solution;
  
};
  
/**
 *  The namespace Cartesian deals only Hamiltonians H=Dxx + V.
 */

  
/*--------------------------  Hamiltonian class  ----------------------------*/

class Hamiltonian {
  public:

  /**
   * The crucial grid parameters
   * (introduced here to check if wf dimensions match grid dimensions)
   */
  double dx[DIM+1];/**< Spatial Step x[j]( i ) */  
  int n[DIM+1];/**< Number of points in the N[j] dimension */  

  temp_data_t temp_data[DIM+1]; /* provide workspace for DIM dimensions */

  /**
   * Constructor just Initializes the workspace
   */
  Hamiltonian( const wavefunction &wf );
  
  /**
   * Propagation Operators
   */

  void X( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );


  /*
   * Propagation of one time step
   */
  void operator()( const complex<double>  time_step, wavefunction &wf, const ABVparam &p)
	  {
		  X( time_step*.5 , wf , p );

	  }

  void X_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X_ECS_E_R(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage );
  void X_ECS_E_R_back_left(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage );
   void X_ECS_E_R_back_right(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage );
  void X_ECS_E_R_unequalabsorb(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage_left, double absorb_percentage_right );

  //void ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  //void X1_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  //void X2_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  //void X3_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  //void X3_ECS_E_R( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
  //   const double electric_field ); 
  //void X3_ECS_A_P( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
  //    const double electric_field ); 


 /**
   * Adaptative mesh
   */

  //void X1_AM( const complex<double>  time_step , wavefunction &wf , const ABVparam &p );
  //void X2_AM(const complex<double>  time_step , wavefunction &wf  , const ABVparam &p );
  
  
  //void X2_Laser_AM( const complex<double>  time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge);
  //void X3_Laser_AM( const complex<double>  time_step , wavefunction &wf , const ABVparam &p , const double field , const gauge_t gauge );

}; // end of class Hamiltonian
  
  /**
   * Calculation of Observables:
   */

  double Obs_Norm( wavefunction &wf );
  double Obs_Energy( wavefunction &wf, ABVparam &p );
  double Obs_Expectation_Value_X( wavefunction &wf );

 
  double Obs_Expectation_Value_Width_X( wavefunction &wf );


  complex<double>  First_Derivative_X( wavefunction &wf , int j , int i );
 
  complex<double>  Second_Derivative_X( wavefunction &wf , int j , int i );
 

/*------------------  Stuff needed for FFT and Masks  -----------------------*/


  void Initialize_Momentum( wavefunction &wf, wavefunction &wf_mom);
  void FFT( wavefunction &wf , wavefunction &wf_transformed);
  void FFT1D( wavefunction &wf , wavefunction &wf_transformed);


  void Initialize( wavefunction &wf, int n_point , double spatial_step, int symmetry_x1 = 0);
  //void Runge_Kutta( wavefunction &wf, Potential_H2_plus_e3D_n2D &pot, double dt, const double Ex, const double Ey );
  //vector<double> H2p_Force(wavefunction &wf , Potential_H2_plus_e3D_n2D &pot, const double Ex, const double Ey);

  void Guess_Function_Gauss( wavefunction &wf , double sigma , double gauss_center , double initial_momentum );
  void Guess_Function_Uniform( wavefunction &wf , double wavefunction_value );
  void Normalize( wavefunction &wf );

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf ); 
  complex<double>  Project( wavefunction &wf , wavefunction &wf2 );
  double Overlap( wavefunction &wf , wavefunction &wf2);
  complex<double> Dipole( wavefunction &wf, wavefunction &wf2 );
  complex<double> Dipole_k( wavefunction &wf, double &k );
  void Mask_Function( wavefunction &wf , double frac_n1_right, double frac_n1_left, double exponent);
  void Mask_Function_Inner( wavefunction &wf , double x_left, double x_right, double exponent);
  void Mask_Function_Part( wavefunction &wf , double x_left, double x_right, double exponent);
  complex<double>  Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );
  Projection1D& Prj_Wavefunction_X(wavefunction &wf, Projection1D& prj,  string title);

  
 }

#endif
