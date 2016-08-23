// $Id: geometry.h 512 2006-11-21 15:22:16Z goerke $
#include <cstdlib>
#ifndef GEOMETRY_H
#define GEOMETRY_H
#include "constant.h"
#include "wavefunction.h"
#include "Projection.h"
#include "potentials.h"

using namespace std; // to save the std:: in front of vector and cout

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

namespace Cartesian_3D{

  
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

  void X1( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X2( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X3( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );

  /*
   * Propagation of one time step
   */
  void operator()( const complex<double>  time_step, wavefunction &wf, const ABVparam &p)
  {
    X1( time_step*.25 , wf , p );
    X2( time_step*.25 , wf , p );
    X3( time_step*.5  , wf , p );
    X2( time_step*.25 , wf , p );
    X1( time_step*.25 , wf , p );
  }

  void X1_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X2_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X3_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );

  void ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X1_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X2_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS_E_R( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double electric_field ); 
  void X3_ECS_A_P( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double electric_field ); 

}; // end of class Hamiltonian
  
  /**
   * Calculation of Observables:
   */

  double Obs_Norm( wavefunction &wf );
  double Obs_Energy( wavefunction &wf, ABVparam &p );
  double Obs_Expectation_Value_X1( wavefunction &wf );
  double Obs_Expectation_Value_X2( wavefunction &wf );
  double Obs_Expectation_Value_X3( wavefunction &wf );
  double Obs_Expectation_Value_Width_X1( wavefunction &wf );
  double Obs_Expectation_Value_Width_X2( wavefunction &wf );
  double Obs_Expectation_Value_Width_X3( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X1( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X2( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X3( wavefunction &wf );
  double Obs_Expectation_Value_r( wavefunction &wf );
  double Obs_Expectation_Value_Inv_r( wavefunction &wf );
  vector<double>  H2p_nuclei( Potential_H2_plus &pot );

  Projection& Prj_Wavefunction_X1( wavefunction &wf , Projection& prj, string title = "PrjX1" );
  Projection& Prj_Wavefunction_X2( wavefunction &wf , Projection& prj, string title = "PrjX2" );
  Projection& Prj_Wavefunction_X3( wavefunction &wf , Projection& prj, string title = "PrjX3" ); 
  //Projection& Prj_Wavefunction_X1_Low( wavefunction &wf, Projection& prj, string title="PrjX1_Low");
  Projection& Prj_Wavefunction_X1_High( wavefunction &wf , Projection& prj , string title = "PrjX1_High" );
  Projection& Prj_Wavefunction_X1_Low( wavefunction &wf , Projection& prj , double limit_atom_size , string title = "PrjX1_Low" );
  vector<complex<double> > Obs_Serialize( wavefunction &wf );
  void Mask_Function( wavefunction &wf , double frac_n1_right , double frac_n2_right , double frac_n3_right , double frac_n2_left , double frac_n3_left , double exponent);
  //vector<double> Ionization_HE( wavefunction &wf , double small_bound , double big_bound );
  vector<double> Ionization_He( wavefunction &wf );
  vector<double> Ionization_H2( wavefunction &wf , vector<double> &parameters );
  void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n1_left, double frac_n2_left, double frac_n3_left , double exponent);
  //vector<double> Ionization_HE( wavefunction &wf, double small_bound, double big_bound);
  vector<double> Ionization_He( wavefunction &wf);
  //double Obs_Expectation_Value_relative_coordinate( wavefunction &wf ); //Works only for the Helium and for the H2
  //double Obs_Expectation_InverseValue_relative_coordinate( wavefunction &wf ); //Works only for the Helium and for the H2
  //vector<double> Obs_HighHarmonicGeneration( wavefunction &wf ); 

  /**
   * Calculation of Derivatives:
   */
  complex<double>  First_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  First_Derivative_X2( wavefunction &wf , int k , int j , int i );
  complex<double>  First_Derivative_X3( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X2( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X3( wavefunction &wf , int k , int j , int i );

/*------------------  Wavefunction related functions  -----------------------*/

  // Attention !!!!!
     
  
  /**
   * A 3D initializer (void).  
   * Simply allocate the arrays and instanciate the number of points and the spatial step.   
   * Allocate the wavefunction (wave) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the potential (v) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the x1 (x1) as an array of (n1+2).   
   */

  void Initialize( wavefunction &wf, vector<int> N_points , vector<double> spatial_steps, int symmetry_x1 = 0, int symmetry_x2 = 0, int symmetry_x3 = 0);
  vector<double> H2p_Force(wavefunction &wf , Potential_H2_plus &pot, const double Ex, const double Ey);
  void Euler( wavefunction &wf, Potential_H2_plus &pot, double dt, const double Ex, const double Ey );
  void Runge_Kutta( wavefunction &wf, Potential_H2_plus &pot, double dt, const double Ex, const double Ey );

  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta );
  void Normalize( wavefunction &wf );

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf ); 
  complex<double>  Project( wavefunction &wf , wavefunction &wf2 );
  double Overlap( wavefunction &wf , wavefunction &wf2);
  complex<double>  Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf );

/*----------------  Classical Nuclear Motion functions  ---------------------*/

  double Calc_Distance( double velocity );
  double Calc_Velocity( const wavefunction &wf , double distance , Potential_H2_planar_parallel_deriv &pot );
  void Newton_Equation_Moving_Nuclei( const wavefunction &wf , vector<double> &initial_values , double time_step , Potential_H2_planar_parallel_deriv &pot );
  
  void Initial_Eng( wavefunction &wf );
  
 }

/**
 *  The namespace Cylindrical deals only Hamiltonians H=Dxx+r*Dx+Drr+ V.
 */

namespace Cylindrical_3D{

  
  // Attention !!!!!
     
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

  void X1( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X2( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X3( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );

  /*
   * Propagation of one time step
   */
  void operator()( const complex<double>  time_step, wavefunction &wf, const ABVparam &p)
  {
    X1( time_step*.25 , wf , p );
    X2( time_step*.25 , wf , p );
    X3( time_step*.5  , wf , p );
    X2( time_step*.25 , wf , p );
    X1( time_step*.25 , wf , p );
  }

  void X1_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X2_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );
  void X3_Laser( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double field, const gauge_t gauge );

  void ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X1_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X2_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X3_ECS_E_R( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double electric_field ); 
  void X3_ECS_A_P( const complex<double>  time_step, wavefunction &wf, const ABVparam &p,
      const double electric_field ); 

  /**
   * Adaptative mesh
   */

  void X2_AM( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X1_AM( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );
  void X3_AM( const complex<double>  time_step, wavefunction &wf, const ABVparam &p );  
  void X3_Laser_AM( const complex<double>  time_step, wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge );

}; // end of class Hamiltonian

  
  /**
   * One time step in Imaginary time, the default value for a in a*(d/dx^2) is a=1/2
   * IMPORTANT, the time step used comes from the Cranck-Nichols equation
   * The auxiliary vector are allocated and initialized inside the function, it uses the tridag (www.nr.com) function implemented in utils.h 
   */  

  /**
   * Calculation of Observables:
   */

  double Obs_Norm( wavefunction &wf );
  double Obs_Energy( wavefunction &wf, ABVparam &p );
  double Obs_Expectation_Value_X1( wavefunction &wf );
  double Obs_Expectation_Value_X2( wavefunction &wf );
  double Obs_Expectation_Value_X3( wavefunction &wf );
  double Obs_Expectation_Value_Width_X1( wavefunction &wf );
  double Obs_Expectation_Value_Width_X2( wavefunction &wf );
  double Obs_Expectation_Value_Width_X3( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X1( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X2( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X3( wavefunction &wf );
  double Obs_Expectation_Value_r( wavefunction &wf );
  double Obs_Expectation_Value_Inv_r( wavefunction &wf );
  vector<double> Obs_Expectation_Value_e_Nuclei_H2p( wavefunction &wf );
  vector<double> ground_dissociation_ionization_H2p( wavefunction &wf );


  Projection& Prj_Wavefunction_X1(wavefunction &wf, Projection& prj, string title="PrjX1");
  Projection& Prj_Wavefunction_X2(wavefunction &wf, Projection& prj, string title="PrjX2");
  Projection& Prj_Wavefunction_X3(wavefunction &wf, Projection& prj, string title="PrjX3"); 
  //Projection& Prj_Wavefunction_X1_Low( wavefunction &wf, Projection& prj, string title="PrjX1_Low");
  Projection& Prj_Wavefunction_X1_High(wavefunction &wf, Projection& prj, string title="PrjX1_High");
  Projection& Prj_Wavefunction_X1_Low(wavefunction &wf, Projection& prj,  double limit_atom_size, string title="PrjX1_Low");
  vector<complex<double> > Obs_Serialize( wavefunction &wf );
  void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n2_left, double frac_n3_left , double exponent );
  //vector<double> Ionization_HE( wavefunction &wf, double small_bound, double big_bound);
  vector<double> Ionization_He( wavefunction &wf);
  vector<double> Ionization_H2( wavefunction &wf , vector<double> &parameters );
  //double Obs_Expectation_Value_relative_coordinate( wavefunction &wf ); //Works only for the Helium and for the H2
  //double Obs_Expectation_InverseValue_relative_coordinate( wavefunction &wf ); //Works only for the Helium and for the H2
  //vector<double> Obs_HighHarmonicGeneration(wavefunction &wf); 

  /**
   * Calculation of Derivatives:
   */
  complex<double>  First_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  First_Derivative_X2( wavefunction &wf , int k , int j , int i );
  complex<double>  First_Derivative_X3( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X2( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X3( wavefunction &wf , int k , int j , int i );

/*------------------  Wavefunction related functions  -----------------------*/

  /**
   * A 3D initializer (void).  
   * Simply allocate the arrays and instanciate the number of points and the spatial step.   
   * Allocate the wavefunction (wave) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the potential (v) as an array of size (n1+2)*(n2+2)*(n3+2).   
   * Allocate the x1 (x1) as an array of (n1+2).   
   */

  void Initialize( wavefunction &wf , vector<int> N_points , vector<double> spatial_steps, int symmetry_x2=0 );
    
  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta );
  void Normalize( wavefunction &wf );

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf); 
  complex<double>  Project( wavefunction &wf , wavefunction &wf2);
  double Overlap( wavefunction &wf , wavefunction &wf2);
  complex<double>  Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf);
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf);
  /**
   * Adaptative mesh
   */

  void SizeGrid_X1(wavefunction &wf, double criteria, int grow);
  void SizeGrid_X2(wavefunction &wf, double criteria, int grow);
  void SizeGrid_X3(wavefunction &wf, double criteria, int grow);
  void SizeGrid(wavefunction &wf, double criteria, int grow);

/*----------------  Classical Nuclear Motion functions  ---------------------*/

  void Newton_Equation_Moving_Nuclei( const wavefunction &wf , vector<double> &initial_values , double time_step , Potential_H2 &pot );
  double Calc_Distance( double velocity );
  double Calc_Velocity( const wavefunction &wf , double distance , Potential_H2_deriv &pot );

  void Initial_Eng( wavefunction &wf );

}




namespace Cylindrical_2D{

  
  // Attention !!!!!
     
  void Initialize( wavefunction &wf , vector<int> N_points , vector<double> spatial_steps );
  void Initialize_Hamiltonian(  wavefunction &wf, const double A_x1, const double A_x2, const double B_x1, const double B_x2, const string potential , const vector<double> Potential_Parameters );
  // void Initialize_Hamiltonian( const wavefunction &wf, const Potential_H2_plus &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_H2_plus_2D &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_H2 &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_H2_plus_with_two_electrons &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_He &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_He_plus &pot );
  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_H_2D &pot );

  //void Initialize_Hamiltonian( const wavefunction &wf, const Potential_Free &pot );

  //void Create_Potential( const wavefunction &wf , const Potential_H2_plus &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_H2_plus_2D &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_H2 &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_H2_deriv &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_H2_plus_with_two_electrons &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_He &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_He_plus &pot );
  //void Create_Potential( const wavefunction &wf , const Potential_H_2D &pot );

  //void Create_Potential( const wavefunction &wf , const Potential_Free &pot );

  void Create_Potential( const wavefunction &wf , const string potential , const vector<double> SAEparam );
  void Create_Potential_deriv( const wavefunction &wf , const string potential , const vector<double> /*Potential_Parameters*/ );
    
  void Guess_Function_Gauss( wavefunction &wf , vector<double> sigmas , vector<double> gauss_centers , vector<double> initial_momenta );
  void Normalize( wavefunction &wf );

  void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf); 
  complex<double>  Project( wavefunction &wf , wavefunction &wf2);
  double Overlap( wavefunction &wf , wavefunction &wf2);
  complex<double>  Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf);
  void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf);
  /**
   * One time step in Imaginary time, the default value for a in a*(d/dx^2) is a=1/2
   * IMPORTANT, the time step used comes from the Cranck-Nichols equation
   * The auxiliary vector are allocated and initialized inside the function, it uses the tridag (www.nr.com) function implemented in utils.h 
   */  

  void Hamiltonian_X1( const complex<double>  time_step , wavefunction &wf );  
  void Hamiltonian_X2( const complex<double>  time_step , wavefunction &wf );  


    /**
     * Propagation of one time step
     */
  inline void Hamiltonian( const complex<double>  time_step , wavefunction &wf )
  {
    Hamiltonian_X1( time_step*.5 , wf );
    Hamiltonian_X2( time_step , wf );
    Hamiltonian_X1( time_step*.5 , wf );
    
  } // end Of Hamiltonian

  void Hamiltonian_X2_Laser( const complex<double>  time_step , wavefunction &wf , const double field , const gauge_t gauge );


  void Hamiltonian_ECS( const complex<double>  time_step , wavefunction &wf );
  void Hamiltonian_X1_ECS( complex<double>  time_step , wavefunction &wf );  
  void Hamiltonian_X2_ECS( complex<double>  time_step , wavefunction &wf );  


  void Initial_Eng( wavefunction &wf );

  /**
   * Calculation of Observables:
   */

  double Obs_Norm( wavefunction &wf );
  double Obs_Energy( wavefunction &wf );
  double Obs_Expectation_Value_X1( wavefunction &wf );
  double Obs_Expectation_Value_X2( wavefunction &wf );

  double Obs_Expectation_Value_Width_X1( wavefunction &wf );
  double Obs_Expectation_Value_Width_X2( wavefunction &wf );

  double Obs_Expectation_Value_Inv_X1( wavefunction &wf );
  double Obs_Expectation_Value_Inv_X2( wavefunction &wf );

  double Obs_Expectation_Value_r( wavefunction &wf );
  double Obs_Expectation_Value_Inv_r( wavefunction &wf );
  vector<double> Obs_Expectation_Value_e_Nuclei_H2p( wavefunction &wf );
  vector<double> ground_dissociation_ionization_H2p( wavefunction &wf );


  Projection& Prj_Wavefunction_X3(wavefunction &wf, Projection& prj, string title="PrjX3"); 
  //Projection& Prj_Wavefunction_X1_Low( wavefunction &wf, Projection& prj, string title="PrjX1_Low");
  Projection& Prj_Wavefunction_X1_High(wavefunction &wf, Projection& prj, string title="PrjX1_High");
  Projection& Prj_Wavefunction_X1_Low(wavefunction &wf, Projection& prj,  double limit_atom_size, string title="PrjX1_Low");
  vector<complex<double> > Obs_Serialize( wavefunction &wf );
  void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n2_right, double frac_n3_right, double frac_n2_left, double frac_n3_left , double exponent );
  //vector<double> Ionization_HE( wavefunction &wf, double small_bound, double big_bound);

  /**
   * Calculation of Derivatives:
   */
  complex<double>  First_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  First_Derivative_X2( wavefunction &wf , int k , int j , int i );

  complex<double>  Second_Derivative_X1( wavefunction &wf , int k , int j , int i );
  complex<double>  Second_Derivative_X2( wavefunction &wf , int k , int j , int i );


  /**
   * Adaptative mesh
   */

  void SizeGrid_X1(wavefunction &wf, double criteria, int grow);
  void SizeGrid_X2(wavefunction &wf, double criteria, int grow);

  void SizeGrid(wavefunction &wf, double criteria, int grow);
  void Hamiltonian_X2_AM(const complex<double>  time_step , wavefunction &wf );
  void Hamiltonian_X1_AM( const complex<double>  time_step , wavefunction &wf );
}










#endif	/* GEOMETRY_H */
