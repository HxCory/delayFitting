//Cartesian1D.cpp created by Jing on 02/16/2010
#include <cstdlib>
#include "Cartesian1D.h"
#include "utils.h"
#include "potentials.h"
#include <math.h>
#include <complex>
#include <fftw3.h>
#include "Projection.h" 
#include "ParameterMap.h" 
#include "wavefunction.h"
#include "constant.h"

namespace Cartesian_1D{

/*--------------------------  Hamiltonian class  ----------------------------*/

Hamiltonian::Hamiltonian(const wavefunction &wf) {
  /**
   *  Note the grid dimensions that were used to initialize the
   * ABV and temp data vectors
   */
  for (int coord_index=1; coord_index<2; ++coord_index) {
    n[coord_index]  = wf.n[coord_index];
    dx[coord_index] = wf.x[coord_index].delta;
  }

  /**
   * Allocate the temporary storage for Tridag
   */
  for ( int coord_index=1; coord_index <2; ++coord_index ) {
    int num = wf.n[coord_index] + 2;
    temp_data_t &d=temp_data[coord_index];

    d.gam.resize(num, complex<double> (0., 0.) );
    d.tridag_low.resize(num, complex<double> (0., 0.) );
    d.tridag_mid.resize(num, complex<double> (0., 0.) );
    d.tridag_upp.resize(num, complex<double> (0., 0.) );
    d.v_1D.resize(num, complex<double> (0., 0.) );
    d.wf_1D.resize(num, complex<double> (0., 0.) );
    d.wf_1D_rightside.resize(num, complex<double> (0., 0.) );
    d.wf_1D_solution.resize(num, complex<double> (0., 0.) );
  }
};

  /**
   * One time step in Real time (NO LASER), the default value for a in a*(d/dx^2) is a=1/2
   * IMPORTANT, the time step used comes from the Cranck-Nichols equation
   * The auxiliary vector are allocated and initialized inside the function, it uses the tridag (www.nr.com) function implemented in utils.h 
   */  

    /**
     * Propagation of one time step in X1 direction
     * Look at doc/abv.tex Eq. (6) in doc/ABV/abv.tex for details.
     */

  void Hamiltonian::X( const complex<double>  time_step , wavefunction &wf , const ABVparam &p )
  {
    temp_data_t &d=temp_data[1]; 
    if(n[1] != wf.n1 || dx[1] != wf.dx1) {
      cerr<<"wavefunction has wrong size for Hamiltonian\n";
      exit(1);
    }
    
    const complex<double>  idt=I*time_step;
    const complex<double>  arg_A = ( idt*.5*wf.one_by_dx1sqr )*p.ABV_A_x;
    const complex<double>  arg_B = 0.;//complex<double> (0.,0.);//( idt*.25*wf.one_by_dx1   )*p.ABV_B_x1;
    complex<double>  arg_V;
    complex<double>  tridag_low_Fast =arg_A-arg_B;
    complex<double>  tridag_upp_Fast =arg_A+arg_B;


	    for( int i = 0 ; i < wf.n1+2 ; ++i )//Copy the wavefunction
	    {
		    /**
		     * Creates a 1D potential and a 1D wavefunction out of the 3D ones
		     */
		    const int index=wf.in2(1,i);
		    d.v_1D[ i ] = p.ABV_V[ index ];
		    d.wf_1D[ i ] = wf.wave[  index ];
	      }
	    
	    
	    for( int i = 1 ; i <= wf.n1 ; ++i )
	    {
		    /**
		     * Definition of the arguments in Eq. (7-9) in doc/ABV/abv.tex
		     */
		    
		    arg_V = ( idt*.5 )*d.v_1D[ i ];
		    
		    /**
		     * Define the 3 diagonals X of the left side of Eq. (6)  in doc/ABV/abv.tex => look at Eqs. (7-9) in doc/ABV/abv.tex
		     */
		    d.tridag_mid[ i ] =1.-2.*arg_A+arg_V;
		    
		    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ]+( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//general case, symmetry==1
		    
		    /**
		     * for (anti-)symmetric wavefunctions one has to take care of the right boundary conditions (phi(0)=+/-phi(1). That leads to modified values.
		     */
		    if (i==1) {
			    if (wf.symmetry_x1==1)
			    {
				    d.tridag_mid[ i ] =1.-arg_A-arg_B+arg_V;
				    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
			    }
			    
			    if (wf.symmetry_x1==-1)
			    {
				    d.tridag_mid[ i ] =1.-3.*arg_A+arg_B+arg_V;
				    d.wf_1D_rightside[ i ]=-tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ]-tridag_upp_Fast*d.wf_1D[ i+1 ];
			    }
		    }
		    
	    }
	    
	    /**
	     * Solve the system of linear equation of Eq. (6) in doc/ABV/abv.tex
	     */
	    Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
	    
	    /**
	     * Copy back the 1D wavefunction to the 3D one
	     */
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
		    wf.wave[ wf.in2(1,i) ] = d.wf_1D_solution[ i ];
	    }

      
    
  } // end of X1_3D

  
  void Hamiltonian::X_Laser(const complex<double>  time_step , wavefunction &wf, const ABVparam &p, const double field, const gauge_t gauge )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
	  
	  const complex<double>  idt=I*time_step;
	  const complex<double>  arg_A = idt*.5*wf.one_by_dx1sqr *p.ABV_A_x;
	  const complex<double>  f_arg_B=time_step*.25*wf.one_by_dx1; 
	  complex<double>  arg_B;
	  complex<double>  arg_V;
	  complex<double>  tridag_upp_Fast;
	  complex<double>  tridag_low_Fast;
	  

		  for( int i = 0 ; i < wf.n1+2 ; i++ )
		  {
			  const int index=wf.in2( 1,i );
			  d.v_1D[ i ]  = p.ABV_V[ index ];
			  d.wf_1D[ i ] = wf.wave[ index ];
		  }
		  
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  
			  if ( gauge == lengthgauge )
			  {
				  arg_B = 0;//f_arg_B       *p.ABV_B_x2;
				  arg_V = idt*0.5*( /*p.ABV_B_x/p.ABV_A_x*0.5**/wf.x1[ i ]*field + d.v_1D[ i ] ); //p.ABV_A_x2*0.5 is coming from /doc/ABV/factors
			  }
			  else if ( gauge == velocitygauge )
			  {
				  arg_B = f_arg_B*( -p.ABV_B_x*field*one_by_lightC_au );
				  arg_V = idt*.5                  *d.v_1D[ i ];
			  }
			  
			  tridag_low_Fast    = arg_A-arg_B;	     
			  d.tridag_mid[ i ] = 1.-2.*arg_A+arg_V;
			  tridag_upp_Fast    = arg_A+arg_B;	     
			  
			  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+2.*arg_A-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];//sym==1
			  
			  if (i==1) 
			  {
				  if (wf.symmetry_x1==1)
				  {
					  d.tridag_mid[ i ] = 1.-arg_A-arg_B+arg_V;
					  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+arg_A+arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
				  }
				  
				  if (wf.symmetry_x1==-1)
				  {
					  d.tridag_mid[ i ] = 1.-3.*arg_A+arg_B+arg_V;
					  d.wf_1D_rightside[ i ]= -tridag_low_Fast*d.wf_1D[ i-1 ] +( 1.+3.*arg_A-arg_B-arg_V )*d.wf_1D[ i ] -tridag_upp_Fast*d.wf_1D[ i+1 ];
				  }
			  }
			  
		  }
		  
		  Tridag_Fast( tridag_low_Fast , d.tridag_mid, tridag_upp_Fast , d.wf_1D_rightside, d.wf_1D_solution, d.gam );
		  
		  
		  for( int i = 1 ; i <= wf.n1 ; i++ )
		  {
			  wf.wave[ wf.in2( 1 , i ) ] = d.wf_1D_solution[ i ];
		  }
		  

  } // end of X1_2D_Laser

  void Hamiltonian::X_ECS_E_R(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
    
    int N_z_right=int(wf.n1*absorb_percentage);
    int N_z_left=int(wf.n1*absorb_percentage);
	
	for( int j = 0 ; j < wf.n1+2 ; j++ )
	  {
	    d.v_1D[ j ]  = p.ABV_V[ wf.in2( 1 , j ) ];
	    d.wf_1D[ j ] = wf.wave[  wf.in2(  1 , j ) ];
	  }
	
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  {
	    complex<double>  arg_A = I*time_step/( 2.*wf.dx1*wf.dx1 )*p.ABV_A_x;
	    complex<double>  arg_B = 0.;//( I*time_step/( 4.*wf.dx2 )    )*ABV_B_x2;
	    //complex<double>  arg_V = ( I*time_step/2.               )*d.v_1D[ j ];
	    complex<double>  arg_V =   I*time_step/2.               *(d.v_1D[ j ]+wf.x1[j]*electric_field);
	    

	    if (wf.symmetry_x1==0) // X1 starts from -X2_max to X2_max
	      {   // ECS is done in the left and right side of X2
		
		if (j>N_z_left && j<wf.n1-N_z_right) // center of X2
		  {
		    d.tridag_low[ j ]     = arg_A-arg_B;	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A+arg_V;
		    d.tridag_upp[ j ]     = arg_A+arg_B;	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		
		if (j<N_z_left )  // left of x2
		  {
		    d.tridag_low[ j ]     = arg_A*exp(-2.*I*eta);	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*eta)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*exp(-2.*I*eta);	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-2.*I*eta)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		if (j==N_z_left) // left discontinuity point
		  {
		    d.tridag_low[ j ]     = arg_A*(2.*exp(-I*eta))/(1.+exp(I*eta));  
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*eta)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*2./(1.+exp(I*eta));	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-I*eta)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
	      } // end if
	    
	
	    
	    if ( j>wf.n1-N_z_right) // right part of X2
	      {
		d.tridag_low[ j ]     = arg_A*exp(-2.*I*eta);	     
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*eta)+arg_V;
		d.tridag_upp[ j ]     = arg_A*exp(-2.*I*eta);
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-2.*I*eta)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	    if (j==wf.n1-N_z_right) // right discontinuity point
	      {
		d.tridag_low[ j ]     = arg_A*2./(1.+exp(I*eta));
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*eta)+arg_V;
		d.tridag_upp[ j ]     = arg_A*(2.*exp(-I*eta))/(1.+exp(I*eta));	     
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-I*eta)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	  }
	Tridag( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , d.gam );
	// Tridag_Fast( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , wf.n2 );
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  wf.wave[ wf.in2( 1 , j ) ] = d.wf_1D_solution[ j ];
  } // end of Hamiltonian_X_ECS  

 void Hamiltonian::X_ECS_E_R_back_left(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
    
	  int N_z_right=int(wf.n1*absorb_percentage);
	  int N_z_left=int(0);
	
	for( int j = 0 ; j < wf.n1+2 ; j++ )
	  {
	    d.v_1D[ j ]  = p.ABV_V[ wf.in2( 1 , j ) ];
	    d.wf_1D[ j ] = wf.wave[  wf.in2(  1 , j ) ];
	  }
	
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  {
	    complex<double>  arg_A = I*time_step/( 2.*wf.dx1*wf.dx1 )*p.ABV_A_x;
	    complex<double>  arg_B = 0.;//( I*time_step/( 4.*wf.dx2 )    )*ABV_B_x2;
	    //complex<double>  arg_V = ( I*time_step/2.               )*d.v_1D[ j ];
	    complex<double>  arg_V =   I*time_step/2.               *(d.v_1D[ j ]+wf.x1[j]*electric_field);
	    

	    if (wf.symmetry_x1==0) // X1 starts from -X2_max to X2_max
	      {   // ECS is done in the left and right side of X2
		
		if (j>N_z_left && j<wf.n1-N_z_right) // center of X2
		  {
		    d.tridag_low[ j ]     = arg_A-arg_B;	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A+arg_V;
		    d.tridag_upp[ j ]     = arg_A+arg_B;	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		
	       	if (j<N_z_left)  // left of x2
		  {
		    d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		if (j==N_z_left) // left discontinuity point
		  {
		    d.tridag_low[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));  
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*2./(1.+exp(I*etab));	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
	      } // end if
	    
	
	    
	    if ( j>wf.n1-N_z_right ) // right part of X2
	      {
		d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	    if (j==wf.n1-N_z_right) // right discontinuity point
	      {
		d.tridag_low[ j ]     = arg_A*2./(1.+exp(I*etab));
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));	     
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	  }
	Tridag( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , d.gam );
	// Tridag_Fast( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , wf.n2 );
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  wf.wave[ wf.in2( 1 , j ) ] = d.wf_1D_solution[ j ];
  } // end of Hamiltonian_X_ECS_back_left

  void Hamiltonian::X_ECS_E_R_back_right(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
    
	  int N_z_right=int(wf.n1);
	  int N_z_left=int(wf.n1*absorb_percentage);
	
	for( int j = 0 ; j < wf.n1+2 ; j++ )
	  {
	    d.v_1D[ j ]  = p.ABV_V[ wf.in2( 1 , j ) ];
	    d.wf_1D[ j ] = wf.wave[  wf.in2(  1 , j ) ];
	  }
	
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  {
	    complex<double>  arg_A = I*time_step/( 2.*wf.dx1*wf.dx1 )*p.ABV_A_x;
	    complex<double>  arg_B = 0.;//( I*time_step/( 4.*wf.dx2 )    )*ABV_B_x2;
	    //complex<double>  arg_V = ( I*time_step/2.               )*d.v_1D[ j ];
	    complex<double>  arg_V =   I*time_step/2.               *(d.v_1D[ j ]+wf.x1[j]*electric_field);
	    

	    if (wf.symmetry_x1==0) // X1 starts from -X2_max to X2_max
	      {   // ECS is done in the left and right side of X2
		
		if (j>N_z_left && j<N_z_right) // center of X2
		  {
		    d.tridag_low[ j ]     = arg_A-arg_B;	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A+arg_V;
		    d.tridag_upp[ j ]     = arg_A+arg_B;	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		
	       	if (j<N_z_left)  // left of x2
		  {
		    d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		if (j==N_z_left) // left discontinuity point
		  {
		    d.tridag_low[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));  
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*2./(1.+exp(I*etab));	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
	      } // end if
	    
	
	    
	    if ( j>N_z_right ) // right part of X2
	      {
		d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	    if (j==N_z_right) // right discontinuity point
	      {
		d.tridag_low[ j ]     = arg_A*2./(1.+exp(I*etab));
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));	     
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	  }
	Tridag( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , d.gam );
	// Tridag_Fast( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , wf.n2 );
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  wf.wave[ wf.in2( 1 , j ) ] = d.wf_1D_solution[ j ];
  } // end of Hamiltonian_X_ECS_back_right

  void Hamiltonian::X_ECS_E_R_unequalabsorb(complex<double>  time_step , wavefunction &wf, const ABVparam &p, double electric_field, double absorb_percentage_left, double absorb_percentage_right )
  {
	  temp_data_t &d=temp_data[1]; 
	  if(n[1] != wf.n1 || dx[1] != wf.dx1) {
		  cerr<<"wavefunction has wrong size for Hamiltonian\n";
		  exit(1);
	  }
    
    int N_z_right=int(wf.n1*absorb_percentage_right);
    int N_z_left=int(wf.n1*absorb_percentage_left);
	
	for( int j = 0 ; j < wf.n1+2 ; j++ )
	  {
	    d.v_1D[ j ]  = p.ABV_V[ wf.in2( 1 , j ) ];
	    d.wf_1D[ j ] = wf.wave[  wf.in2(  1 , j ) ];
	  }
	
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  {
	    complex<double>  arg_A = I*time_step/( 2.*wf.dx1*wf.dx1 )*p.ABV_A_x;
	    complex<double>  arg_B = 0.;//( I*time_step/( 4.*wf.dx2 )    )*ABV_B_x2;
	    //complex<double>  arg_V = ( I*time_step/2.               )*d.v_1D[ j ];
	    complex<double>  arg_V =   I*time_step/2.               *(d.v_1D[ j ]+wf.x1[j]*electric_field);
	    

	    if (wf.symmetry_x1==0) // X1 starts from -X2_max to X2_max
	      {   // ECS is done in the left and right side of X2
		
		if (j>N_z_left && j<wf.n1-N_z_right) // center of X2
		  {
		    d.tridag_low[ j ]     = arg_A-arg_B;	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A+arg_V;
		    d.tridag_upp[ j ]     = arg_A+arg_B;	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		
		if (j<N_z_left )  // left of x2
		  {
		    d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
		if (j==N_z_left) // left discontinuity point
		  {
		    d.tridag_low[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));  
		    d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		    d.tridag_upp[ j ]     = arg_A*2./(1.+exp(I*etab));	     
		    d.wf_1D_rightside[ j ]=
		      -d.tridag_low[j]*d.wf_1D[ j-1 ]
		      +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		      -d.tridag_upp[j]*d.wf_1D[ j+1 ];
		  }
	      } // end if
	    
	
	    
	    if ( j>wf.n1-N_z_right) // right part of X2
	      {
		d.tridag_low[ j ]     = arg_A*exp(-2.*I*etab);	     
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-2.*I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*exp(-2.*I*etab);
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-2.*I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	    if (j==wf.n1-N_z_right) // right discontinuity point
	      {
		d.tridag_low[ j ]     = arg_A*2./(1.+exp(I*etab));
		d.tridag_mid[ j ]     = 1.-2.*arg_A*exp(-I*etab)+arg_V;
		d.tridag_upp[ j ]     = arg_A*(2.*exp(-I*etab))/(1.+exp(I*etab));	     
		d.wf_1D_rightside[ j ]=
		  -d.tridag_low[j]*d.wf_1D[ j-1 ]
		  +( 1.+2.*arg_A*exp(-I*etab)-arg_V )*d.wf_1D[ j ]
		  -d.tridag_upp[j]*d.wf_1D[ j+1 ];
	      }
	  }
	Tridag( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , d.gam );
	// Tridag_Fast( d.tridag_low , d.tridag_mid , d.tridag_upp , d.wf_1D_rightside , d.wf_1D_solution , wf.n2 );
	for( int j = 1 ; j <= wf.n1 ; j++ )
	  wf.wave[ wf.in2( 1 , j ) ] = d.wf_1D_solution[ j ];
  } // end of Hamiltonian_X_ECS  

/*------------------  Wavefunction related functions  -----------------------*/

  /**
   * Initialize wavefunction
   */
	void Initialize( wavefunction &wf , int n_point , double spatial_step , int symmetry_x1)
	{
		wf.verbose=false;
		/**
		 * Initialize symmetries of the wavefunction.  
		 */
		wf.symmetry_x1=symmetry_x1;
	
                wf.n1 = n_point;
		
		/**
		 * Allocate the wavefunction and the potential
		 */
		wf.wave.resize(( wf.n1+2 )*( 0+2 ), 0.);
		
		
		/**
		 * Allocate the grid
		 */
		if (abs(symmetry_x1) == 1) wf.x1.type=HalfAxis;
		else wf.x1.type=FullAxis;
		wf.x1.Init("x", n_point, spatial_step, wf.x1.type);
		wf.dp1 = 2*pi/(wf.dx1*wf.n1);
		wf.p1.Init("p_x", n_point, wf.dp1, FFTWAxis);
		

    
		/**
		 * Initialize the spatial grid helper constants.
		 */
		wf.one_by_dx1=1. / wf.dx1;     // wf.dx1=x[1].delta (which itself is set in Axis.Init)
		wf.one_by_dx1sqr=1. / ( wf.dx1*wf.dx1 );
	
	
	} // Initialize_wavefunction
	
	void Initialize_Momentum( wavefunction &wf , wavefunction &wf_mom)
	{
		
		wf_mom.verbose=false;
		
		/**
		 * Initialize symmetries of the wavefunction.  
		 */
		wf_mom.symmetry_x1=wf.symmetry_x1;
	
		//symmetry=0: X starts from -X_max, no symmetry or antisymmetry for wavefunction necessary
		//symmetry=1: X starts from 0, wavefunction is symmetric
		//symmetry=-1: X starts from 0, wavefunction is antisymmetric
		
		/**
		 * Allocate the grid
		 */
		// !!!!!!!!!!! only symmetry =0 !!!!!!!!!!!!
		wf_mom.x1.Init("p_x", wf_mom.n1, wf_mom.dx1, FullAxis);

		wf_mom.wave.resize(( wf_mom.n1+2 )*( 0+2 ), 0.);
		
		
		/**
		 * Allocate the grid
		 */
		wf_mom.x1.resize(wf_mom.n1+2, 0.);
		/**
		 * Define the real grid coordinates
		 */
		
		//    if(wf.symmetry_x1==1 || wf.symmetry_x1==-1)
		/**
		 * symmetric(=1) or antisymmetric(=-1) wavefunction; only half the grid is needed
		 */
		



    // !!!!!!!!!!! only symmetry =0 !!!!!!!!!!!!

		for( int i = 1 ; i <= int(wf_mom.n1/2) +1 ; i++ )
		{	
			wf_mom.x1[ i ] =   (i-1)*wf_mom.dx1;
		}  
		for( int i = int(wf_mom.n1/2)+2 ; i <= wf_mom.n1 ; i++ )
		{	
			wf_mom.x1[ i ] =   (i-1-wf_mom.n1)*wf_mom.dx1;
		}  
		
		
	} // Initialize_momentum
	
	
	/**
	 * Initialize a Gaussian function as a initial wavepacket for imaginary time propagation
	 */
	void Guess_Function_Gauss( wavefunction &wf , double sigma , double gauss_center , double initial_momentum )
	{
		

			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				double arg_x1 = ( wf.x1[ i ]-gauss_center )/( 2.*sigma );
				
				complex<double>  arg_p1 = I*initial_momentum*wf.x1[ i ];
				
				wf.wave[ wf.in2(1,i) ] = exp( -arg_x1*arg_x1 )*exp( arg_p1);
			}

		
		
	} // end of Guess_Function_Gauss
	
        /*another guess function created by Jing on 11/06/2009*/
	void Guess_Function_Uniform( wavefunction &wf , double wavefunction_value )
	{
		
		

			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				wf.wave[ wf.in2(1,i) ] = wavefunction_value;
			}

		
		
	} // end of Guess_Function_Uniform  
  
    /**
     * Normalize of the wavefunction phi' = phi/sqrt( <phi|phi> )
     */
	void Normalize( wavefunction &wf )
	{
		double norm_wf = Obs_Norm( wf );
		//	if(wf.verbose==true)
		//	cout << "Normalizing: Input Norm = " << norm << "\n";
		
		for( int index = 1 ; index <= wf.n1 ; index++ ) 
		{
			wf.wave[ wf.in2(1,index) ] /= norm_wf; 
		}
		
		//norm = Obs_Norm( wf );
		
		//	if(wf.verbose==true)
		//	cout << "Normalizing: Output Norm = " << norm << "\n";
		
	} // end of Normalize3D
	
	
    /**
     * Calculate the norm ( = sqrt( <phi|phi> ) = sqrt( int( x1*dx1*dx2* ( conj( phi )*phi ) ) ) )
     */
	double Obs_Norm( wavefunction &wf )
	{
		double obs=0.;
		

			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += norm( wf.wave[ wf.in2(1,i) ] );
			}

		
		obs *= wf.dx1;
		obs= sqrt( obs );
		return obs;
		
	} // end of Obs_Norm
	

	
	/**
	 * Calculate the energy ( = <phi|H|phi> with the Hamiltonian in Eq. (4) in doc/ABV/abv.tex
	 */
	double Obs_Energy( wavefunction &wf, ABVparam &p)
	{
		double obs=0.;
		

			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += real( conj( wf.wave[ wf.in2(1,i) ] )*(
						     wf.wave[ wf.in2(1,i) ]*( p.ABV_V[ wf.in2(1,i) ] )
						     +Second_Derivative_X( wf , 1 , i )*p.ABV_A_x
						     ) );
			//cout<<Second_Derivative_X( wf , 1 , i )*p.ABV_A_x<<endl;
			}

		obs *= wf.dx1;
		return obs;
		
	} // end of Obs_Energy
	


	/**
	 * Calculate the expectation value <x1> = <phi|x1|phi>
	 */
	double Obs_Expectation_Value_X( wavefunction &wf )
	{
		double obs = 0.;
		

			for( int i = 1 ; i <= wf.n1 ; i++ )
			{
				obs += wf.x1[ i ]
					*norm( wf.wave[ wf.in2(1,i) ] );
			}

	
		obs *= wf.dx1;
		return obs;
		
		
	} // end of Obs_Expectation_Value_X1


    /**
     * Calculate dx1 = <phi|sqrt( <x1^2> -<x1>^2 )|phi>
     */
  double Obs_Expectation_Value_Width_X( wavefunction &wf )
  {
    double obs;
    double exp_x1 = Obs_Expectation_Value_X( wf );
    double exp_x1_2 = 0.;
    
 

	    for( int i = 1 ; i <= wf.n1 ; i++)
	      {
		exp_x1_2 += wf.x1[ i ]*wf.x1[ i ]
				  *norm( wf.wave[ wf.in2(1,i) ] );
	      }

      
    exp_x1_2 *=wf.dx1;
    obs = sqrt( exp_x1_2-exp_x1*exp_x1 );
    return obs;
    
  } // end of Obs_Expectation_Value_Width_X1

   	
    /**
     * Calculate first derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
	complex<double>  First_Derivative_X( wavefunction &wf , int j , int i )
	{
		complex<double>  deriv;
		
		if (wf.symmetry_x1==1 && i ==1 )//symmetry
		{
			
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx1 );
			
		}
		if (wf.symmetry_x1==-1 && i ==1 )//antisymmetry
		{
			
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]+wf.wave[ wf.in2(j,i) ] ) / ( 2.0*wf.dx1 );
			
		}
		else
		{
			deriv = ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(  j , i-1 ) ] ) / ( 2.0*wf.dx1 );
		}//sym==1
		return deriv;
		
	} // end of First_Derivative_X1
	
	

 
  


    /**
     * Calculate second derivative in x1-direction - look at Eq. (5) in doc/ABV/abv.tex
     * note the boundary condition phi( dx1/2 ) = phi( -dx1/2 ) which leads to phi( index = 0 ) = phi( index = 1 )
     */
  complex<double>  Second_Derivative_X( wavefunction &wf , int j , int i )
  {

    if (i==1&&fabs(wf.symmetry_x1)==1) {
    if (wf.symmetry_x1==1)//symmetry
      {
	return ( wf.wave[ wf.in2(  j , i+1 ) ]-wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx1sqr;
      }
    if (wf.symmetry_x1==-1)//antisymmetry
      {
	return ( wf.wave[ wf.in2( j , i+1 ) ]-3.*wf.wave[ wf.in2(j,i) ] ) * wf.one_by_dx1sqr;
      }
    //if (wf.symmetry_x1==1)
    } else {

      return ( wf.wave[ wf.in2(  j , i+1 ) ]-2.*wf.wave[ wf.in2(j,i) ] + wf.wave[ wf.in2(  j , i-1 ) ] ) * wf.one_by_dx1sqr;//sym==1

    }
      
  } // end of Second_Derivative_X1
  



	void PlaceWaveFunction( wavefunction &wf , wavefunction &small_wf)
	{
		
		//int shifter_n1=(wf.n1-small_wf.n1)/2;
		for (int i=1; i<2; i++) {
			if(small_wf.n[i] > wf.n[i]) {
				cerr<<"cannot load bigger wavfunction into small grid"<<endl;
				exit(1);
			}
		}
		
		int shifter_n1=0;
		
		
		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;
		}    
		
		

			for(int i=1;i<=small_wf.n1;i++)
			{
				wf.wave[wf.in2(1,shifter_n1+i)]=small_wf.wave[small_wf.in2(1,i)];
			}
		
}


	complex<double>  Project( wavefunction &wf , wavefunction &wf2)
	{
		
		complex<double>  projection=complex<double> (0.,0.);
		

			for(int i=1;i<=wf.n1;i++)
			{
				projection+= conj(wf.wave[wf.in2(1,i)])*wf2.wave[wf2.in2(1,i)];
			}
		projection*=wf.dx1;
		return projection;
	}
	
	/*
	double Overlap( wavefunction &wf , wavefunction &wf2)
	{
		double overlap_;
		complex<double>  proj = Project ( wf , wf2 );
		overlap_ = sqrt( real( proj*conj( proj ) ) );
		return overlap_;
	}
	*/
  
  //Dipole Operators, by Cory, 5/10/2016
  
  complex<double> Dipole( wavefunction &wf, wavefunction &wf2){
    complex<double> Mu (0.0, 0.0);
    for(int i=1; i<=wf.n1; i++){
      Mu+=-conj(wf.wave[wf.in2(1,i)])*complex<double>(wf.x1[i], 0.0)*wf2.wave[wf2.in2(1,i)];
    }
    Mu*=wf.dx1;
    return Mu;
  }

  complex<double> Dipole_k( wavefunction &wf, double &k){
    //Calculates dipole between bound state and a plane wave continuum
    complex<double> Mu_k (0.0, 0.0);
    for(int i=1; i<=wf.n1; i++){
      Mu_k+=-(wf.wave[wf.in2(1,i)])*complex<double>(wf.x1[i], 0.0)*exp(complex<double>(0.0, -k*wf.x1[i]));
    }
    Mu_k*=wf.dx1;
    return Mu_k;
  }
	complex<double>  Project_Diff_Sizes( wavefunction &wf , wavefunction &small_wf)
	{
		int shifter_n1=0;

		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;
		}    
		
		
		complex<double>  projection=complex<double> (0.,0.);
		
		

			for(int i=1;i<=small_wf.n1;i++)
			{
				projection+= conj(small_wf.wave[small_wf.in2(1,i)])
					*wf.wave[wf.in2(1,shifter_n1+i)];
			}
		projection*=small_wf.dx1;
		return projection;
	}
	
	
	void ProjectOUT_Diff_Sizes( wavefunction &wf , wavefunction &small_wf)
	{
		int shifter_n1=wf.n1-small_wf.n1;
	
		
		if (wf.symmetry_x1==0)
		{
			shifter_n1=(wf.n1-small_wf.n1)/2;
		}    

		
		complex<double>  projection=complex<double> (0.,0.);   
		

			for(int i=1;i<=small_wf.n1;i++)
			{
				projection+= conj(small_wf.wave[small_wf.in2(1,i)])
					*wf.wave[wf.in2(1,shifter_n1+i)];
			}
		projection*=small_wf.dx1;
		
		

			for(int i=1;i<=small_wf.n1;i++)
			{
				wf.wave[wf.in2(1,shifter_n1+i)]-=projection*small_wf.wave[small_wf.in2(1,i)];
			}
		
		//    return projection;
	}
	
  Projection1D& Prj_Wavefunction_X(wavefunction &wf, Projection1D& prj,  string title)
  {
    prj.title=title;
    prj.axis_check(wf.x1);
    double sum;
	    for( int i = 1 ; i <= wf.n1+1 ; i++) { 
		    sum = norm(wf.wave[ wf.in2(1,i) ]);
		    prj.data[1*(wf.n1+2)+i]=sum;
	    }
    
    return prj;
  }
        //Mask function centered at x=0
	void Mask_Function_Inner(wavefunction &wf, double x_left, double x_right, double exponent )
	{
		
		
		double argument_x1_right;

		double argument_x1_left;
		
		double mask_x1_right;

		double mask_x1_left;
		
		
			for(int i=1;i<=wf.n1;i++)
			{
				argument_x1_right=(pi/2.)*(x_right-wf.x1[i])/(x_right+1.e-20);
			
				
				argument_x1_left=(pi/2.)*(-x_left+wf.x1[i])/(-x_left+1.e-20);
			
				
				mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
			
				
				mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
			
				
				
				if (wf.symmetry_x1==0)
				{
					if (wf.x1[i]>=x_left&&wf.x1[i]<0)		  
						wf.wave[wf.in2(1,i)]*=mask_x1_left;		  
				}//if sym==1 or sym==-1, we needn't this part.
				
				if (wf.x1[i]<=x_right&&wf.x1[i]>=0)
				{
				        //cout<<wf.wave[wf.in2(j,i)]<<"\t";
					wf.wave[wf.in2(1,i)]*=mask_x1_right;
					//cout<<mask_x1_right<<"\t"<<wf.wave[wf.in2(j,i)]<<endl;
				}
				
			}
	}  
  
	void Mask_Function(wavefunction &wf, double frac_n1_right, double frac_n1_left, double exponent )
	{
		if ( frac_n1_right>1. || frac_n1_right<0. )
		{
			cout<< "frac_n1_right is wrong";
			exit(1);
		}
		
		
		if ( frac_n1_left>1. || frac_n1_left<0. )
		{
			cout<< "frac_n1_left is wrong";
			exit(1);
		}
		
		double mask_start_x1_right=wf.x1[int((wf.n1)*(1.-frac_n1_right))];
		double mask_start_x1_left=wf.x1[int((wf.n1)*frac_n1_left)+1];
		
		
		double argument_x1_right;

		double argument_x1_left;
		
		double mask_x1_right;

		double mask_x1_left;
		
		
			for(int i=1;i<=wf.n1;i++)
			{
				argument_x1_right=(pi/2.)*(wf.x1[i]-mask_start_x1_right)/(wf.x1[wf.n1]-mask_start_x1_right+1.e-20);
			
				
				argument_x1_left=(pi/2.)*(wf.x1[i]-mask_start_x1_left)/(wf.x1[1]-mask_start_x1_left+1.e-20);
			
				
				mask_x1_right=pow(fabs(cos(argument_x1_right)),exponent);
			
				
				mask_x1_left=pow(fabs(cos(argument_x1_left)),exponent);
			
				
				
				if (wf.symmetry_x1==0)
				{
					if (i< int(wf.n1*frac_n1_left))		  
						wf.wave[wf.in2(1,i)]*=mask_x1_left;		  
				}//if sym==1 or sym==-1, we needn't this part.
				
				if (i> int(wf.n1*(1.-frac_n1_right)))
				{
				        //cout<<wf.wave[wf.in2(j,i)]<<"\t";
					wf.wave[wf.in2(1,i)]*=mask_x1_right;
					//cout<<mask_x1_right<<"\t"<<wf.wave[wf.in2(j,i)]<<endl;
				}
				
			}
	}
 
	void Mask_Function_Part(wavefunction &wf, double x_left, double x_right, double exponent )
	{
		
		double argument_x1;
		double mask_x1;

		
			for(int i=1;i<=wf.n1;i++)
			{
				argument_x1=(pi/2.)*(-wf.x1[i]+x_right)/(x_right-x_left+1.e-20);
				mask_x1=pow(fabs(cos(argument_x1)),exponent);
			
				if (wf.x1[i]>=x_left&&wf.x1[i]<=x_right)
				{
					wf.wave[wf.in2(1,i)]*=mask_x1;
				}
				
			}
	}


	double Overlap( wavefunction &wf , wavefunction &wf2)
	{
		double overlap_;
		complex<double>  proj = Project ( wf , wf2 );
		//overlap_ = sqrt( real( proj*conj( proj ) ) );
		 overlap_ =  real( proj*conj( proj ) );
		return overlap_;
	}

  //FFT1D written by Cory 9/6/2013

  void FFT1D(wavefunction &wf, wavefunction &wf_transformed){
    fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wf.n1);
    fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * wf.n1);
    fftw_plan p1 = fftw_plan_dft_1d(wf.n1, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for(int i=1; i<=wf.n1; i++){
      int fft_index = i-1;
      //in[fft_index] = reinterpret_cast<fftw_complex>(wf.wave[wf.in2(1,i)]);
      // memcpy(in[fft_index], wf.wave[wf.in2(1,i)], sizeof(fftw_complex));
      in[fft_index][0] = real(wf.wave[wf.in2(1,i)]);
      in[fft_index][1] = imag(wf.wave[wf.in2(1,i)]);
    }
    fftw_execute(p1); 
    double fac = 2*pi/(wf_transformed.n1*wf_transformed.dx1);
    for(int i=1; i<=wf.n1; i++){
      int fft_index = i-1;
      //wf_transformed.wave[wf_transformed.in2(1,i)] = reinterpret_cast<fftw_complex*>(fac*(out[fft_index]));
      //memcpy(wf_transformed.wave[wf_transformed.in2(1,i)], fac*(out[fft_index], sizeof(fftw_complex)));
      wf_transformed.wave[wf_transformed.in2(1,i)] = fac*(out[fft_index][0] + I*out[fft_index][1]);
    }
         fftw_destroy_plan(p1);
         fftw_free(in);
	 fftw_free(out);
    }


} // end of namespace Cartesian_1D
