// $Id: Laser.h 315 2006-07-18 21:41:15Z arvid $
#include <cstdlib>
#include <cmath> 
#ifndef POTENTIALS_H
#define POTENTIALS_H

using namespace std;

namespace Cartesian_3D{

  struct ABVparam {
  
    /**
     * Declare ABV coefficients.
     */
    double ABV_A_x1;
    double ABV_B_x1;		// different from Cylindrical3D !!     Everything has to be changed and checked accordingly!
  
    double ABV_A_x2;
    double ABV_B_x2;
  
    double ABV_A_x3;
    double ABV_B_x3;
  
    /**
     * Declare potentials.
     */
    vector<double> ABV_V;
    vector<double> ABV_V_deriv;
  
  };


  /***********************************************************************************************/
  /***********************************************************************************************/
  // Potential SAE

  struct Potential_SAE_noble : ABVparam
  {
    double SAEparam_0;
    double SAEparam_1;
    double SAEparam_2;
    double SAEparam_3;
    double SAEparam_4;
    double SAEparam_5;
    double SAEparam_6;
    double mass;
    double electron_charge; //Charge in ATOMIC UNITS
    //Potential taken from the Paper JPB Lin 38, (2005) 2593
    Potential_SAE_noble(const wavefunction &wf) {
      SAEparam_0=1.;
      SAEparam_1=1.231;
      SAEparam_2=0.662;
      SAEparam_3=-1.325;
      SAEparam_4=1.236;
      SAEparam_5=-0.231;
      SAEparam_6=0.480;
      mass=1;
      electron_charge=-1.;
      updateABV(wf);
    };// Here the default values are for Helium
    
    Potential_SAE_noble(const wavefunction &wf,vector<double> SAE_param, double _mass, double _electron_charge ) 
    {

      SAEparam_0= SAE_param[0];
      SAEparam_1= SAE_param[1];
      SAEparam_2= SAE_param[2];
      SAEparam_3= SAE_param[3];
      SAEparam_4= SAE_param[4];
      SAEparam_5= SAE_param[5];
      SAEparam_6= SAE_param[6];
      mass=_mass;
      electron_charge=_electron_charge;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = electron_charge/mass;
      ABV_B_x2 = electron_charge/mass;
      ABV_B_x3 = electron_charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      //Potential taken from the Paper JPB Lin 38, (2005) 2593
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		double r=sqrt( wf.x1[ i ]*wf.x1[ i ]+ wf.x2[ j ]*wf.x2[ j ]+wf.x3[ k ]*wf.x3[ k ]);
		ABV_V[ wf.in3( k , j , i ) ] = -(SAEparam_0+SAEparam_1*exp(-r*SAEparam_2)+SAEparam_3*r*exp(-r*SAEparam_4)+SAEparam_5*exp(-r*SAEparam_6))/r;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };
  
  
  //user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential: ABVparam
  {
    double mass;
    double electron_charge;
    int n1, n2, n3;
    
    user_input_potential(const vector<double> &potential_array, int _n1, int _n2, int _n3) 
    {
      mass=1.0;
      electron_charge=-1.0;   
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = electron_charge/mass;
      ABV_B_x2 = electron_charge/mass;
      ABV_B_x3 = electron_charge/mass;
      n1=_n1;
      n2=_n2;
      n3=_n3;
 	          
      ABV_V.resize(( n1+2 )*( n2+2 )*( n3+2 ), 0.);
      
      for( int k = 1 ; k <= n3 ; k++ )
      {
	for( int j = 1 ; j <= n2 ; j++ )
	  {
	    for( int i = 1 ; i <= n1 ; i++ )
	      {
		ABV_V[ k*(n2+2)*(n1+2)+j*(n1+2)+i ] = potential_array[(i-1)*n2*n3+(j-1)*n3+(k-1)];
		ABV_V[ k*(n2+2)*(n1+2)+j*(n1+2)+i ] /= 3.;		
	      }
	  }
      }
     };
	  
  };
  

  /***********************************************************************************************/
  /***********************************************************************************************/
  // Potential SAE Muller FOAM

  struct Potential_SAE_muller : ABVparam
  {
    double SAEparam_0;
    double mass;
    double electron_charge; //Charge in ATOMIC UNITS

    Potential_SAE_muller(const wavefunction &wf) {
      SAEparam_0=-2.5;
      mass=1.;
      electron_charge=-1.;
      updateABV(wf);
    };
    
    Potential_SAE_muller(const wavefunction &wf,double SAE_param, double _mass, double _electron_charge ) 
    {

      SAEparam_0= SAE_param;
      mass=_mass;
      electron_charge=_electron_charge;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 =  electron_charge/mass;
      ABV_B_x2 =  electron_charge/mass;
      ABV_B_x3 =  electron_charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);


     
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		double r=sqrt( wf.x1[ i ]*wf.x1[ i ]+ wf.x2[ j ]*wf.x2[ j ]+wf.x3[ k ]*wf.x3[ k ]);
		ABV_V[ wf.in3( k , j , i ) ] = -1./r-(2.+1./r)*exp(-4.*r);
		complex<double>  temp= ABV_V[ wf.in3( k , j , i ) ];
		if (real(temp) < SAEparam_0)
		  {
		    ABV_V[ wf.in3( k , j , i ) ] = SAEparam_0;
		  }
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };
 

  struct Potential_H2_plus_e3D_n2D: ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double nuc_x1,nuc_x2,nuc_y1,nuc_y2;
    double deriv_a;
    double deriv_b;
    double deriv_c;
    double deriv_d;
    double v_x1;
    double v_y1;
    double v_x2;
    double v_y2;

    Potential_H2_plus_e3D_n2D(const wavefunction &wf) 
    {
      charge = -1.;
      charge_nucleus_1=charge_nucleus_2 = 1.;
      mass   = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      softening_nn = 0.01;
      softening_en = 0.001;
      nuc_x1=-1;
      nuc_x2=1;
      nuc_y1=nuc_y2=0;
      deriv_a=deriv_b=deriv_c=deriv_d=0;
      v_x1=v_x2=v_y1=v_y2=0;
    };

    Potential_H2_plus_e3D_n2D(const wavefunction &wf,double _charge, double _charge_nucleus_1,double _charge_nucleus_2,
			      double _mass, double _mass_nucleus_1,double _mass_nucleus_2,double _softening_en,double _softening_nn, 
			      vector<double> nuc_position,vector<double> deriv, vector<double> nuc_velocity)
    {      
      charge= _charge;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass =_mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      nuc_x1=nuc_position[0];
      nuc_y1=nuc_position[1];
      nuc_x2=nuc_position[2];
      nuc_y2=nuc_position[3];

      deriv_a=deriv[0];
      deriv_b=deriv[1];
      deriv_c=deriv[2];
      deriv_d=deriv[3];

      v_x1=nuc_velocity[0];
      v_y1=nuc_velocity[1];
      v_x2=nuc_velocity[2];
      v_y2=nuc_velocity[3];

      softening_en = _softening_en;
      softening_nn = _softening_nn;

      updateABV(wf);
    };

    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = charge/mass;
      ABV_B_x2 = charge/mass;
      ABV_B_x3 = charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      // the two soft core parameters are 0.1 and 0.0109
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		/**************   3D  ******************/
	 	ABV_V[ wf.in3( k , j , i ) ] =  		  
 		  charge*charge_nucleus_1/sqrt(wf.x3[k]*wf.x3[k]+(wf.x1[i]-nuc_x1)*(wf.x1[i]-nuc_x1)+
						       (wf.x2[j]-nuc_y1)*(wf.x2[j]-nuc_y1)+softening_en ) +
 		  charge*charge_nucleus_2/sqrt(wf.x3[k]*wf.x3[k]+(wf.x1[i]-nuc_x2)*(wf.x1[i]-nuc_x2)+(wf.x2[j]-nuc_y2)*(wf.x2[j]-nuc_y2)+softening_en );
 		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
		
		/**************   2D   **************
		 //softening_en=0.31.
		ABV_V[ wf.in3( k , j , i ) ] = 
		  charge*charge_nucleus_1/sqrt((wf.x1[i]-nuc_x1)*(wf.x1[i]-nuc_x1)+
		  				       (wf.x2[j]-nuc_y1)*(wf.x2[j]-nuc_y1)+softening_en)+
		  charge*charge_nucleus_2/sqrt((wf.x1[i]-nuc_x2)*(wf.x1[i]-nuc_x2)+
						       (wf.x2[j]-nuc_y2)*(wf.x2[j]-nuc_y2)+softening_en);
		
		
		ABV_V[ wf.in3( k , j , i ) ] /=2.;*/

	      }
	  }
      }
    }


  };



  /***********************************************************************************************/
  /***********************************************************************************************/

  struct Potential_H2_2D : ABVparam
  {
    double charge_1;
    double charge_2;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass_1;
    double mass_2;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double internuclear_distance;
    double softening_en;
    double softening_ee;

    /* Constructor to initialize the default values: */
    Potential_H2_2D(const wavefunction &wf) 
    {
      /* what to do here ?? */
      charge_1=charge_2 = -1.;
      charge_nucleus_1=charge_nucleus_2 = 1.;
      mass_1 = mass_2 = 1.;
      mass_nucleus_1 = mass_nucleus_2 = 1836.;
      internuclear_distance=1.4;
      softening_en = 0.7;
      softening_ee = 1.2375;
      updateABV(wf);
    };

    Potential_H2_2D(const wavefunction &wf,double _charge_1,double _charge_2, double _charge_nucleus_1,double _charge_nucleus_2,
			double _mass_1, double _mass_2, double _mass_nucleus_1,double _mass_nucleus_2,double _internuclear_distance, 
			double _softening_en,double _softening_ee  ) 
    {      
      charge_1= _charge_1;
      charge_2 = _charge_2 ;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass_1 = _mass_1;
      mass_2 = _mass_2;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      internuclear_distance= _internuclear_distance;
      softening_en = _softening_en;
      softening_ee = _softening_ee;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */

      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 );
      double center_of_mass = mass_1 + mass_2;

      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/(mass_1+mass_2);
      double center_of_mass_charge =(charge_1+charge_2);
      
      ABV_A_x1 = 0.;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;
      ABV_B_x3 = 0.;
      

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = charge_1*charge_2/sqrt( wf.x2[ j ]*wf.x2[ j ] +softening_ee)
		  +charge_1*charge_nucleus_1/sqrt( ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		  +charge_2*charge_nucleus_1/sqrt( ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		  +charge_1*charge_nucleus_2/sqrt( ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		  +charge_2*charge_nucleus_2/sqrt( ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 2.;
	      }
	  }
      }
    }
  };



  struct Potential_H2_planar : ABVparam
  {
    double charge_1;
    double charge_2;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass_1;
    double mass_2;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double internuclear_distance;
    double angle_nuclei;
    double softening_en;
    double softening_ee;

    /* Constructor to initialize the default values: */
    Potential_H2_planar(const wavefunction &wf) 
    {
      /* what to do here ?? */
      charge_1=charge_2 = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      mass_1 = mass_2 = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      internuclear_distance=1.4;
      softening_en = 0.375;
      softening_ee = 0.;
      angle_nuclei=0.;
      updateABV(wf);
    };

    Potential_H2_planar(const wavefunction &wf,double _charge_1,double _charge_2, double _charge_nucleus_1,double _charge_nucleus_2,
			double _mass_1, double _mass_2, double _mass_nucleus_1,double _mass_nucleus_2,double _internuclear_distance, 
			double _angle_nuclei,double _softening_en,double _softening_ee  ) 
    {      
      charge_1= _charge_1;
      charge_2 = _charge_2 ;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass_1 =_mass_1;
      mass_2 = _mass_2;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      internuclear_distance= _internuclear_distance;
      angle_nuclei=_angle_nuclei;
      softening_en = _softening_en;
      softening_ee = _softening_ee;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */

      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 );
      double center_of_mass = mass_1 + mass_2;

      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/(mass_1+mass_2);
      double center_of_mass_charge =(charge_1+charge_2);
      
      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;
      ABV_B_x1 = reduced_charge/reduced_mass;
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = charge_1*charge_2/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ]+ softening_ee )
		  +charge_1*charge_nucleus_1/sqrt( 0.25*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_2*charge_nucleus_1/sqrt( 0.25*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_1*charge_nucleus_2/sqrt( 0.25*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_2*charge_nucleus_2/sqrt( 0.25*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };



  struct Potential_H2_plus_with_two_electrons_planar : ABVparam
  {
    double charge_1;
    double charge_2;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass_1;
    double mass_2;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double internuclear_distance;
    double angle_nuclei;
    double softening_en;

    Potential_H2_plus_with_two_electrons_planar(const wavefunction &wf) 
    {
      charge_1=charge_2 = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      mass_1 = mass_2 = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      internuclear_distance=1.4;
      softening_en = 0.375;
      angle_nuclei=0.;
      updateABV(wf);
    };

    Potential_H2_plus_with_two_electrons_planar(const wavefunction &wf,double _charge_1,double _charge_2, double _charge_nucleus_1,double _charge_nucleus_2,
			double _mass_1, double _mass_2, double _mass_nucleus_1,double _mass_nucleus_2,double _internuclear_distance, 
			double _angle_nuclei,double _softening_en  ) 
    {      
      charge_1= _charge_1;
      charge_2 = _charge_2 ;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass_1 =_mass_1;
      mass_2 = _mass_2;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      internuclear_distance= _internuclear_distance;
      angle_nuclei=_angle_nuclei;
      softening_en = _softening_en;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */

      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 );
      double center_of_mass = mass_1 + mass_2;

      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/(mass_1+mass_2);
      double center_of_mass_charge =(charge_1+charge_2);
      
      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;
      ABV_B_x1 = reduced_charge/reduced_mass;
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = charge_1*charge_nucleus_1/sqrt( 0.25*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_2*charge_nucleus_1/sqrt( 0.25*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]+internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_1*charge_nucleus_2/sqrt( 0.25*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  +charge_2*charge_nucleus_2/sqrt( 0.25*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )*( wf.x1[ i ]-internuclear_distance*sin( angle_nuclei ) )+
		     ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance*cos( angle_nuclei ) )+softening_en )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };
  
  struct Potential_H_like_e3D : ABVparam
  {
    double charge;
    double charge_nucleus;
    double screened_charge_factor;
    double mass;
    

    /* Constructor to initialize the default values: */
    Potential_H_like_e3D(const wavefunction &wf) 
    {
      charge=-1;
      charge_nucleus=1.;
      screened_charge_factor=1;
      mass=1.;
      updateABV(wf);
    };

    Potential_H_like_e3D(const wavefunction &wf, double _charge, double _charge_nucleus,double _screened_charge_factor,double _mass) 
    {
      charge=_charge;
      charge_nucleus=_charge_nucleus;
      screened_charge_factor=_screened_charge_factor;
      mass=_mass;      
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = charge/mass;
      ABV_B_x2 = charge/mass;
      ABV_B_x3 = charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      // the two soft core parameters are 0.1 and 0.0109
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = -1/sqrt( wf.x1[i]*wf.x1[i]+ wf.x2[j]*wf.x2[j]+wf.x3[k]*wf.x3[k] );
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };

  struct Potential_H_like_e3minus2D : ABVparam
  /*
   * Reduced dimensionality Howto:
   * 1. set ABV_A_x and ABV_B_x to zero for the dummy axes
   * 2. in updateABV:
   *    a) devide ABV_V only by the reduced dimension (e.g. 3-2=1) 
   *    b) remove unused axis from the length (x^2 +y^2 +z^2) in the potential
   * 3. in Obs_Energy_3minus2D (Cartesian3D.cpp) multiply p.ABV_V only with
   *    the reduced dimension.
   * 4. in the cfg set a) n_points = ( 1, 1, Number )
   *                   b) spatial_step = ( 1, 1, delta )
   * 5. adjust the softening_en to match the energy of the state of interest
   */
  {
    double charge;
    double charge_nucleus;
    double softening_en;
    double mass;
    
    /* Constructor to initialize the default values: */
    Potential_H_like_e3minus2D(const wavefunction &wf)
    {
      charge=-1;
      charge_nucleus=1.;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.92; //tested for dx=0.1 to 0.3 a.u.
      mass=1.;
      updateABV(wf);
    };

    Potential_H_like_e3minus2D(const wavefunction &wf, double _charge, double _charge_nucleus, double _softening_en, double _mass) 
    {
      charge=_charge;
      charge_nucleus=_charge_nucleus;
      softening_en = _softening_en;
      mass=_mass;      
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 = 0;
      ABV_A_x2 = 0;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = 0;
      ABV_B_x2 = 0;
      ABV_B_x3 = charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      // the two soft core parameter is 0.0109
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = charge*charge_nucleus/sqrt( wf.x3[k]*wf.x3[k] + softening_en );
	      }
	  }
      }
    }
  };

  /*
 struct H_Lorentz : ABVparam
  {
    double charge;
    double charge_nucleus;
    double screened_charge_factor;
    double mass;
    double mass_nucleus;


    Potential_H_like_e3D(const wavefunction &wf) {
     

      constructABV(wf);
    };

    void constructABV(const wavefunction &wf) {
    
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;
      ABV_B_x3 = 0.;

      // following code is from code from Initialize_grid //
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      // following code is from code from Initialize_Potential //
      // the two soft core parameters are 0.1 and 0.0109
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = +charge*charge_nucleus*screened_charge_factor/sqrt( wf.x1[i]*wf.x1[i]+ wf.x2[j]*wf.x2[j]+wf.x3[k]*wf.x3[k] )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }


 void updateABV(double time,double gamma, double impact_param,double charge_projectile, const wavefunction &wf) {
      // following code is from code from Create_Potential //
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;
      ABV_B_x3 = 0.;

      // following code is from code from Initialize_grid //
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      // following code is from code from Initialize_Potential //
      // the two soft core parameters are 0.1 and 0.0109
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = +charge*charge_nucleus*screened_charge_factor/sqrt( wf.x1[i]*wf.x1[i]+ wf.x2[j]*wf.x2[j]+wf.x3[k]*wf.x3[k] )+charge*charge_projectile*gamma/sqrt(gamma*gamma*(wf.x1[i]+x0-v_projectile*time)*(wf.x1[i]+x0-v_projectile*time)+(wf.x2[j]-impact_param)* (wf.x2[j]-impact_param)+wf.x3[k]*wf.x3[k]  )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }

  };

  */



  
  struct Potential_Free : ABVparam
  {
    double mass;
    double charge;
    /* Constructor to initialize the default values: */
    Potential_Free(const wavefunction &wf) 
    {
      mass=1.;
      charge=-1.;
      updateABV(wf);
    };

    Potential_Free(const wavefunction &wf,double _mass, double _charge ) 
    {
      mass=_mass;
      charge=_charge;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_A_x3 =-0.5/mass;
      ABV_B_x1 = charge/mass;
      ABV_B_x2 = charge/mass;
      ABV_B_x3 = charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = 0.;
	      }
	  }
      }
    }
  };
  
  
  
} // end of namespace Cartesian3D





namespace Cylindrical_3D
{
  struct ABVparam {
  
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x1;
    vector<double> ABV_B_x1;
  
    double ABV_A_x2;
    double ABV_B_x2;
  
    double ABV_A_x3;
    double ABV_B_x3;
  
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;
    vector<double> ABV_V_deriv;

  };
  /***********************************************************************************************/
  /***********************************************************************************************/
  //Example potential in Cylindrical 3D

  struct Potential_Example: ABVparam
  {

    /* Little commet on what this potential mean: This potential is and example and does nothing */

    double parameter1;   /*Parameters from the potential itsels i.e. V(parameter1,parameter2 ; x1,x2,x3) */
    double parameter2;

    /*Charge of the particle in coordinate 1 example; electron 1n 3D x1=rho charge1=charge_electron1=-1 
      Another example x1=rho in Dresden model charge1=reduced_charge_of_two_electrons=(m1*q1-m2*q2)/(m1+m2)=0; */
    
    double charge1;   
    double charge2;
    double charge3;

    /*Mass of the particle in coordinate 1 example; electron 1n 3D x1=rho mass1=mass_electron1=1 
      Another example x1=rho in Dresden model mass1=reduced_mass_of_two_electrons=(m1*m2)/(m1+m2)=1/2; */
    
    double mass1;
    double mass2;
    double mass3;

    /* Constructor to initialize the default values: */
    Potential_Example(const wavefunction &wf) {
      parameter1=1;
      parameter2=1;
      charge1=-1; /* ...etc */
      mass1=1; /* ...etc */
      
      /************************************/
      /* THIS PART NEVER CHANGES in Cyl3D */

      ABV_A_x1 =-0.5/mass1;
      ABV_A_x2 =-0.5/mass2;
      ABV_A_x3 =-0.5/mass3;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = charge2/mass2;
      ABV_B_x3 = charge3/mass3;
      /*********************************/
      updateABV(wf);
    }
    /* Constructor to define paramters from the cfg file: */

    Potential_Example(const wavefunction &wf, double input1,double q1,double m1 /* as many as you want */) {
      parameter1=input1;
      
      charge1=q1; /* ...etc */
      mass1=m1; /* ...etc */
      
      /************************************/
      /* THIS PART NEVER CHANGES in Cyl3D */
      ABV_A_x1 =-0.5/mass1;
      ABV_A_x2 =-0.5/mass2;
      ABV_A_x3 =-0.5/mass3;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = charge2/mass2;
      ABV_B_x3 = charge3/mass3;
      /*********************************/
      updateABV(wf);
    }

    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
     

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      // softening_en=0.0109
      // softening_nn=0.1
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] =parameter1*wf.x1[i]; /* whatever */
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  
  }; //End Example potential in Cylindrical 3D

  /***********************************************************************************************/
  /***********************************************************************************************/

  struct Potential_H2_plus_e2D_R1D : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e2D_R1D(const wavefunction &wf) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 0.0109;
      //softening_en = 0.001;
      softening_nn = 0.1;
      //softening_nn = 0.001;
      updateABV(wf);
    };
    
    Potential_H2_plus_e2D_R1D(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      updateABV(wf);
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      double reduced_mass =  mass_nucleus_1*mass_nucleus_2/( mass_nucleus_1 + mass_nucleus_2 );
      double reduced_charge =(mass_nucleus_2*charge_nucleus_1-mass_nucleus_1*charge_nucleus_2)/( mass_nucleus_1 + mass_nucleus_2 );
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = charge/mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      // softening_en=0.0109
      // softening_nn=0.1
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = screened_charge_factor*charge_nucleus_1*charge_nucleus_2*( 1./sqrt( wf.x2[ j ]*wf.x2[ j ]+softening_nn ) )
		  +charge*charge_nucleus_1/sqrt( wf.x1[ i ]*wf.x1[ i ]+( wf.x3[ k ]+0.5*wf.x2[ j ])*( wf.x3[ k ]+0.5*wf.x2[ j ])+softening_en )
		  +charge*charge_nucleus_2/sqrt( wf.x1[ i ]*wf.x1[ i ]+( wf.x3[ k ]-0.5*wf.x2[ j ])*( wf.x3[ k ]-0.5*wf.x2[ j ])+softening_en )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };

  /***********************************************************************************************/
  /***********************************************************************************************/
  
  struct Potential_H2 : ABVparam
  {
    double charge_1;
    double charge_2;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass_1;
    double mass_2;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double internuclear_distance;
    double softening_en;
    double softening_ee;

    /* Constructor to initialize the default values: */
    Potential_H2(const wavefunction &wf) 
    {
      /* what to do here ?? */
      charge_1=charge_2 = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      mass_1 = mass_2 = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      internuclear_distance=1.4;
      softening_en = 0.375;
      softening_ee = 0.;
      updateABV(wf);
    };

    Potential_H2(const wavefunction &wf,double _charge_1,double _charge_2, double _charge_nucleus_1,double _charge_nucleus_2,
			double _mass_1, double _mass_2, double _mass_nucleus_1,double _mass_nucleus_2,double _internuclear_distance, 
			double _softening_en,double _softening_ee  ) 
    {      
      charge_1= _charge_1;
      charge_2 = _charge_2 ;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass_1 =_mass_1;
      mass_2 = _mass_2;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      internuclear_distance= _internuclear_distance;
      softening_en = _softening_en;
      softening_ee = _softening_ee;
      updateABV(wf);
    };

    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 );
      double center_of_mass = mass_1 + mass_2;
      
      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 );
      double center_of_mass_charge =  charge_1+charge_2;

      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;

      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = charge_1*charge_2/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ] )
		  +charge_1*charge_nucleus_1/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
							   ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		  +charge_2*charge_nucleus_1/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
							   ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		  +charge_1*charge_nucleus_2/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
							   ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		  +charge_2*charge_nucleus_2/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
							   ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		  ;
		ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };
  
  struct Potential_H2_deriv : Potential_H2
  {


    /* No new variables here */

    /* Constructor to initialize the default values: */
    Potential_H2_deriv(const wavefunction &wf) :
      /* base class constructor would get automatically called
       * but it needs an argument.. */
      Potential_H2(wf)
    {
      /* nothing additional to do here */
      updateABV(wf);
    };

    void updateABV(const wavefunction &wf) {
      // call the update for the parent class:
      Potential_H2::updateABV(wf);

      // and initialize the ABV_deriv data

      /* following code is from code from Initialize_grid */
      ABV_V_deriv.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.); 

      /* following code is from code from Initialize_Potential_deriv */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V_deriv[ wf.in3( k , j , i ) ] =
		  charge_1*charge_nucleus_1*0.5*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )/
		  ( sqrt(
			 ( ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 ) )
		  +
		  charge_2*charge_nucleus_1*0.5*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )/
		  ( sqrt(
			 ( ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 ) )
		  -
		  charge_1*charge_nucleus_2*0.5*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )/
		  ( sqrt(
			 ( ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 ) )
		  -
		  charge_2*charge_nucleus_2*0.5*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )/
		  ( sqrt(
			 ( ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 *( ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+wf.x1[ i ]*0.25*wf.x1[ i ]+softening_en )
			 ) )
		  ;
		
		ABV_V_deriv[ wf.in3( k , j , i ) ] /= 3.;
		
	      }
	  }
      }
    }
  };

  /***********************************************************************************************/
  /***********************************************************************************************/
  
  struct Potential_H2_plus_with_two_electrons : ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double internuclear_distance;
    double softening_en;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_with_two_electrons(const wavefunction &wf) 
    {
      /* what to do here ?? */
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      internuclear_distance=1.4;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 0.415;

      updateABV(wf);
    };

    Potential_H2_plus_with_two_electrons(const wavefunction &wf, double _charge, double _charge_nucleus_1, double _charge_nucleus_2,
			double _mass, double _mass_nucleus_1,double _mass_nucleus_2,double _internuclear_distance, 
			double _softening_en ) 
    {


      charge= _charge;
      charge_nucleus_1= _charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      mass =_mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      internuclear_distance= _internuclear_distance;
      softening_en = _softening_en;
      updateABV(wf);
    };

    void updateABV(const wavefunction &wf) 
    {
      /* following code is from code from Create_Potential */
      
      double reduced_mass =  mass*mass/( mass + mass );
      double center_of_mass = mass + mass;
    
      double reduced_charge =  (mass*charge-mass*charge)/( mass + mass );
      double center_of_mass_charge =  charge+charge;
      
      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;
    
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
	{
	  ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
	}
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      
      
      
      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
	{
	  for( int j = 1 ; j <= wf.n2 ; j++ )
	    {
	      for( int i = 1 ; i <= wf.n1 ; i++ )
		{
		  ABV_V[ wf.in3( k , j , i ) ] =  charge*charge_nucleus_1/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
										( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		    +charge*charge_nucleus_1/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
						   ( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]+0.5*internuclear_distance )+softening_en )
		    +charge*charge_nucleus_2/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
						   ( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]+0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		    +charge*charge_nucleus_2/sqrt( 0.25*wf.x1[ i ]*wf.x1[ i ]+
						   ( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )*( wf.x3[ k ]-0.5*wf.x2[ j ]-0.5*internuclear_distance )+softening_en )
		  ;
		  ABV_V[ wf.in3( k , j , i ) ] /= 3.;
		}
	    }
	}
    }
  };
  
  /***********************************************************************************************/
  /***********************************************************************************************/
  
  //Potential for helium in the Dresden model, see C Ruiz PRL 2006
  
  struct Potential_He_dresden : ABVparam
  {
    /* Two electron in 3D the dresden model x1=rho x2=zr x3=zs, Helium problem */
    
    double charge_1;
    double charge_2;
    double charge_nucleus;
    double mass_1;
    double mass_2;
    double softening_en;
    
    /* Constructor to initialize the default values: */
    Potential_He_dresden( wavefunction &wf) {
      charge_1=charge_2 = -1 ;
      charge_nucleus = 2;
      mass_1 = mass_2 = 1;
      softening_en = 0.135; //valid for dx=0.3 a.u.

      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
      
      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */

      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      
      updateABV(wf);

      wf.log << "***************************************************\n";
      wf.log << "Default Constructor Potential_He_dresden"<<"\n";
      wf.log << "***************************************************\n";
      wf.log << "charge 1 & 2  "<< charge_1 << " "<<charge_2 <<"\n";
      wf.log << "charge nucleus  "<< charge_nucleus <<"\n";
      wf.log << "mas 1& 2  "<< mass_1 << " "<<mass_2<<"\n";
      wf.log << " ABV_A_x1 "<< ABV_A_x1 <<"\n";
      wf.log << " ABV_A_x2 "<< ABV_A_x2 <<"\n";
      wf.log << " ABV_A_x3 "<< ABV_A_x3 <<"\n";
      wf.log << " ABV_B_x1 "<< ABV_B_x1 <<"\n";
      wf.log << " ABV_B_x2 "<< ABV_B_x2 <<"\n";
      wf.log << " ABV_B_x3 "<< ABV_B_x3 <<"\n";
      wf.log << "softening_en " <<softening_en<<"\n";;
      wf.log << "***************************************************\n";


    };

    
    Potential_He_dresden(wavefunction &wf, double q1,double q2,double m1,double m2, double Q,double asqure_en) {
      charge_1 = q1;
      charge_2 = q2;
      charge_nucleus = Q;
      mass_1 = m1;
      mass_2 = m2;
 
      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
      
      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */

      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      

      updateABV(wf);

      wf.log << "***************************************************\n";
      wf.log << "Constructor Potential_He_dresden"<<"\n";
      wf.log << "***************************************************\n";
      wf.log << "charge 1 & 2  "<< charge_1 << " "<<charge_2 <<"\n";
      wf.log << "charge nucleus  "<< charge_nucleus <<"\n";
      wf.log << "mas 1& 2  "<< mass_1 << " "<<mass_2<<"\n";
      wf.log << " ABV_A_x1 "<< ABV_A_x1 <<"\n";
      wf.log << " ABV_A_x2 "<< ABV_A_x2 <<"\n";
      wf.log << " ABV_A_x3 "<< ABV_A_x3 <<"\n";
      wf.log << " ABV_B_x1 "<< ABV_B_x1 <<"\n";
      wf.log << " ABV_B_x2 "<< ABV_B_x2 <<"\n";
      wf.log << " ABV_B_x3 "<< ABV_B_x3 <<"\n";
      wf.log << "softening_en " <<softening_en<<"\n";;
      wf.log << "***************************************************\n";
    };
   

    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
     
      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
     	for( int j = 1 ; j <= wf.n2 ; j++ )
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      ABV_V[ wf.in3( k , j , i ) ] = charge_1*charge_2/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ] )
		+charge_1*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]/4.+( wf.x3[ k ]+wf.x2[ j ]/2. )*( wf.x3[ k ]+wf.x2[ j ]/2. )+softening_en )
		+charge_2*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]/4.+( wf.x3[ k ]-wf.x2[ j ]/2. )*( wf.x3[ k ]-wf.x2[ j ]/2. )+softening_en );
	      ABV_V[ wf.in3( k , j , i ) ] /= 3.;
	    }
    
    }
  };
  
 /**********************************************************************************************************************************/
  //Potential for helium in the Dresden model, see C Ruiz PRL 2006

  struct Potential_He_plus_dresden : ABVparam
  {

    /* Two electron in 3D the dresden model x1=rho x2=zr x3=zs, Helium problem */

    double charge_1;
    double charge_2;
    double charge_nucleus;
    double mass_1;
    double mass_2;
    double softening_en;

    /* Constructor to initialize the default values: */
    Potential_He_plus_dresden( wavefunction &wf) {
      charge_1=charge_2 = -1 ;
      charge_nucleus = 2;
      mass_1 = mass_2 = 1;
      softening_en = 0.135; //valid for dx=0.3 a.u.

      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
      
      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */

      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      
      updateABV(wf);

      cout << "***************************************************\n";
      cout << "Default Constructor Potential_He_dresden"<<"\n";
      cout << "***************************************************\n";
      cout << " charge 1 & 2   =  "<< charge_1 << " "<<charge_2 <<"\n";
      cout << " charge nucleus =  "<< charge_nucleus <<"\n";
      cout << " mas 1 & 2      =  "<< mass_1 << " "<<mass_2<<"\n";
      cout << " ABV_A_x1 =    "<< ABV_A_x1 <<"\n";
      cout << " ABV_A_x2 =    "<< ABV_A_x2 <<"\n";
      cout << " ABV_A_x3 =    "<< ABV_A_x3 <<"\n";
      cout << " ABV_B_x1*x1 = "<< ABV_A_x1 <<"\n";
      cout << " ABV_B_x2 =    "<< ABV_B_x2 <<"\n";
      cout << " ABV_B_x3 =    "<< ABV_B_x3 <<"\n";
      cout << "softening_en   " <<softening_en<<"\n";;
      cout << "***************************************************\n";


    };

    Potential_He_plus_dresden(wavefunction &wf, double q1,double q2,double m1,double m2, double Q,double asqure_en) {
      charge_1 = q1;
      charge_2 = q2;
      charge_nucleus = Q;
      mass_1 = m1;
      mass_2 = m2;
 
      double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
      
      double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
      double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */

      ABV_A_x1 =-0.5/reduced_mass;
      ABV_A_x2 =-0.5/reduced_mass;
      ABV_A_x3 =-0.5/center_of_mass;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }
      ABV_B_x2 = reduced_charge/reduced_mass;
      ABV_B_x3 = center_of_mass_charge/center_of_mass;
      

      updateABV(wf);
      
      cout << "***************************************************\n";
      cout << "Default Constructor Potential_He_dresden"<<"\n";
      cout << "***************************************************\n";
      cout << " charge 1 & 2   =  "<< charge_1 << " "<<charge_2 <<"\n";
      cout << " charge nucleus =  "<< charge_nucleus <<"\n";
      cout << " mas 1 & 2      =  "<< mass_1 << " "<<mass_2<<"\n";
      cout << " ABV_A_x1 =    "<< ABV_A_x1 <<"\n";
      cout << " ABV_A_x2 =    "<< ABV_A_x2 <<"\n";
      cout << " ABV_A_x3 =    "<< ABV_A_x3 <<"\n";
      cout << " ABV_B_x1*x1 = "<< ABV_A_x1 <<"\n";
      cout << " ABV_B_x2 =    "<< ABV_B_x2 <<"\n";
      cout << " ABV_B_x3 =    "<< ABV_B_x3 <<"\n";
      cout << "softening_en   " <<softening_en<<"\n";;
      cout << "***************************************************\n";

     
    };


    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
     
      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
     	for( int j = 1 ; j <= wf.n2 ; j++ )
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      ABV_V[ wf.in3( k , j , i ) ] =
		charge_1*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]/4.+( wf.x3[ k ]+wf.x2[ j ]/2. )*( wf.x3[ k ]+wf.x2[ j ]/2. )+softening_en )
		+charge_2*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]/4.+( wf.x3[ k ]-wf.x2[ j ]/2. )*( wf.x3[ k ]-wf.x2[ j ]/2. )+softening_en );
	      ABV_V[ wf.in3( k , j , i ) ] /= 3.;
	    }
    
    }
  }; //En Potential for helium in the Dresden model, see C Ruiz PRL 2006
  
 
  
   /*******************************************************************************************/
   /* Free potential, in Cylindrical 3D */

  struct Potential_Free : ABVparam
  {
    double mass_x1;
    double mass_x2;
    double mass_x3;
    
    double charge_x1;
    double charge_x2;
    double charge_x3;

    /* Constructor to initialize the default values:  Cylindrical mass 1 and charge1 for every coordinate */
    Potential_Free(const wavefunction &wf) {
      
      mass_x1=1;
      mass_x2=1;
      mass_x3=1;
      charge_x1 = -1; /* There is really no way to couple a laser to x1 in this namespace, anymay ... */
      charge_x2 = -1;
      charge_x3 = -1;

      ABV_A_x1 =-0.5/mass_x1;
      ABV_A_x2 =-0.5/mass_x2;
      ABV_A_x3 =-0.5/mass_x3;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }

      ABV_B_x2 = charge_x2/mass_x2;
      ABV_B_x3 = charge_x3/mass_x3;

      updateABV(wf);

       cout << "\n***************************************************\n";
      cout << "Constructor Potential_Free"<<"\n";
      cout << "***************************************************\n";
      cout << " charge 1 & 2 & 3 = "<< charge_x1 << " "<<charge_x2 << " "<<charge_x3 <<  "\n";
      cout << " mass   1 & 2 & 3 = "<< mass_x1 << " "<<mass_x2<< " "<<mass_x2<< "\n";
      cout << " ABV_A_x1   = "<< ABV_A_x1 <<"\n";
      cout << " ABV_A_x2   = "<< ABV_A_x2 <<"\n";
      cout << " ABV_A_x3   = "<< ABV_A_x3 <<"\n";
      cout << " ABV_B_x1*x1= "<< ABV_A_x1 <<"\n";
      cout << " ABV_B_x2   = "<< ABV_B_x2 <<"\n";
      cout << " ABV_B_x3   = "<< ABV_B_x3 <<"\n";
      cout << "***************************************************\n";


    }

    /* Constructor to change the values */
    Potential_Free(const wavefunction &wf, double mass1,double charge1, double mass2,double charge2, double mass3,double charge3 )
    {
      mass_x1=mass1;
      mass_x2=mass2;
      mass_x3=mass3;
      charge_x1 = charge1; /* There is really no way to couple a laser to x1 in this namespace, anymay ... */
      charge_x2 = charge2;
      charge_x3 = charge3;

      ABV_A_x1 =-0.5/mass_x1;
      ABV_A_x2 =-0.5/mass_x2;
      ABV_A_x3 =-0.5/mass_x3;

      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = ABV_A_x1/wf.x1[ i ];
        }

      ABV_B_x2 = charge_x2/mass_x2;
      ABV_B_x3 = charge_x3/mass_x3;
    
      updateABV(wf);

      cout << "\n***************************************************\n";
      cout << "Constructor Potential_Free"<<"\n";
      cout << "***************************************************\n";
      cout << " charge 1 & 2 & 3 = "<< charge_x1 << " "<<charge_x2 << " "<<charge_x3 <<  "\n";
      cout << " mass   1 & 2 & 3 = "<< mass_x1 << " "<<mass_x2<< " "<<mass_x2<< "\n";
      cout << " ABV_A_x1   = "<< ABV_A_x1 <<"\n";
      cout << " ABV_A_x2   = "<< ABV_A_x2 <<"\n";
      cout << " ABV_A_x3   = "<< ABV_A_x3 <<"\n";
      cout << " ABV_B_x1*x1= "<< ABV_A_x1 <<"\n";
      cout << " ABV_B_x2   = "<< ABV_B_x2 <<"\n";
      cout << " ABV_B_x3   = "<< ABV_B_x3 <<"\n";
      cout << "***************************************************\n";



    }
    
    void updateABV(const wavefunction &wf) {
     
      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 )*( wf.n3+2 ), 0.);

      /* following code is from code from Initialize_Potential */
      for( int k = 1 ; k <= wf.n3 ; k++ )
      {
	for( int j = 1 ; j <= wf.n2 ; j++ )
	  {
	    for( int i = 1 ; i <= wf.n1 ; i++ )
	      {
		ABV_V[ wf.in3( k , j , i ) ] = 0.;
	      }
	  }
      }
    }
  };

  /* End Free potential, in Cylindrical 3D */

} // end of Namespace Cylindrical3D




  namespace Cylindrical2D {
  
  struct ABVparam {
    
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x1;
    vector<double> ABV_B_x1;
    
    double ABV_A_x2;
    double ABV_B_x2;
    
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;
    vector<double> ABV_V_deriv;
    
  };
  
  struct Potential_H_2D : ABVparam
  {
    double charge;
    double charge_nucleus;
    double mass;
    double mass_nucleus;
    double softening_en;
    
    /* Constructor to initialize the default values: */
    Potential_H_2D(const wavefunction &wf) {
      /* what to do here ?? */
      charge = -1 ;
      charge_nucleus = 1.0;
      mass = 1;
      mass_nucleus = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      //softening_en = 0.000073;
      softening_en = 0.001;
      //softening_en = 0.0065;
      
      updateABV(wf);
    };

    Potential_H_2D(const wavefunction &wf, double _charge,double _charge_nucleus,double _mass , double _mass_nucleus, double _softening_en ) {
      /* what to do here ?? */
      charge = _charge ;
      charge_nucleus =_charge_nucleus;
      mass = _mass;
      mass_nucleus = _mass_nucleus;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = _softening_en ;
      
      updateABV(wf);
    };
    
    
    void updateABV(const wavefunction &wf) {
      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = -0.5/( mass*wf.x1[ i ]);
        }
      ABV_B_x2 = 0.;
      
      /* following code is from code from Initialize_grid */
      ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
      /* following code is from code from Initialize_Potential */
      for( int j = 1 ; j <= wf.n2 ; j++ )
	{
	  for( int i = 1 ; i <= wf.n1 ; i++ )
	    {
	      ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ] +softening_en )
		;
		ABV_V[ wf.in2( j , i ) ] /= 2.;
	    }
	}
    }
  };
  

  //user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential: ABVparam
  {
    double mass;
    double electron_charge;
    int n1, n2;
    
    user_input_potential(const wavefunction &wf, const vector<double> &potential_array, int _n1, int _n2) 
    {
      mass=1.0;
      electron_charge=1.0;   
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = -0.5/( mass*wf.x1[ i ]);
        }
      ABV_B_x2 = 0.;
      
      n1=_n1;
      n2=_n2;
 	          
      ABV_V.resize(( n1+2 )*( n2+2 ), 0.);
      
	for( int j = 1 ; j <= n2 ; j++ )
	  {
	    for( int i = 1 ; i <= n1 ; i++ )
	      {
		ABV_V[ j*(n1+2)+i ] = potential_array[(i-1)*n2+j-1];
		ABV_V[ j*(n1+2)+i ] /= 2.;		
		//cout<<ABV_V[ j*(n1+2)+i ]<<endl;
	      }
	  }
     };
	  
  }; 
  //Potential for H2+ for fixed internuclear distance R

  struct Potential_H2plus_fixed_R : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  double mass_nucleus;
	  double softening_en;
	  double internuclear_distance;
	  
	  /* Constructor to initialize the default values: */
	  Potential_H2plus_fixed_R(const wavefunction &wf, const double _R) {
		    
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;
		  
		  charge_nucleus=1;
		  mass_nucleus=1837;

		  
		  internuclear_distance=_R;
		  
		  updateABV(wf);
	  };

	  /* Constructor to input values: */
	  Potential_H2plus_fixed_R(const wavefunction &wf, const double charge_ion, const double soft_en, const double _R) {
		
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;

		  /* These parameters can be changed */
		  charge_nucleus=charge_ion;		
		  softening_en=1;
		  
		  internuclear_distance=_R;
		  
		  /* The rest stays the same */
		  mass_nucleus=1837;
		  
		  updateABV(wf);
	  };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = -0.5/( mass*wf.x1[ i ]);
        }
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt((wf.x2[ j ]+internuclear_distance/2.0)*(wf.x2[ j ]+internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[ i ])+charge*charge_nucleus/sqrt((wf.x2[ j ]-internuclear_distance/2.0)*(wf.x2[ j ]-internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[ i ])+charge_nucleus*charge_nucleus/internuclear_distance;	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };//End Potential for H2+ with a fixed R 


  //Potential for H2+ for fixed internuclear distance R with softcore

  struct Potential_H2plus_fixed_R_softcore : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  double mass_nucleus;
	  double softening_en;
	  double softening_nn;
	  double internuclear_distance;
	  
	  /* Constructor to initialize the default values: */
	  Potential_H2plus_fixed_R_softcore(const wavefunction &wf, const double _R) {
		    
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;
		  
		  charge_nucleus=1;
		  mass_nucleus=1837;
		  softening_en=0.0;
		  softening_nn=0.0;

		  
		  internuclear_distance=_R;
		  
		  updateABV(wf);
	  };

	  /* Constructor to input values: */
	  Potential_H2plus_fixed_R_softcore(const wavefunction &wf, const double charge_ion, const double soft_en, const double soft_nn, const double _R) {
		
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;

		  /* These parameters can be changed */
		  charge_nucleus=charge_ion;		
		  softening_en=soft_en;
		  softening_nn=soft_nn;
		  internuclear_distance=_R;
		  
		  /* The rest stays the same */
		  mass_nucleus=1837;
		  
		  updateABV(wf);
	  };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = -0.5/( mass*wf.x1[ i ]);
        }
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
			    ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt((wf.x2[ j ]+internuclear_distance/2.0)*(wf.x2[ j ]+internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[ i ] + softening_en)+charge*charge_nucleus/sqrt((wf.x2[ j ]-internuclear_distance/2.0)*(wf.x2[ j ]-internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[i] + softening_en);//+charge_nucleus*charge_nucleus/sqrt((internuclear_distance*internuclear_distance)+softening_nn);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };//End Potential for H2+ with a fixed R with softcore



 //Screwed up by Cory, 8/29/14

 /* struct Potential_H2plus_fixed_R_with_softcore_parameters: ABVparam
  {
        double charge  ;
        double charge_nucleus_1;
        double charge_nucleus_2;
	double charge_proton
        double screened_charge_factor;
        double mass;
        double mass_nucleus_1;
        double mass_nucleus_2;
	double mass_proton
        double softening_en;
        double softening_nn;
        double internuclear_distance;*/
	
	  
	  /* Constructor to initialize the default values: */
  /*	Potential_H2plus_fixed_R_with_softcore_parameters(const wavefunction &wf, const double _L, const double _Rr, const double _Rl, const double _Rf, const double _Rc) {
		    
		  /* Electron charge and mass */
	  /*		  charge=-1;
		  mass=1;
		  
		  mass_proton=1837;
		  charge_proton=1;

		  charge_nucleus_1=1;
		  charge_nucleus_2=1;
		  mass_nucleus_1=10837;
		  mass_nucleus_2=10837;
		  
		  internuclear_distance=_L;
		  
		  Rl=_Rl;
		  Rf=_Rf;
		  Rc=_Rc;
		    
		  updateABV(wf);
		  };*/

	  /* Constructor to input values: */
	  /*  Potential_H2plus_fixed_R_with_softcore_parameters(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
                               double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
                               double _softening_en,double _softening_nn,const double _R)
     {
       charge = _charge ;
       charge_nucleus_1=_charge_nucleus_1;
       charge_nucleus_2 = _charge_nucleus_2;
       screened_charge_factor=_screened_charge_factor ;
       mass = _mass;
       mass_nucleus_1 = _mass_nucleus_1;
       mass_nucleus_2 = _mass_nucleus_2;
       softening_en = _softening_en;
       softening_nn = _softening_nn;
       internuclear_distance=_R;
       updateABV(wf);
     };



	  void updateABV(const wavefunction &wf) {
	  */
      /* following code is from code from Create_Potential */
	  /*  ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1.resize( wf.n1+2 , 0);
      for( int i = 1 ; i <= wf.n1 ; i++ )
        {
          ABV_B_x1[ i ] = -0.5/( mass*wf.x1[ i ]);
        }
	ABV_B_x2 = 0.;*/	  
	          
		  
		  /* following code is from code from Initialize_grid */
	  /*  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);*/
		  
		  /* following code is from code from Initialize_Potential */
		  /*	  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+internuclear_distance/2.0)*(wf.x2[ j ]+internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[ i ]+softening_en)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-internuclear_distance/2.0)*(wf.x2[ j ]-internuclear_distance/2.0)+wf.x1[ i ]*wf.x1[ i ]+softening_en)+charge_nucleus_1*charge_nucleus_2/sqrt(internuclear_distance*internuclear_distance+softening_nn);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }*/
  //  };//End Potential for H2+ with a fixed R with soft-core parameters

  
  
} // end of Namespace Cylindrical2D


namespace Cartesian_2D {

  struct ABVparam {
  
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x1;
    double ABV_B_x1;
  
    double ABV_A_x2;
    double ABV_B_x2;
  
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;
    vector<double> ABV_V_deriv;

  };


  /***********************************************************************************************/
  /***********************************************************************************************/
  
  //Potential for Hydorgen like  2D

  struct Potential_H_2D : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  double mass_nucleus;
	  double softening_en;
	  
	  /* Constructor to initialize the default values: */
	  Potential_H_2D(const wavefunction &wf) {
		    
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;
		  
		  charge_nucleus=1;
		  mass_nucleus=1837;
		  //softening_en=1;
		  softening_en=0.0;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;
		  
		  updateABV(wf);
	  };

	  /* Constructor to input values: */
	  Potential_H_2D(const wavefunction &wf, const double charge_ion, const double soft_en) {
		
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;

		  /* These parameters can be changed */
		  charge_nucleus=charge_ion;		
		  softening_en=1;
		  
		  /* The rest stays the same */
		  mass_nucleus=1837;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;
		  
		  updateABV(wf);
	  };


	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt( wf.x1[ i ]*wf.x1[ i ]+wf.x2[ j ]*wf.x2[ j ] +softening_en );	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };//End Potential Hydrogen 2D

  /***********************************************************************************************/
  /***********************************************************************************************/


  //user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential: ABVparam
  {
    double mass;
    double electron_charge;
    int n1, n2;
    
    user_input_potential(const vector<double> &potential_array, int _n1, int _n2) 
    {
      mass=1.0;
      electron_charge=-1.0;   
      ABV_A_x1 =-0.5/mass;
      ABV_A_x2 =-0.5/mass;
      ABV_B_x1 = electron_charge/mass;
      ABV_B_x2 = electron_charge/mass;
      n1=_n1;
      n2=_n2;
 	          
      ABV_V.resize(( n1+2 )*( n2+2 ), 0.);
      
	for( int j = 1 ; j <= n2 ; j++ )
	  {
	    for( int i = 1 ; i <= n1 ; i++ )
	      {
		ABV_V[ j*(n1+2)+i ] = potential_array[(i-1)*n2+j-1];
		ABV_V[ j*(n1+2)+i ] /= 2.;		
	      }
	  }
     };
	  
  };
  
  //Potential for H2+ for fixed internuclear distance R

  struct Potential_H2plus_fixed_R : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  double mass_nucleus;
	  double softening_en;
	  double internuclear_distance;
	  
	  /* Constructor to initialize the default values: */
	  Potential_H2plus_fixed_R(const wavefunction &wf, const double _R) {
		    
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;
		  
		  charge_nucleus=1;
		  mass_nucleus=1837;

		  
		  internuclear_distance=_R;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = 0.0;
		  ABV_B_x2 = 0.0;
		  
		  updateABV(wf);
	  };

	  /* Constructor to input values: */
	  Potential_H2plus_fixed_R(const wavefunction &wf, const double charge_ion, const double soft_en, const double _R) {
		
		  /* Electron charge and mass */
		  charge=-1;
		  mass=1;

		  /* These parameters can be changed */
		  charge_nucleus=charge_ion;		
		  softening_en=1;
		  
		  internuclear_distance=_R;
		  
		  /* The rest stays the same */
		  mass_nucleus=1837;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = 0.0;
		  ABV_B_x2 = 0.0;
		  
		  updateABV(wf);
	  };


	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus/sqrt((wf.x1[ i ]+internuclear_distance/2.0)*(wf.x1[ i ]+internuclear_distance/2.0)+wf.x2[ j ]*wf.x2[ j ])+charge*charge_nucleus/sqrt((wf.x1[ i ]-internuclear_distance/2.0)*(wf.x1[ i ]-internuclear_distance/2.0)+wf.x2[ j ]*wf.x2[ j ])+charge_nucleus*charge_nucleus/internuclear_distance;	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };//End Potential for H2+ with a fixed R

  /***********************************************************************************************/
  /***********************************************************************************************/
  
//potential for H2+ using 2D model, created by Jing on 02/15/2010. here x1 is R, x2 is z. 
    struct Potential_H2_plus_e1D_R1D : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double softening_nn;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_R1D(const wavefunction &wf) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = 1.0;
      softening_en2=1.0;
      softening_nn = 0.03;
      
      updateABV(wf);
    };
    
    Potential_H2_plus_e1D_R1D(const wavefunction &wf, double _softening_en1, double _softening_en2, double _softening_nn  ) 
    {
      charge = -1 ;
      charge_nucleus_1=1;
      charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = 1836;
      mass_nucleus_2 = 1836;
      softening_en1 = _softening_en1;

      softening_en2 = _softening_en2;
      softening_nn = _softening_nn;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en2)+charge_nucleus_1*charge_nucleus_2/sqrt(wf.x1[i]*wf.x1[i]+softening_nn);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };
  
  
//Morse potential for H2+ using 2D model, created by Jing on 07/28/2010. here x1 is R, x2 is z. 
    struct Potential_H2_plus_e1D_R_Morse : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double softening_nn;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_R_Morse(const wavefunction &wf, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _re, double _mass) 
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = _mass;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;
      softening_nn = 0.03;
      De=_De;
      a=_a;
      Re=_re;
      
      updateABV(wf);
    };
    
    Potential_H2_plus_e1D_R_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn,double _De,double _a,double _Re  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en1 = _softening_en;
      softening_nn = _softening_nn;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(wf.x1[i]-Re)),2.0)-1.0);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };


    struct Potential_H2_plus_e1D_R_Morse_StaticField : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double softening_nn;
    double De;
    double a;
    double Re;
    double E_static;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_R_Morse_StaticField(const wavefunction &wf, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _re, double _mass, double _E_static) 
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = _mass;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;
      softening_nn = 0.03;
      De=_De;
      a=_a;
      Re=_re;
      E_static=_E_static;
      
      updateABV(wf);
    };
    
    Potential_H2_plus_e1D_R_Morse_StaticField(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn,double _De,double _a,double _Re, double _E_static) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en1 = _softening_en;
      softening_nn = _softening_nn;
      De=_De;
      a=_a;
      Re=_Re;
      E_static=_E_static;
      
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(wf.x1[i]-Re)),2.0)-1.0)+charge*E_static*wf.x2[ j ];	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };



//Morse potential for D2+ using 2D model, created by Jing on 11/30/2010. here x1 is R, x2 is z. 
    struct Potential_D2_plus_e1D_R_Morse : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_D2_plus_e1D_R_Morse(const wavefunction &wf, double Z1, double Z2, double _De, double _a, double _re) 
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836*4.0;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      De=_De;
      a=_a;
      Re=_re;
      
      updateABV(wf);
    };
    
    Potential_D2_plus_e1D_R_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn,double _De,double _a,double _Re  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en)+De*(pow(1.0-exp(-a*(wf.x1[i]-Re)),2.0)-1.0);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };
 
  
//Morse potential for H2+ using 2D model, created by Jing on 07/28/2010. here x1 is R, x2 is z. 
    struct Potential_H2_plus_e1D_R_Morse_new : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_R_Morse_new(const wavefunction &wf) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      De=0.2;
      a=1.2;
      Re=3.0;
      
      updateABV(wf);
    };
    
    Potential_H2_plus_e1D_R_Morse_new(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn,double _De,double _a,double _Re  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en)+De*(pow(1.0-exp(-a*(wf.x1[i]-Re)),2.0)-1.0);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };

  
//potential for D2+ using 2D model, created by Jing on 03/08/2010. here x1 is R, x2 is z. 
    struct Potential_D2_plus_e1D_R1D : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;

    /* Constructor to initialize the default values: */
    Potential_D2_plus_e1D_R1D(const wavefunction &wf) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 3671.5;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      
      updateABV(wf);
    };
    
    Potential_D2_plus_e1D_R1D(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en)+charge*charge_nucleus_2/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en)+charge_nucleus_1*charge_nucleus_2/sqrt(wf.x1[i]*wf.x1[i]+softening_nn);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };

 
  //Potential for Free electron in 2D
  struct Free_2D : ABVparam
  {
	  double charge;
	  double mass;
	  
	  /* Constructor to initialize the default values: */
	  Free_2D(const wavefunction &wf) {
		  
		  charge=-1;
		  mass=1;
		  updateABV(wf);

		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;

		  
	  };

	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =0.;
				  ABV_V[ wf.in2( j , i ) ] /= 2.;		// Just to remember . ..
			  }
		  }
	  }
  };
  
  /***********************************************************************************************/
  /***********************************************************************************************/


//Morse potential for N2+ using 2D model, created by Jing on 09/17/2010. here x1 is R, x2 is z. 
    struct Potential_N2_plus_e1D_R_Morse : ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_N2_plus_e1D_R_Morse(const wavefunction &wf, double _screened_charge_factor, double _De, double _a, double _re) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=_screened_charge_factor;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      De=_De;
      a=_a;
      Re=_re;
      
      updateABV(wf);
    };
    
    Potential_N2_plus_e1D_R_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn,double _De,double _a,double _Re  ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
    };


	  void updateABV(const wavefunction &wf) {

      /* following code is from code from Create_Potential */
      ABV_A_x1 =-0.5/(mass_nucleus_1*0.5);
      ABV_A_x2 =-0.5/mass;
      
      ABV_B_x1 = 0.;
      ABV_B_x2 = 0.;	  
	          
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] =  charge*charge_nucleus_1*screened_charge_factor/sqrt((wf.x2[ j ]+wf.x1[i]/2.0)*(wf.x2[ j ]+wf.x1[i]/2.0)+softening_en)+charge*charge_nucleus_2*screened_charge_factor/sqrt((wf.x2[ j ]-wf.x1[i]/2.0)*(wf.x2[ j ]-wf.x1[i]/2.0)+softening_en)+De*(pow(1.0-exp(-a*(wf.x1[i]-Re)),2.0)-1.0);	
      				  ABV_V[ wf.in2( j , i ) ] /= 2.;	
			  }
		  }
	  }
  };


  //Potential for Helium in 1D in Jacobi Coordinates
  struct Helium_2D_Jacobi : ABVparam
  {
	  double charge_1,charge_2;
	  double charge_nucleus;
	  double mass_1;
	  double mass_2;
	  double mass_nucleus;
	  double softening_en;
	  double softening_ee;
	  
	  /* Constructor to initialize the default values: */
	  Helium_2D_Jacobi(const wavefunction &wf) {
		  
		  charge_1=charge_2=-1.;
		  charge_nucleus=2;
		  mass_1=mass_2=1;
		 
		  softening_en=1.;
		  softening_ee=1.;
		  
		  double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
		  
		  double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */
		  
		  ABV_A_x1 =-0.5/center_of_mass;
		  ABV_A_x2 =-0.5/reduced_mass;
		  
		  ABV_B_x1 = center_of_mass_charge/center_of_mass;
		  ABV_B_x2 = reduced_charge/reduced_mass;	 
		  
		  updateABV(wf);
	  };

	  /* Constructor to initialize the default values: */
	  Helium_2D_Jacobi(const wavefunction &wf, const double q1, const double q2, const double m1, const double m2, 
			   const double soft_en, const double soft_ee, const double Q ) {
		  
		  charge_1=q1;
		  charge_2=q2;
		  charge_nucleus=Q;
		  mass_1=m1;
		  mass_2=m2;
	
		  softening_en=soft_en;
		  softening_ee=soft_ee;
		  
		  double reduced_mass =  mass_1*mass_2/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass = mass_1 + mass_2;                  /* Helper variables */
		  
		  double reduced_charge =  (mass_2*charge_1-mass_1*charge_2)/( mass_1 + mass_2 ); /* Helper variables */
		  double center_of_mass_charge =  charge_1+charge_2;                              /* Helper variables */
		  
		  ABV_A_x1 =-0.5/center_of_mass;
		  ABV_A_x2 =-0.5/reduced_mass;
		  
		  ABV_B_x1 = center_of_mass_charge/center_of_mass;
		  ABV_B_x2 = reduced_charge/reduced_mass;	 
		  
		  updateABV(wf);
	  };

	  void updateABV(const wavefunction &wf) {
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  for( int j = 1 ; j <= wf.n2 ; j++ )
		  {
			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( j , i ) ] = 
					  charge_1*charge_2/sqrt( wf.x2[ j ]*wf.x2[ j ]+softening_ee )
					  +charge_1*charge_nucleus/sqrt( ( wf.x1[ i ]+wf.x2[ j ]/2. )*( wf.x1[ i ]+wf.x2[ j ]/2. )+softening_en )
					  +charge_2*charge_nucleus/sqrt( ( wf.x1[ i ]-wf.x2[ j ]/2. )*( wf.x1[ i ]-wf.x2[ j ]/2. )+softening_en );
				  
				  ABV_V[ wf.in2( j , i ) ] /= 2.;		
			  }
		  }
	  }
  };
  
  /***********************************************************************************************/
  /***********************************************************************************************/

  //Potential for Circular molecules in 2D
  struct Molecule_2D_Circular : ABVparam
  {
	  
	  vector<double> charge_nucleus;
	  vector<double> mass_nucleus;
	  
	  double charge;
	  double mass;
	  double softening_en;
	  double radius;
	  int  number_ions;
	  /* Constructor to initialize the default values: */
	 Molecule_2D_Circular  (const wavefunction &wf, const int no_centers, double radii) {
		 
		 number_ions=no_centers;
                 charge_nucleus.resize( number_ions+2, 1.);
                 mass_nucleus.resize( number_ions+2, 1837.);
		 
		 charge=-1.;
		 mass=1;
		 softening_en=1.;
		 
		 radius=radii;
		 
		 ABV_A_x1 =-0.5/mass;
		 ABV_A_x2 =-0.5/mass;
		 
		 ABV_B_x1 = charge/mass;
		 ABV_B_x2 = charge/mass;	 
		 
		 updateABV(wf,radius);
	 };
	  
	 void updateABV(const wavefunction &wf,const double radii) {
		 
		 /* following code is from code from Initialize_grid */
		 ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		 
		 /* following code is from code from Initialize_Potential */
		 
		 radius=radii;
		 for(int ion=1;ion<=number_ions;ion++)
		 {
			 double coord_x=radius*cos(ion*2.*pi/number_ions);
			 double coord_y=radius*sin(ion*2.*pi/number_ions);
			 
			 for( int j = 1 ; j <= wf.n2 ; j++ )
			 {
				 for( int i = 1 ; i <= wf.n1 ; i++ )
				 {
					 ABV_V[ wf.in2( j , i ) ] +=
						 charge*charge_nucleus[ion]/sqrt( 
							 (wf.x1[ i ]-coord_x )*(wf.x1[ i ]-coord_x )+
							 (wf.x2[ j ]-coord_y )*(wf.x2[ j ]-coord_y )+ softening_en );
					 
				 }
			 }
		 }
		 
		 for( int j = 1 ; j <= wf.n2 ; j++ )
			 for( int i = 1 ; i <= wf.n1 ; i++ )
				 ABV_V[ wf.in2( j , i ) ] /= 2.;		
		 
		 
	 }
  };

  /***********************************************************************************************/
  /***********************************************************************************************/
  //Potential for Polyatomic Molecules in 2D

  struct Molecule_2D: ABVparam
  {
	  
	  vector<double> charge_nucleus;
	  vector<double> mass_nucleus;
	  vector<double> x_coord_nucleus;
	  vector<double> y_coord_nucleus;
	  
	  double charge;
	  double mass;
	  double softening_en;
	  int  number_ions;
	  
	 /* Constructor to initialize the default values: */
	  Molecule_2D(const wavefunction &wf, const int no_centers) {
		  
		  number_ions=no_centers;
		  charge_nucleus.resize(  number_ions+2, -1.);
		  mass_nucleus.resize(    number_ions+2, 1837.);
		  x_coord_nucleus.resize( number_ions+2, 0.);
		  y_coord_nucleus.resize( number_ions+2, 0.);

		  //BY DEFAULT we make a circular colecule, but there is a more general constructor
		  double radius=4.;
		  for(int ion=1;ion<=number_ions;ion++)
		  {
			  x_coord_nucleus[ion]=radius*cos(ion*2.*pi/number_ions);
			  y_coord_nucleus[ion]=radius*sin(ion*2.*pi/number_ions);
		  }
		  
		  /*
		    Parameter for the elctron
		  */
		  
		  charge=-1.;
		  mass=1;
		  softening_en=1.;
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;	 
		  
		  updateABV(wf);
	  };
	  
	  Molecule_2D(const wavefunction &wf, const int no_centers, 
		      const vector<double> charges, const vector<double> masses,
		      const vector<double> x_position, const vector<double> y_position ) {
		  
		  number_ions=no_centers;
		  charge_nucleus.resize(  number_ions+2, -1.);
		  mass_nucleus.resize(    number_ions+2, 1837.);
		  x_coord_nucleus.resize( number_ions+2, 0.);
		  y_coord_nucleus.resize( number_ions+2, 0.);
		  
		  if(charges.size()!=no_centers)
		  {
			  cerr<<"vector charges has wrong size for the potential\n";
			  exit(1);
		  }
		  if(masses.size()!=no_centers)
		  {
			  cerr<<"vector masses has wrong size for the potential\n";
			  exit(1);
		  }
		  if(x_position.size()!=no_centers)
		  {
			  cerr<<"vector x position has wrong size for the potential\n";
			  exit(1);
		  }
		  if(y_position.size()!=no_centers)
		  {
			  cerr<<"vector y position has wrong size for the potential\n";
			  exit(1);
		  }
		  
		  
		  
		 for(int ion=1;ion<=number_ions;ion++)
		 {
			 x_coord_nucleus[ion]=x_position[ion-1];
			 y_coord_nucleus[ion]=y_position[ion-1];
			 charge_nucleus[ion]= charges[ion-1];
			 mass_nucleus[ion]= masses[ion-1];
		 }
		 
		 /*
		   Parameter for the elctron
		 */
		 
		 charge=-1.;
		 mass=1;
		 softening_en=1.;
		 
		 ABV_A_x1 =-0.5/mass;
		 ABV_A_x2 =-0.5/mass;
		 
		 ABV_B_x1 = charge/mass;
		 ABV_B_x2 = charge/mass;	 
		 
		 updateABV(wf);
	  };
	  
	  
	  void updateABV(const wavefunction &wf) {
		  
		  
		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		  
		  /* following code is from code from Initialize_Potential */
		  
		  for(int ion=1;ion<=number_ions;ion++)
		  {
			  
			  for( int j = 1 ; j <= wf.n2 ; j++ )
			  {
				  for( int i = 1 ; i <= wf.n1 ; i++ )
				  {
					  ABV_V[ wf.in2( j , i ) ] +=
						  charge*charge_nucleus[ion]/sqrt( 
							  (wf.x1[ i ]-x_coord_nucleus[ion] )*(wf.x1[ i ]-x_coord_nucleus[ion] )+
							  (wf.x2[ j ]-y_coord_nucleus[ion] )*(wf.x2[ j ]-y_coord_nucleus[ion] )+ softening_en );
					  
				  }
			  }
		  }
		  
		  for( int j = 1 ; j <= wf.n2 ; j++ )
			  for( int i = 1 ; i <= wf.n1 ; i++ )
				  ABV_V[ wf.in2( j , i ) ] /= 2.;		
		  
		  
	  }
	  
	  //This one updates the molecule very fast, just update the position, but not re-allocate the memory
	  void updateABV_fast(const wavefunction &wf, const vector<double> x_position, const vector<double> y_position  ) {
		  
		  if(x_position.size()!=number_ions)
		  {
			  cerr<<"vector x position has wrong size for the potential\n";
			  exit(1);
		  }
		  if(y_position.size()!=number_ions)
		  {
			  cerr<<"vector y position has wrong size for the potential\n";
			  exit(1);
		  }

		  for(int ion=1;ion<=number_ions;ion++)
		  {
			  x_coord_nucleus[ion]=x_position[ion-1];
			  y_coord_nucleus[ion]=y_position[ion-1];
			  
		  }
		  /* following code is from code from Initialize_Potential */
		  
		  for(int ion=1;ion<=number_ions;ion++)
		  {
			  for( int j = 1 ; j <= wf.n2 ; j++ )
			  {
				  for( int i = 1 ; i <= wf.n1 ; i++ )
				  {
					  ABV_V[ wf.in2( j , i ) ] +=
						  charge*charge_nucleus[ion]/sqrt( 
							  (wf.x1[ i ]-x_coord_nucleus[ion] )*(wf.x1[ i ]-x_coord_nucleus[ion] )+
							  (wf.x2[ j ]-y_coord_nucleus[ion] )*(wf.x2[ j ]-y_coord_nucleus[ion] )+ softening_en );
				  }
			  }
		  }
		  
		  for( int j = 1 ; j <= wf.n2 ; j++ )
			  for( int i = 1 ; i <= wf.n1 ; i++ )
				  ABV_V[ wf.in2( j , i ) ] /= 2.;		
		  
		  
	  }
  };


  /***********************************************************************************************/
  /***********************************************************************************************/
  //Potential for Benzenec Molecules in 2D
  
  struct Benzene_2D : ABVparam
  {
	  double charge;
	  double charge_nucleus;
	  double mass;
	  
	  double mass_nucleus;
	  double softening_en;
	  double softening_ee;
	  double radius;
	  
	  /* Constructor to initialize the default values: */
	  Benzene_2D(const wavefunction &wf) {
		 		  
		  charge=-1.;
		  charge_nucleus=1.;
		  mass=1;
		  mass_nucleus=1837;
		  softening_en=1.;
		  softening_ee=1.;
		  radius=1.98488; //See wikipedia
		  
		  
		  ABV_A_x1 =-0.5/mass;
		  ABV_A_x2 =-0.5/mass;
		  
		  ABV_B_x1 = charge/mass;
		  ABV_B_x2 = charge/mass;	 
		  
		  updateABV(wf);
	  };
	  
	  void updateABV(const wavefunction &wf) {
		 		  
		  /* following code is from code from Initialize_grid */
		  ABV_V.resize(( wf.n1+2 )*( wf.n2+2 ), 0.);
		 
		 /* following code is from code from Initialize_Potential */

		 double n1_x1=radius*cos(1.*2.*pi/6.);
		 double n1_x2=radius*cos(2.*2.*pi/6.);
		 double n1_x3=radius*cos(3.*2.*pi/6.);
		 double n1_x4=radius*cos(4.*2.*pi/6.);
		 double n1_x5=radius*cos(5.*2.*pi/6.);
		 double n1_x6=radius*cos(6.*2.*pi/6.);

		 double n1_y1=radius*sin(1.*2.*pi/6.);
		 double n1_y2=radius*sin(2.*2.*pi/6.);
		 double n1_y3=radius*sin(3.*2.*pi/6.);
		 double n1_y4=radius*sin(4.*2.*pi/6.);
		 double n1_y5=radius*sin(5.*2.*pi/6.);
		 double n1_y6=radius*sin(6.*2.*pi/6.);

		 for( int j = 1 ; j <= wf.n2 ; j++ )
		 {
			 for( int i = 1 ; i <= wf.n1 ; i++ )
			 {
				 ABV_V[ wf.in2( j , i ) ] =
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x1 )*(wf.x1[ i ]-n1_x1 )+(wf.x2[ j ]-n1_y1 )*(wf.x2[ j ]-n1_y1 )+ softening_en )+
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x2 )*(wf.x1[ i ]-n1_x2 )+(wf.x2[ j ]-n1_y2 )*(wf.x2[ j ]-n1_y2 )+ softening_en )+
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x3 )*(wf.x1[ i ]-n1_x3 )+(wf.x2[ j ]-n1_y3 )*(wf.x2[ j ]-n1_y3 )+ softening_en )+
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x4 )*(wf.x1[ i ]-n1_x4 )+(wf.x2[ j ]-n1_y4 )*(wf.x2[ j ]-n1_y4 )+ softening_en )+
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x5 )*(wf.x1[ i ]-n1_x5 )+(wf.x2[ j ]-n1_y5 )*(wf.x2[ j ]-n1_y5 )+ softening_en )+
					 charge*charge_nucleus/sqrt( (wf.x1[ i ]-n1_x6 )*(wf.x1[ i ]-n1_x6 )+(wf.x2[ j ]-n1_y6 )*(wf.x2[ j ]-n1_y6 )+ softening_en );
	
				

				 ABV_V[ wf.in2( j , i ) ] /= 2.;		// wrong! devided by two!!
			 }
		 }
	 }
  };

 
 /***********************************************************************************************/
  /***********************************************************************************************/
 





  
} // end of Namespace Cartesian2D


namespace Spherical_1D {

  struct ABVparam {
  
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x;
    double ABV_B_x;
  
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;

  };
  
    struct user_input_potential: ABVparam
  {
    double mass;
    
    user_input_potential(const vector<double> &potential_array) 
    {

      mass=1.0;   
      ABV_A_x =-0.5/mass;     
      ABV_B_x = 0.0;
 	          
      ABV_V.resize(( potential_array.size()+2 )*( 0+2 ), 0.0);
	
      //cout<<"as a check"<<endl;	  
      for( int i = 1 ; i <= potential_array.size() ; i++ )
      {
      ABV_V[ (potential_array.size()+2)*1+i  ] =  potential_array[i-1];
      //cout<<ABV_V[ (potential_array.size()+2)*1+i  ]<<endl;	
      }
      //cout<<endl;
     };
	  
  };
  
}


namespace Cartesian_1D {

  struct ABVparam {
  
    /**
     * Declare ABV coefficients.  
     */
    double ABV_A_x;
    double ABV_B_x;
  
    /**
     * Declare potentials.  
     */
    vector<double> ABV_V;

  };

  
//potential for H2+ using 2D model, created by Jing on 02/15/2009. here x1 is z, R is fixed. 
    struct Potential_H2_plus_e1D : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double softening_nn;
    double R;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = 1.0;

      softening_en2 = 1.0;
      softening_nn = 0.03;
      R=_R;
      
      updateABV(wf);
    };

    Potential_H2_plus_e1D(const wavefunction &wf, double _R, double _softening_en1, double _softening_en2, double _softening_nn)
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;

      softening_en2 = _softening_en2;

      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en2)+charge_nucleus_1*charge_nucleus_2/sqrt(R*R+softening_nn);		
			  }
		  }
	  
  };

 //BO Potential for Metiu system AHB2+ w/ A and B fixed, solve for H+/"moving nucleus" and 1e-
 //Created by Cory, 8/29/2014
 struct AHB2BO : ABVparam
  {
	double charge;
        double charge_nucleus_1;
        double charge_nucleus_2;
	double charge_proton;
        double mass;
        double mass_nucleus_1;
        double mass_nucleus_2;
	double mass_proton;
	double Rl;
	double Rf;
	double Rr;
	double Rp;
	double L;
	  /* Constructor to initialize the default values: */
	AHB2BO(const wavefunction &wf, const double _L, const double _Rr, const double _Rl, const double _Rf, const double _Rp) {   
	  /* Electron charge and mass */
	    charge=-1;
	    mass=1;
	    /*Charge and mass of moving proton*/  
	    mass_proton=1837;
	    charge_proton=1;
	    /*Charge and mass of fixed ions, not really needed
	    charge_nucleus_1=1;
	    charge_nucleus_2=1;
	    mass_nucleus_1=183700;
	    mass_nucleus_2=183700;
	    */
	    L=_L;
	    Rp=_Rp; //moving nuclear BO fixed distance
	    Rl=_Rl;
	    Rf=_Rf;
	    Rr=_Rr;
	    
	    updateABV(wf);
	  };

	  /* Constructor to input values: */
	/*AHB2BO(const wavefunction &wf, double _L, double _Rr, double _Rl, double _Rf, double _R){
	    charge=-1;
	    mass=1;
	    
	    mass_proton=1837;
	    charge_proton=1;
	    
	   
	    R=_R; //moving nuclear BO fixed distance
	    Rl=_Rl;
	    Rf=_Rf;	  			 
	    Rr=_Rr;
	    updateABV(wf);
	    };*/

	void updateABV(const wavefunction &wf) {
	  
	  
	  ABV_A_x =-0.5/mass;
	  
	  ABV_B_x = 0.;
	  
	  
	  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
	  
	  
	  for( int i = 1 ; i <= wf.n1 ; i++ ){
	    ABV_V[ wf.in2( 1 , i ) ] =  (-erf(abs(Rp-wf.x1[i])/Rf)/(abs(Rp-wf.x1[i]))) -  (erf(abs(wf.x1[i]-L/2)/Rr)/abs(wf.x1[i]-L/2)) -  (erf(abs(wf.x1[i]+L/2)/Rl)/abs(wf.x1[i]+L/2)) + (1/abs(L/2-Rp)) + (1/abs(L/2+Rp)); 
	  }
	}
  };//End Metiu/Gross potential

//potential for H2+ using 2D model, created by Jing on 02/15/2009. here x1 is z, R is fixed. 
    struct Potential_HD_plus_e1D : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double softening_nn;
    double R;

    /* Constructor to initialize the default values: */
    Potential_HD_plus_e1D(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = 1836;
      mass_nucleus_2 = 1836*2.0;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = 1.0;

      softening_en2 = 1.0;
      softening_nn = 0.03;
      R=_R;
      
      updateABV(wf);
    };

    Potential_HD_plus_e1D(const wavefunction &wf, double _R, double _softening_en1, double _softening_en2, double _softening_nn)
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;
      mass = 1;
      mass_nucleus_1 = 1836;
      mass_nucleus_2 = 1836*2.0;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;

      softening_en2 = _softening_en2;

      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]-R/3.0)*(wf.x1[ i ]-R/3.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]+R*2.0/3.0)*(wf.x1[ i ]+R*2.0/3.0)+softening_en2)+charge_nucleus_1*charge_nucleus_2/sqrt(R*R+softening_nn);		
			  }
		  }
	  
  };


//morse potential for H2+ using 2D model, created by Jing on 07/28/2010. here x1 is z, R is fixed. 
    struct Potential_H2_plus_e1D_Morse : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    //double softening_nn;
    double R;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_Morse(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = softening_en2=1.0;
      //softening_nn = 0.03;
      R=_R;
      De=0.2;
      Re=3.0;
      a=0.5;
      
      updateABV(wf);
    };

    Potential_H2_plus_e1D_Morse(const wavefunction &wf, double _R, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _Re)
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      screened_charge_factor=1;
      mass = 1;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;

      //softening_nn = _softening_nn;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(R-Re)),2.0)-1.0);
			  }
		  }
	  
  };


//morse potential for H2+ using 2D model, created by Jing on 07/28/2010. here x1 is z, R is fixed. 
    struct Potential_H2_plus_e1D_Morse_new : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double R;
    double De;
    double a;
    double Re, epsilon;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_Morse_new(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 0.5;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      //softening_nn = 0.03;
      R=_R;
      De=0.2;
      Re=3.0;
      a=1.2;
      epsilon=0.0;
      
      updateABV(wf);
    };

    Potential_H2_plus_e1D_Morse_new(const wavefunction &wf, double _R, double _softening_en, double _De, double _a, double _Re, double _epsilon)
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = _softening_en;
      //softening_nn = _softening_nn;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      epsilon=_epsilon;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en)+De*(pow(1.0-exp(-a*(R-Re*(1+epsilon*wf.x1[i]))),2.0)-1.0);
			  }
		  }
	  
  };


  
//potential for revised H2+ using 2D model, created by Jing on 07/28/2010. here x1 is z, R is fixed. 
    struct Potential_H2_plus_e1D_revised : ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double R;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_revised(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      R=_R;
      
      updateABV(wf);
    };

    Potential_H2_plus_e1D_revised(const wavefunction &wf, double _R, double q, double Q, double a, double b)
    {
      charge = q ;
      charge_nucleus_1=charge_nucleus_2 = Q;
      screened_charge_factor=1;
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = a;
      softening_nn = b;
      R=_R;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_revised(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en)+charge_nucleus_1*charge_nucleus_2/sqrt(R*R+softening_nn);		
			  }
		  }
	  
  };
  
  
//user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential_M1 : ABVparam
  {
    double mass;
    
    user_input_potential_M1(const vector<double> &potential_array) 
    {

      //mass=1836.0/2.0;
      mass=1.0;   
      ABV_A_x =-0.5/mass;     
      ABV_B_x = 0.0;
 	          
      ABV_V.resize(( potential_array.size()+2 )*( 0+2 ), 0.0);
	
      //cout<<"as a check"<<endl;	  
      for( int i = 1 ; i <= potential_array.size() ; i++ )
      {
      ABV_V[ (potential_array.size()+2)*1+i  ] =  potential_array[i-1];
      //cout<<ABV_V[ (potential_array.size()+2)*1+i  ]<<endl;	
      }
      //cout<<endl;
     };
	  
  };
  
  
  //user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential_M4: ABVparam
  {
    double mass;
    
    user_input_potential_M4(const vector<double> &potential_array) 
    {

      mass=1836.0*4.0/2.0;   
      ABV_A_x =-0.5/mass;     
      ABV_B_x = 0.0;
 	          
      ABV_V.resize(( potential_array.size()+2 )*( 0+2 ), 0.0);
	
      //cout<<"as a check"<<endl;	  
      for( int i = 1 ; i <= potential_array.size() ; i++ )
      {
      ABV_V[ (potential_array.size()+2)*1+i  ] =  potential_array[i-1];
      //cout<<ABV_V[ (potential_array.size()+2)*1+i  ]<<endl;	
      }
      //cout<<endl;
     };
	  
  };
  
  
  //user defined input potential, created by Jing on 02/22/2010
    struct user_input_potential: ABVparam
  {
    double mass;
    
    user_input_potential(const vector<double> &potential_array) 
    {

      mass=1.0;   
      ABV_A_x =-0.5/mass;     
      ABV_B_x = 0.0;
 	          
      ABV_V.resize(( potential_array.size()+2 )*( 0+2 ), 0.0);
	
      //cout<<"as a check"<<endl;	  
      for( int i = 1 ; i <= potential_array.size() ; i++ )
      {
      ABV_V[ (potential_array.size()+2)*1+i  ] =  potential_array[i-1];
      //cout<<ABV_V[ (potential_array.size()+2)*1+i  ]<<endl;	
      }
      //cout<<endl;
     };
	  
  };
  
  
//morse potential for N2+ using 2D model, created by Jing on 09/17/2010. here x1 is z, R is fixed. 
    struct Potential_N2_plus_e1D_Morse : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en;
    double softening_nn;
    double R;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_N2_plus_e1D_Morse(const wavefunction &wf, double _screened_charge_factor, double _R, double _De, double _a, double _Re) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=_screened_charge_factor;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = 1.0;
      softening_nn = 0.03;
      R=_R;
      De=_De;
      Re=_a;
      a=_Re;
      
      updateABV(wf);
    };

    Potential_N2_plus_e1D_Morse(const wavefunction &wf, double _R, double _screened_charge_factor, double _softening_en, double _softening_nn, double _De, double _a, double _Re)
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=_screened_charge_factor;
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1*screened_charge_factor/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en)+charge*charge_nucleus_2*screened_charge_factor/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en)+De*(pow(1.0-exp(-a*(R-Re)),2.0)-1.0);
			  }
		  }
	  
  };
 

//coulomb potential for N2+ using 2D model, created by Jing on 11/17/2010. here x1 is z, R is fixed. 
    struct Potential_N2_plus_e1D_Coulomb : ABVparam
  {
    double e_charge;
    double Z1_eff;
    double mass;
    double a_ec;
    double R;

    /* Constructor to initialize the default values: */
    Potential_N2_plus_e1D_Coulomb(const wavefunction &wf, double _Z1_eff, double _a_ec, double _R) 
    {
      e_charge=-1.0;
      Z1_eff=_Z1_eff;     
      mass=1.0;
      a_ec=_a_ec;
      R=_R;
      
      updateABV(wf);
    };

	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  e_charge*Z1_eff/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+a_ec*a_ec)+e_charge*Z1_eff/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+a_ec*a_ec);
			  }
		  }
	  
  };
  
  
  
    struct Potential_HD_plus_e1D_Morse : ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double R;
    double De;
    double a;
    double Re;

    /* Constructor to initialize the default values: */
    Potential_HD_plus_e1D_Morse(const wavefunction &wf, double _R) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;       
      mass = 1;
      mass_nucleus_1 = 1836;
      mass_nucleus_2 = 1836*2;
      softening_en1 = softening_en2=1.0;
      R=_R;
      De=0.2;
      Re=3.0;
      a=0.5;
      
      updateABV(wf);
    };

    Potential_HD_plus_e1D_Morse(const wavefunction &wf, double _R, double _mass1, double _mass2, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _Re)
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      mass = 1;
      mass_nucleus_1 = _mass1;
      mass_nucleus_2 = _mass2;
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(R-Re)),2.0)-1.0);
			  }
		  }
	  
  };
  
  
  
    struct Potential_HD_plus_e1D_Morse_StaticField : ABVparam
  {
    double charge;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    double R;
    double De;
    double a;
    double Re;
    double E_static;

    /* Constructor to initialize the default values: */
    Potential_HD_plus_e1D_Morse_StaticField(const wavefunction &wf, double _R, double _E_static) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;       
      mass = 1;
      mass_nucleus_1 = 1836;
      mass_nucleus_2 = 1836*2;
      softening_en1 = softening_en2=1.0;
      R=_R;
      De=0.2;
      Re=3.0;
      a=0.5;
      E_static=_E_static;
      
      updateABV(wf);
    };

    Potential_HD_plus_e1D_Morse_StaticField(const wavefunction &wf, double _R, double _mass1, double _mass2, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _Re, double _E_static)
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      mass = 1;
      mass_nucleus_1 = _mass1;
      mass_nucleus_2 = _mass2;
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      E_static=_E_static;
      
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(R-Re)),2.0)-1.0)+charge*E_static*wf.x1[ i ];
			  }
		  }
	  
  };



    struct Potential_H2_plus_e1D_Morse_StaticField : ABVparam
  {
    double charge  ;
    double charge_nucleus_1;
    double charge_nucleus_2;
    double screened_charge_factor;
    double mass;
    double mass_nucleus_1;
    double mass_nucleus_2;
    double softening_en1, softening_en2;
    //double softening_nn;
    double R;
    double De;
    double a;
    double Re;
    double E_static;

    /* Constructor to initialize the default values: */
    Potential_H2_plus_e1D_Morse_StaticField(const wavefunction &wf, double _R, double _E_static) 
    {
      charge = -1 ;
      charge_nucleus_1=charge_nucleus_2 = 1;
      screened_charge_factor=1;         
      mass = 1;
      mass_nucleus_1 = mass_nucleus_2 = 1836;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = softening_en2=1.0;
      //softening_nn = 0.03;
      R=_R;
      De=0.2;
      Re=3.0;
      a=0.5;
      E_static=_E_static;
      
      updateABV(wf);
    };

    Potential_H2_plus_e1D_Morse_StaticField(const wavefunction &wf, double _R, double Z1, double Z2, double _softening_en1, double _softening_en2, double _De, double _a, double _Re, double _E_static)
    {
      charge = -1 ;
      charge_nucleus_1 = Z1;
      charge_nucleus_2 = Z2;
      screened_charge_factor=1;
      mass = 1;
      /* Softening Parameter for the electron-nuclei interaction */
      softening_en1 = _softening_en1;
      softening_en2 = _softening_en2;

      //softening_nn = _softening_nn;
      R=_R;
      De=_De;
      a=_a;
      Re=_Re;
      E_static=_E_static;
      
      updateABV(wf);
   };

    /*
    Potential_H2_plus_e1D_Morse(const wavefunction &wf,double _charge,double _charge_nucleus_1,double _charge_nucleus_2,
			      double _screened_charge_factor,double _mass,double _mass_nucleus_1,double _mass_nucleus_2,
			      double _softening_en,double _softening_nn, double _R ) 
    {
      charge = _charge ;
      charge_nucleus_1=_charge_nucleus_1;
      charge_nucleus_2 = _charge_nucleus_2;
      screened_charge_factor=_screened_charge_factor ;         
      mass = _mass;
      mass_nucleus_1 = _mass_nucleus_1;
      mass_nucleus_2 = _mass_nucleus_2;
      softening_en = _softening_en;
      softening_nn = _softening_nn;
      R=_R;
      updateABV(wf);
    };
    */


	  void updateABV(const wavefunction &wf) {


      ABV_A_x =-0.5/mass;
      
      ABV_B_x = 0.;
 
	          
		  ABV_V.resize(( wf.n1+2 )*( 0+2 ), 0.);
		  

			  for( int i = 1 ; i <= wf.n1 ; i++ )
			  {
				  ABV_V[ wf.in2( 1 , i ) ] =  charge*charge_nucleus_1/sqrt((wf.x1[ i ]+R/2.0)*(wf.x1[ i ]+R/2.0)+softening_en1)+charge*charge_nucleus_2/sqrt((wf.x1[ i ]-R/2.0)*(wf.x1[ i ]-R/2.0)+softening_en2)+De*(pow(1.0-exp(-a*(R-Re)),2.0)-1.0)-charge*E_static*wf.x1[ i ];
			  }
		  }
	  
  };

  
} // end of Namespace Cartesian1D








#endif	/* POTENTIALS_H */


/*main() {
 
  wavefunction ground;
  Potential_H2 pot(ground);
  internuclear_distance = 4;

  Hamiltonian H(ground);
  loop {
    H.X1(dt, ground, pot);
  }
}
*/
