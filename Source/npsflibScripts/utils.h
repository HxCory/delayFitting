#ifndef UTILS_H
#include <cstdlib>
#define UTILS_H
#include <math.h>
#include <stdio.h>
#include <vector>
#include <complex>
//#include <fstream>
//#include <iostream>
// shaohao2009.04
#include "constant.h"   // for value of pi
#include <fftw3.h>


using namespace std;

inline void nrerror ( char error_text[] )
{
    fprintf( stderr, "numerical recipes runtime error..\n" );
    fprintf( stderr, "%s\n", error_text );
    fprintf( stderr,"... now exiting system..\n" );
}

/**
 *  The tridag routine to solve the threediagonal system of linear equation in Eq. (6) in doc/ABV/abv.tex
 */
inline void Tridag( const vector<complex<double> > &tridag_lower_diagonal , const vector<complex<double> > &tridag_middle_diagonal , const vector<complex<double> >&tridag_upper_diagonal, const vector<complex<double> > &wf_1D_rightside ,
	     vector<complex<double> > &wf_1D_tmp2 , vector<complex<double> > &gam )
{   
  const unsigned int length=gam.size() -2 ;

/**
 *  forth substitution
 */

  complex<double>  bet=tridag_middle_diagonal[ 1 ];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length; ++index )
    {
      gam[ index ] = tridag_upper_diagonal[ index-1 ]/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal[ index ]*gam[ index ];
      /*      if ( bet == 0.0 )
	{
	 cout<<"Error: bet==0 in tridag! at "<<index<<"\n";
         cout<<"tridag_middle_diagonal="<<tridag_middle_diagonal[ index ]<<"\n";
         cout<<"tridag_lower_diagonal="<<tridag_lower_diagonal[ index ]<<"\n";
         cout<<"gam="<<gam[ index ]<<"\n";
	 } */

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal[ index ]*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
/**
 *  back substitution
 */
  for( unsigned int index = length-1 ; index >= 1 ; --index )
    {
      wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
    }
  
  //  free( gam );
}// end of tridag


inline void Tridag_Fast( const complex<double>  &tridag_lower_diagonal , const vector<complex<double> > &tridag_middle_diagonal , const complex<double>  &tridag_upper_diagonal, const vector<complex<double> > &wf_1D_rightside ,
		  vector<complex<double> > &wf_1D_tmp2 , vector<complex<double> > &gam )
{   
  const unsigned int length=gam.size() - 2;

/*  if ( tridag_middle_diagonal[ 1 ] == 0. )
    {
      nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    }
*/

  complex<double>  bet=tridag_middle_diagonal[1];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length ; ++index )
    {
      gam[ index ] = tridag_upper_diagonal/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal*gam[ index ];
/*      if (bet==0.0)
	{
	nerror("Error: bet==0 in tridag!");
	}
*/
      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
  for( unsigned int index = length-1 ; index >= 1 ; --index )
    {
      wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
    }
  
 
}// end of tridag

inline void Tridag_FastAM( const complex<double>  &tridag_lower_diagonal , const vector<complex<double> > &tridag_middle_diagonal , const complex<double>  &tridag_upper_diagonal, const vector<complex<double> > &wf_1D_rightside ,
		  vector<complex<double> > &wf_1D_tmp2 , vector<complex<double> > &gam , int length)
{   


  //  if ( tridag_middle_diagonal[ 1 ] == 0. )
  //   nerror( "Error: tridag_middle_diagonal[1]==0 in tridag!" );
    


  complex<double>  bet=tridag_middle_diagonal[1];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length ; ++index )
    {
      gam[ index ] = tridag_upper_diagonal/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal*gam[ index ];
      //      if (bet==0.0)
      //	 nerror("Error: bet==0 in tridag!");
	

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
  for( unsigned int index = length-1 ; index >= 1 ; --index )
       wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
 
  
  
}// end of tridag



/**
 *  The tridag routine to solve the threediagonal system of linear equation for ADAPATATIVE MESH
 */
inline void TridagAM( const vector<complex<double> > &tridag_lower_diagonal , const vector<complex<double> > &tridag_middle_diagonal , const vector<complex<double> >&tridag_upper_diagonal, const vector<complex<double> > &wf_1D_rightside ,
		      vector<complex<double> > &wf_1D_tmp2 , vector<complex<double> > &gam, int length )
{   


/**
 *  forth substitution
 */

  complex<double>  bet=tridag_middle_diagonal[ 1 ];
  wf_1D_tmp2[ 1 ] = wf_1D_rightside[ 1 ]/bet;
  for( unsigned int index = 2 ; index <= length; ++index )
    {
      gam[ index ] = tridag_upper_diagonal[ index-1 ]/bet;
      bet  = tridag_middle_diagonal[ index ];
      bet -= tridag_lower_diagonal[ index ]*gam[ index ];
      /*      if ( bet == 0.0 )
	{
	 cout<<"Error: bet==0 in tridag! at "<<index<<"\n";
         cout<<"tridag_middle_diagonal="<<tridag_middle_diagonal[ index ]<<"\n";
         cout<<"tridag_lower_diagonal="<<tridag_lower_diagonal[ index ]<<"\n";
         cout<<"gam="<<gam[ index ]<<"\n";
	 } */

      wf_1D_tmp2[ index ] = ( wf_1D_rightside[ index ]-tridag_lower_diagonal[ index ]*wf_1D_tmp2[ index-1 ] )/bet;
    }
  
/**
 *  back substitution
 */
  for( unsigned int index = length-1 ; index >= 1 ; --index )
       wf_1D_tmp2[ index ] -= gam[ index+1 ]*wf_1D_tmp2[ index+1 ];
 
}// end of tridag











/*
void tridag(complex<double>  a,vector<complex<double> > b,complex<double>  c,vector<complex<double> > r,vector<complex<double> > &u, int n)
{
    int j;
    complex<double>  bet, *gam;
    gam=(complex<double>  *)malloc((n+2)*sizeof (complex<double> ));
    if(!gam) cout<<"\nerror de colocacion en la asignacion complex<double> ";
    if(b[1]==0.0) nrerror ("error 1 en tridag");
    u[1]=r[1]/(bet=b[1]);
    for(j=2;j<=n;j++)
    {
		gam[j]=c/bet;
		bet=b[j]-a*gam[j];
		if (bet==0.0) nrerror("error 2 en tridag");
		u[j]=(r[j]-a*u[j-1])/bet;
    }
    for(j=(n-1);j>=1;j--)
		u[j]-=gam[j+1]*u[j+1];
	
    //for(j=1;j<=n;j++)
    //  printf("%d %e\n",j,real(u[j]));
    free(gam);
    //    printf("I am here\n");
	
}
*/

inline void skip_comments(std::istream & is)
{
  bool comments=true;
  char c;
  std::string linebuffer;
  while (comments) {
    is.get(c);
    if ('#'==c)
      getline(is,linebuffer);
    else {
        is.unget();
        comments = false;
    }
  }
}

inline void set_precision(std::ostream &os, size_t p) {
  os.precision(p);
  os.setf(std::ios::scientific, std::ios::floatfield);
  os.setf(std::ios::showpoint);
}


/**
 *  convert from integer to string
 */

inline string Conversion_To_String( int i )
{
  ostringstream ost;
  ost << i;
  return ost.str();
} // to_string integer


/**
 *  convert from double to string
 */
 
inline string Conversion_To_String( double i )
{
  ostringstream ost;
  ost << i;
  return ost.str();
} // to_string double
 
/*
// shaohao2009.04

// FFT functions based on fftw3 libs
void fft_real_1D (const int n, vector<double> &inf, vector<double> &outf)
{
  fftw_plan p_real_1d = fftw_plan_r2r_1d(n, &inf[0], &outf[0], FFTW_RODFT11, FFTW_ESTIMATE);

  fftw_execute(p_real_1d);

  fftw_destroy_plan(p_real_1d);

}

// fft for complex<double>  1D
void fft_complex_1D (const int &n, vector<complex<double> > &input_f, vector<complex<double> > &output_f)
{
  fftw_complex *in;
  fftw_complex *out;

  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

// plan first, fill later
  fftw_plan p_complex_1d = fftw_plan_dft_1d(n, &in[0], &out[0], FFTW_FORWARD, FFTW_ESTIMATE);

// fill input
  for(int i=0; i<n; i++){
    in[i][0]=real(input_f[i]);     // for fftw3, shaohao 2009.01
    in[i][1]=imag(input_f[i]);
  }

// excute
  fftw_execute(p_complex_1d);

// fill output
  for(int i=0; i<n; i++){
    output_f[i] = out[i][0] + I*out[i][1];
  }

// destroy and free
  fftw_destroy_plan(p_complex_1d);
  free(in);  free(out);

}


// ft for complex<double>  1D
enum dump_type { no_dump, exp_dump, poly_dump };
enum ft_type { odd, even, forward, backward };

void ft_complex_1D (vector<complex<double> > &input, vector<double> &input_x, const int n_inp, const double inpstep, vector<complex<double> > &output, vector<double> &out_x, const double outmin, const double outmax, const double outstep, const dump_type dump, const ft_type type, const double damp_f)
{


  vector<double> input_r(n_inp), input_i(n_inp);
  for( int i=0; i < n_inp; i++)
  {
    input_r[i] = real(input[i]);
    input_i[i] = imag(input[i]);
  }

  int n_out =int((outmax-outmin)/outstep);
  output.resize(n_out);  out_x.resize(n_out);
  vector<double> output_r(n_out), output_i(n_out);
  vector<double> dump_factor(n_out);
  double x_odd = 0., x_even = 0.;
// const std::complex<double>  I=std::complex<double> (0.,1.);
  complex<double>  y = complex<double> (0.,0.);

  for( int j=0; j < n_out; j++)
  {

    out_x[j] = outmin + j * outstep;
    output_r[j] =0.; output_i[j] =0.;
    output[j] = complex<double> (0.,0.);

    for( int i=0; i < n_inp; i++)
    {
// dump type  
       if( dump == no_dump )
         dump_factor[i] = 1.0;
       else if( dump == exp_dump )
         dump_factor[i] = exp(-i * inpstep * damp_f);
       else if( dump == poly_dump )
         dump_factor[i] = 1.0 - 3.0*(i/n_inp)*(i/n_inp) + 2.0*(i/n_inp)*(i/n_inp)*(i/n_inp);
       else
         cout<<"dump factor error!"<<endl;

// ft type
       if( type == odd){
          x_odd = sin(out_x[j] * input_x[i]);
          output_r[j] += input_r[i] * x_odd * dump_factor[i];
          output_i[j] += input_i[i] * x_odd * dump_factor[i];
      }
      else if( type == even){
          x_even = cos(out_x[j] * input_x[i]);
          output_r[j] += input_r[i] * x_even * dump_factor[i];
          output_i[j] += input_i[i] * x_even * dump_factor[i];
      }
     else if( type == forward){
          x_odd = sin(out_x[j] * input_x[i]);
          x_even = cos(out_x[j] * input_x[i]);
          output_r[j] += x_even * input_r[i] + x_odd * input_i[i];
          output_i[j] += x_even * input_i[i] - x_odd * input_r[i];
       }
      else if( type == backward){
          x_odd = sin(out_x[j] * input_x[i]);
          x_even = cos(out_x[j] * input_x[i]);
          output_r[j] += x_even * input_r[i] - x_odd * input_i[i];
          output_i[j] += x_even * input_i[i] + x_odd * input_r[i];
       }

//       else if( type == forward)
//          y = exp ( -I * out_x[j] * input_x[i]);
//       else if( type == backward)
//          y = exp ( I * out_x[j] * input_x[i]);
       else
          cout<<"ft type error!"<<endl;

    }  // end of iteration of input data
// normalization

    if ( type == odd or even){
       output_r[j] = 2.0*output_r[j] / n_inp;
       output_i[j] = 2.0*output_i[j] / n_inp;
       output[j] = output_r[j] + I * output_i[j];
    }
    else if ( type == forward or backward ) {
       output_r[j] = output_r[j] / sqrt(n_inp);
       output_i[j] = output_i[j] / sqrt(n_inp);
       output[j] = output_r[j] + I * output_i[j];
    }
    else
        cout<<"ft type error!"<<endl;

  }  // end of iteration of output data

}
*/

#endif	/* UTILS_H */
