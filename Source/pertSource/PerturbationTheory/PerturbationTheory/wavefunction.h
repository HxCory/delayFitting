// $Id: geometry.h 229 2006-06-16 13:00:51Z arvid $
#include <cstdlib>
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H
#include "constant.h"
#include "ParameterMap.h"

void touchdir(string name);  // defined in OutputStuff.cpp , used in wf.save()

using namespace std; // to save the std:: in front of vector and cout

/**
 *  The class wavefunction contains all the information related to the wavefunction.
 */

template<class T>
class vec : public vector<T> {
#ifdef DEBUG
  public:
  const T& operator[](const size_t index) const {
    return vector<T>::at(index);
  }
  T& operator[](const size_t index) {
    return vector<T>::at(index);
  }
#endif
};

// shaohao2009.06
enum AxisType { FullAxis, HalfAxis, FFTWAxis, FFTWHalfAxis, PartAxis, FFTWPartAxis };  /**< Axis Type: Full, Half or FFTW-ordering */

class Axis : public vec<double> {
  public:
  string label;
  double delta;   /**< Spatial Step */  
  AxisType type;  /**< Axis Type: Full, Half or FFTW-ordering */

  // convenience functions:
  // Copy constructor:
  Axis(const Axis &rhs) :
    vec<double>(rhs), label(rhs.label), delta(rhs.delta), type(rhs.type) {}

  // Assignment operator
  Axis & operator=(const Axis &rhs) {
    if (this == &rhs) return *this;   // Gracefully handle self assignment
    // call the operator= of the base class (vector<double>) to copy its data
    this->vec<double>::operator=(rhs);
    // finally copy all additional data
    label=rhs.label;
    delta=rhs.delta;
    type=rhs.type;
    return *this;
  }

  // Initialization
  void Init(string _label, int N_points, double _delta, AxisType _type = FullAxis) {

    label=_label;
    delta=_delta;
    type=_type;

// for one-side absorber
        double a_percent = 0.03/10.; // 0.1 for test, 0.03/10. , for real calculation;
        int N_right = 10000; // 20 for test, 10000 for real calculation   // only for the case of H atom in 2D Cylindrical coordinates with 10000 grids on the right.
        int N_a = int(N_points*a_percent);
        int N_left = N_points - N_a - N_right;
        double l_percent = (double(N_left) + double(N_a))/double(N_points);

    this->resize(N_points+2, 0.); // convention: 'zero' boundary points for derivatives

    switch (type) {
      case FullAxis:
        /**
         * no symmetry; the whole grid is needed
         */
        for( int i = 1 ; i <= N_points ; i++ ) {
          (*this)[ i ] = ( -( N_points+1 )/2.+i )*delta;
        }
        break;
      case HalfAxis:
        /**
         * symmetric(=1) or antisymmetric(=-1) wavefunction; only half the grid is needed
         */
        for( int i = 1 ; i <= N_points ; i++ ) {
          (*this)[ i ] = ( i-0.5 )*delta;
        }
        break;
      case FFTWAxis:
        /**
         * FFT ordering
         * (mangle the axis ordering to match the convention of fftw: from plus to minus)
         */
        for( int i = 1 ; i <= int(N_points/2) +1 ; i++ ) {
          (*this)[ i ] =   (i-1)*delta;
        }
        for( int i = int(N_points/2)+2 ; i <= N_points ; i++ ) {
          (*this)[ i ] =   (i-1-N_points)*delta;
        }
        break;
      case FFTWHalfAxis:
        /**
         * FFT ordering
         * (mangle the axis ordering to match the convention of fftw)
         */
        for( int i = 1 ; i <= N_points ; i++ ) {
          (*this)[ i ] =   (i-1)*delta;
        }
        break;
// shaohao2009.06
      case PartAxis:
        /**
         * no symmetry; the grids on the left is less than those of the right
         */
         cout<<"---sh---"<< N_right <<'\t'<< N_left <<'\t'<< N_a <<'\t'<< N_points << endl;
         cout<<"---sh---"<< a_percent <<'\t'<< l_percent << endl;
        for( int i = 1 ; i <= N_points ; i++ ) {
          (*this)[ i ] = ( -( N_points+1 )*l_percent + i )*delta;
        }
        break;
      case FFTWPartAxis:
        /**
         * no symmetry; the grids on the left is less than those of the right; FFT order
         */
        for( int i = 1 ; i <= N_right +1 ; i++ ) {
          (*this)[ i ] =   (i-1)*delta;
        }
        for( int i = N_right+2 ; i <= N_points ; i++ ) {
          (*this)[ i ] =   (i-1-N_points)*delta;
        }
        break;
    }

  };

  // Default constructor : do nothing, Init later
  Axis() {};
};

class wavefunction { 
  public:
  string name;
  /// the data of the wavefunction
  // Note: we allocate 7 units of each and we don't use the first one
  //       to start intuitive counting at 1
  int n[DIM+1];/**< Number of points in the N[j] dimension */  
  Axis x[DIM+1];
  Axis p[DIM+1];

  vec<complex<double> > wave;/**< the Wavefunction, a complex<double>  1D array, size (n1+2)*(n2+2)*(n3+2)... */

  /// references for compatibility and short notation
  double &dx1;/**< Spatial Step x1( i ) */  
  double &dx2;/**< Spatial Step x2( i ) */  
  double &dx3;/**< Spatial Step x3( i ) */  
  double &dx4;/**< Spatial Step x4( i ) */
  double &dx5;/**< Spatial Step x5( i ) */  
  double &dx6;/**< Spatial Step x6( i ) */    
  double &dp1;/**< Spatial Step p1( i ) */  
  double &dp2;/**< Spatial Step p2( i ) */  
  double &dp3;/**< Spatial Step p3( i ) */  
  double &dp4;/**< Spatial Step p4( i ) */
  double &dp5;/**< Spatial Step p5( i ) */  
  double &dp6;/**< Spatial Step p6( i ) */    
  double one_by_dx1, one_by_dx1sqr;  // helper variables
  double one_by_dx2, one_by_dx2sqr;
  double one_by_dx3, one_by_dx3sqr;

  int &n1;  /**< Number of points in the N1 dimension */  
  int &n2;  /**< Number of points in the N2 dimension */  
  int &n3;  /**< Number of points in the N3 dimension */  
  int &n4;  /**< Number of points in the N4 dimension */  
  int &n5;  /**< Number of points in the N5 dimension */  
  int &n6;  /**< Number of points in the N6 dimension */  
  Axis &x1; /**< vector x1( i ) */  
  Axis &x2; /**< vector x2( i ) */
  Axis &x3; /**< vector x3( i ) */
  Axis &x4; /**< vector x4( i ) */
  Axis &x5; /**< vector x5( i ) */
  Axis &x6; /**< vector x6( i ) */
  Axis &p1; /**< vector p1( i ) */  
  Axis &p2; /**< vector p2( i ) */
  Axis &p3; /**< vector p3( i ) */
  Axis &p4; /**< vector p4( i ) */
  Axis &p5; /**< vector p5( i ) */
  Axis &p6; /**< vector p6( i ) */
  
  bool verbose; /** Tune the level of output */  
  int symmetry_x1; /** kind of symmetry of the wavefunction with respect to x1 */
  int symmetry_x2; /** kind of symmetry of the wavefunction with respect to x2 */
  int symmetry_x3; /** kind of symmetry of the wavefunction with respect to x3 */
  /**
   * A index function in 2D to acces the [j][i] number in the array  
   * Transform two indexs in one, allowing the use of a 1D arrays.
   * in2 is only useful if n1>1 and n2>1
   */

  int in2( int a , int b ) const
  {
    // Storage Order is Column-major !
    //return a*( n1+2 )+b;
    int index=a;
    index*=n1+2;
    index+=b;
    return index;
  }

  /**
   * A index function in 3D to acces the [k][j][i] number in the array  
   * Transform three indexs in one, allowing the use of a 1D arrays.
   * in2 is only useful if n1>1 , n2>1 and n3>1
   */

  int in3( int a , int b , int c ) const
  {
    // Storage Order is Column-major !
    //return a*( n2+2 )*( n1+2 )+b*( n1+2 )+c; 
    int index=a;
    index*=n2+2;
    index+=b;
    index*=n1+2;
    index+=c;
    return index;
  }

  ofstream log; //The log file
  int left_x1, right_x1; // The limits of the grid,. useful for Adaptatuve grid
  int left_x2, right_x2;
  int left_x3, right_x3;

  void dump(const string & filename) {
    // binary opening only means, that line ends are not converted.
    ofstream ofs;
    ofs.open(filename.c_str(), ios::out | ios::binary);
    if(ofs.fail()) {
      cerr<<"Could not open output file "<< filename<<endl;
      exit(1);
    }
    size_t width=sizeof(complex<double> );
    for(size_t i=0; i<wave.size(); i++)
      ofs.write((char*)&wave[i], width);
    ofs.close();
    if(! ofs.fail() ) cout<<"Wrote "<<filename<<endl;
  }
  bool restore(const string & filename) {
    ifstream ifs;
    ifs.open(filename.c_str(), ios::in | ios::binary);
    if(ifs.fail()) {
      cerr<<"Could not open imput file "<< filename<<endl;
    }
    size_t width=sizeof(complex<double> );
    for(size_t i=0; i<wave.size(); i++)
      ifs.read((char*)&wave[i], width);
    ifs.close();
    if(! ifs.fail() ) {
      cout<<"Read "<<filename<<endl;
      return true;
    } else return false;
  }

  void save(const string &_label) {
    ParameterMap main_parameters;
    string dirname = *main_parameters["wavefunctionFolder"];

    touchdir(dirname.c_str());
    string basename=dirname+"/"+_label;
    string txtname=basename+".txt";
    ofstream ofs(txtname.c_str());
    name=_label;
    ofs<<"wfname = "<<_label<<endl;
    ofs<<"n_points = ("<<n1<<", "<<n2<<", "<<n3<<")"<<endl;
    ofs<<"spatial_step = ("<<dx1<<", "<<dx2<<", "<<dx3<<")"<<endl;
    ofs<<"symmetry_x2 = ("<<symmetry_x2<<")"<<endl;
    ofs<< main_parameters;
    //ofs<<"R_kz = "<<*main_parameters["R_kz"]<<endl;

    dump(basename+".dat");
  }

  bool load(const string &_label) {
    ParameterMap main_parameters;
    string dirname = *main_parameters["wavefunctionFolder"];

    string basename=dirname+"/"+_label;
    string txtname=basename+".txt";
   

    try {   // need a try-catch block here because ConfigFileParser throws
            // an exception if the file cannot be opened.
      ConfigFileParser metadata( txtname.c_str() );
      
      {
        FromStringToVector<int > v(metadata["n_points"]);
        if ( (n1!=0 && n1!=v[0]) /*||
             (n2!=0 && n2!=v[1]) ||
             (n3!=0 && n3!=v[2])*/  ) {
          cerr <<"wf.load: n_points does not match!"<<endl;
          exit(1);
        }  
	FromStringToVector<double > vv(metadata["spatial_step"]);
	if ( (dx1!=0 && dx1!=vv[0]) /*||
             (dx2!=0 && dx2!=vv[1]) ||
             (dx3!=0 && dx3!=vv[2])*/  ) {
          cerr <<"wf.load: dx does not match!"<<endl;
          exit(1);
        } 
	//FromStringToVector<int> vvv(metadata["symmetry_x2"]); 
	//if ( (symmetry_x2!=vvv[0])  ) {
          //cerr <<"wf.load: symmetry does not match!"<<endl;
          //exit(1);
        //} 

else {
          n1=v[0]; n2=v[1], n3=v[2];
        }
      }

      {
        FromStringToVector<double > v(metadata["spatial_step"]);
        if ( (dx1!=0 && dx1!=v[0]) /*||
             (dx2!=0 && dx2!=v[1]) ||
             (dx3!=0 && dx3!=v[2])*/ ) {
          cerr <<"wf.load: spatial_step does not match!"<<endl;
          exit(1);
        } else {
          dx1=v[0]; dx2=v[1], dx3=v[2];
        }
      }

      /*
      double r_old=*main_parameters["R_kz"];
      double r_new=FromStringTo<double>(metadata["R_kz"]);
      if (r_new != r_old) {
        cerr <<"wf.load: r_new != r_old"<<endl;
        exit(1);
      }
      */

      name=metadata["wfname"];

      return restore(basename+".dat");
    } catch (ConfigFileParser::FileNotFound) {
      return false;
    }
  }

  /// Constructor needed to initialize the shorthand references
  wavefunction() :
    dx1(x[1].delta), dx2(x[2].delta), dx3(x[3].delta), dx4(x[4].delta), dx5(x[5].delta), dx6(x[6].delta),
    dp1(p[1].delta), dp2(p[2].delta), dp3(p[3].delta), dp4(p[4].delta), dp5(p[5].delta), dp6(p[6].delta),
    n1(n[1]), n2(n[2]), n3(n[3]), n4(n[4]), n5(n[5]), n6(n[6]),
    x1(x[1]), x2(x[2]), x3(x[3]), x4(x[4]), x5(x[5]), x6(x[6]),
    p1(p[1]), p2(p[2]), p3(p[3]), p4(p[4]), p5(p[5]), p6(p[6])
    {}

  wavefunction(const wavefunction &rhs) : /* Copy constructor */
    wave(rhs.wave),
    verbose(rhs.verbose),
    symmetry_x1(rhs.symmetry_x1),
    symmetry_x2(rhs.symmetry_x2),
    symmetry_x3(rhs.symmetry_x3),
    name(rhs.name),
    n1(rhs.n1), n2(rhs.n2), n3(rhs.n3), n4(rhs.n4), n5(rhs.n5), n6(rhs.n6),
    dx1(rhs.dx1), dx2(rhs.dx2), dx3(rhs.dx3), dx4(rhs.dx4), dx5(rhs.dx5), dx6(rhs.dx6),
    one_by_dx1(rhs.one_by_dx1),
    one_by_dx2(rhs.one_by_dx2),
    one_by_dx3(rhs.one_by_dx3),
    one_by_dx1sqr(rhs.one_by_dx1sqr),
    one_by_dx2sqr(rhs.one_by_dx2sqr),
    one_by_dx3sqr(rhs.one_by_dx3sqr),
    dp1(rhs.dp1), dp2(rhs.dp2), dp3(rhs.dp3), dp4(rhs.dp4), dp5(rhs.dp5), dp6(rhs.dp6),
    x1(rhs.x1), x2(rhs.x2), x3(rhs.x3), x4(rhs.x4), x5(rhs.x5), x6(rhs.x6),
    p1(rhs.p1), p2(rhs.p2), p3(rhs.p3), p4(rhs.p4), p5(rhs.p5), p6(rhs.p6)
  {
    for(int i=0; i<DIM+1; i++) {
      n[i]=rhs.n[i];
      x[i]=rhs.x[i];
      p[i]=rhs.p[i];
    }
  };

  wavefunction &operator*=(const double rhs) {
    const int off=wave.size();
    for(int i=0; i<off; i++) {
      wave[i]*=rhs;
    }
    return *this;
  };

  wavefunction &operator+=(const wavefunction &rhs) {
    for(int i=0; i<DIM+1; i++) {
      if(rhs.x[i].size() != x[i].size()) {
        cerr<<"Cannot add wavefunctions of different size (x["<<i<<"]"<<endl;
        exit(1);
      }
      if(rhs.p[i].size() != p[i].size()) {
        cerr<<"Cannot add wavefunctions of different size (p["<<i<<"]"<<endl;
        exit(1);
      }
    }
    const int off=wave.size();
    for(int i=0; i<off; i++) {
      wave[i]+=rhs.wave[i];
    }
    return *this;
  };

  wavefunction &operator-=(const wavefunction &rhs) {
    for(int i=0; i<DIM+1; i++) {
      if(rhs.x[i].size() != x[i].size()) {
        cerr<<"Cannot add wavefunctions of different size (x["<<i<<"]"<<endl;
        exit(1);
      }
      if(rhs.p[i].size() != p[i].size()) {
        cerr<<"Cannot add wavefunctions of different size (p["<<i<<"]"<<endl;
        exit(1);
      }
    }
    const int off=wave.size();
    for(int i=0; i<off; i++) {
      wave[i]-=rhs.wave[i];
    }
    return *this;
  };

}; 

inline wavefunction operator+(const wavefunction &lhs, const wavefunction &rhs)
{
  wavefunction r(lhs);
  return r+=rhs;
};


inline wavefunction operator-(const wavefunction &lhs, const wavefunction &rhs)
{
  wavefunction r(lhs);
  return r-=rhs;
};


#endif /* WAVEFUNCTION_H */
