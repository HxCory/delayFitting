// $Id: Projection.h 557 2007-06-12 16:45:37Z silvio $
#include <cstdlib>
#ifndef PROJECTION_H
#define PROJECTION_H
#include "wavefunction.h"

using std::cout;
using std::cerr;
using std::string;
using std::endl;

class Projection {
public:
  size_t size, rows, cols;
  vec<double> data;
  Axis axis_row, axis_col;
  /// the "title" is set by Output/MatOutput to recognize the Projection
  // and prevent it to be written into the wrong file.
  string title;
  /// Default constructor for late initialization
  Projection() : rows(0), cols(0) {}
  void axis_check(const Axis &axis_1, const Axis &axis_2) {
    size_t _rows = axis_1.size();
    size_t _cols = axis_2.size();
    // check if rows or cols have changed
    if(rows != _rows || cols != _cols) {
      // if they were zero until now, allow them to be changed
      if (0 == rows && 0 == cols) {
        rows = _rows; cols = _cols;
        size=rows*cols;
        data.resize(rows*cols,0.);
        axis_row = axis_1;
        axis_col = axis_2;
      } else { // report an error
        cerr<<"Error: Projection "<< title <<" is ("<<rows<<"x"<<cols<<")"<<endl;
        cerr<<"       Cannot use it to store ("<<_rows<<"x"<<_cols<<") data!"<<endl;
        exit(1);
      }
    }
  }
  
  void Wigner_axis_check(const Axis &axis_1) {
    size_t _rows = axis_1.size();
    // check if rows have changed
    if(rows != _rows) {
      // if it was zero until now, allow it to be changed
      if (0 == rows) {
        rows = _rows; cols=_rows;
        size=rows*cols;
        data.resize(size,0.);
        axis_row = axis_1;
        // setup Momentum axis for Wigner map
        int N = cols-2; // calculation without the border
        double delta=pi/(axis_1.delta*N);
        axis_col.Init("p_" + axis_row.label, N, delta, FFTWAxis);

      } else { // report an error
        cerr<<"Error: Projection "<< title <<" is ("<<rows<<"x"<<cols<<")"<<endl;
        cerr<<"       Cannot use it to store ("<<_rows<<"x"<<_rows<<") data!"<<endl;
        exit(1);
      }
    }
  }
  
  Projection& operator-=(Projection &rhs) {
    if (rhs.data.size() != size) {
      cerr <<"Trying to subtract projections of unequal size"<<endl;
      exit(1);
    }
    for (int i=0; i<size; ++i) {
      data[i]-=rhs.data[i];
    }
    return *this;
  };
  Projection& operator+=(Projection &rhs) {
    if (rhs.data.size() != size) {
      cerr <<"Trying to add projections of unequal size"<<endl;
      exit(1);
    }
    for (int i=0; i<size; ++i) {
      data[i]+=rhs.data[i];
    }
    return *this;
  };
};

class Projection1D {
public:
  size_t size;
  vec<double> data;
  Axis axis;
  /// the "title" is set by Output/MatOutput to recognize the Projection
  // and prevent it to be written into the wrong file.
  string title;
  /// Default constructor for late initialization
  Projection1D() : size(0) {}
  void axis_check(const Axis & axis_) {
    size_t _size = axis_.size();
    // check if size have changed
    if(size != _size) {
      // if it was zero until now, allow it to be changed
      if (0 == size) {
        size = _size;
        data.resize(size,0.);
        axis=axis_;
      } else { // report an error
        cerr<<"Error: Projection "<< title <<" is ("<<size<<")"<<endl;
        cerr<<"       Cannot use it to store ("<<_size<<") data!"<<endl;
        exit(1);
      }
    }
  }
  
  Projection1D& operator-=(Projection1D &rhs) {
    if (rhs.data.size() != size) {
      cerr <<"Trying to subtract projections of unequal size"<<endl;
      exit(1);
    }
    for (int i=0; i<size; ++i) {
      data[i]-=rhs.data[i];
    }
    return *this;
  };
  Projection1D& operator+=(Projection1D &rhs) {
    if (rhs.data.size() != size) {
      cerr <<"Trying to add projections of unequal size"<<endl;
      exit(1);
    }
    for (int i=0; i<size; ++i) {
      data[i]+=rhs.data[i];
    }
    return *this;
  };
};

#endif /* PROJECTION_H */
