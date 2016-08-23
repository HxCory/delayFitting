// $Id: Observable.h 546 2007-05-15 14:15:40Z arvid $
#include <cstdlib>
#ifndef OBSERVABLE_H
#define OBSERVABLE_H
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <complex>
#include <iomanip>
#include "ParameterMap.h"
#include "Projection.h"

using namespace std;

extern Parameter<string> outputFolder; // it's defined outside

/**
 * The next line defines Output Object as a template.
 */
template<class T> class Output;

/**
 * The next few lines provide a specializiation of Output for double
 */
template<>
class Output<double>{
  string label;
  ofstream txt_ofs;            // the text output file for the time, "# keyword = value" comments and descriptive comments
  public:
  Output(string _label) : label(_label) {};
  void openFiles() {
    string basename=outputFolder()+label;
    string txtname=basename+".txt";
    txt_ofs.open(txtname.c_str());
    if(txt_ofs.fail()) {
      cerr<<"Could not open output file "<< txtname<<endl;
      exit(1);
    }
    txt_ofs.setf(ios::showpoint|ios::scientific);
    ParameterMap pMap;
    txt_ofs << pMap;                           // write all parameters in the header.
    txt_ofs << "#  time " << label << endl;
  }
  ~Output() {
    txt_ofs.close();
  }
  inline void save(double time, double val) {
    if(!txt_ofs.is_open()) openFiles();   // lazy opening,just before first write
    txt_ofs << time << "\t" << val << endl;     // both into the text file
  }
};

/**
 * The next few lines provide a specializiation of Output for Projection*
 */
template<>
class Output<Projection &> {
  size_t sample;
  string title;
  string basename;
  ofstream ofs;
  ofstream t_ofs;
  public:
  Output(string _title) : title(_title) {
    basename=outputFolder()+title;
  };


  void openFile() {
    {
      string fullname=basename+".dat";
      // binary opening only means, that line ends are not converted.
      ofs.open(fullname.c_str(), ios::out | ios::binary);
      if(ofs.fail()) {
        cerr<<"Could not open output file "<< fullname<<endl;
        exit(1);
      }
      ofs.setf(ios::showpoint|ios::scientific);
    }

    // write Metadata to separate txt files
    // open t_file
    {
      string t_name=basename+".txt";
      t_ofs.open(t_name.c_str());
      if(t_ofs.fail()) {
        cerr<<"Could not open output file "<< t_name<<endl;
        exit(1);
      }
      t_ofs.setf(ios::showpoint|ios::scientific);
    }

    ParameterMap pMap;
    t_ofs << pMap;                           // write all parameters in the header.

    ostringstream tmp; tmp << title <<"_t";
    t_ofs << "# "<< endl;
    t_ofs << "# " << tmp.str() << endl;
  }


  void write_axis(Projection &prj) {
    prj.title=title;          // remember to which file the projection was saved
    // write Metadata to separate txt files
    // 1. store row-Axis
    {
      string t_name=basename+"_"+prj.axis_row.label+".txt";
      t_ofs.open(t_name.c_str());
      if(t_ofs.fail()) {
        cerr<<"Could not open output file "<< t_name<<endl;
        exit(1);
      }
      t_ofs.setf(ios::showpoint|ios::scientific);
      //ostringstream tmp; tmp << title <<"_"<<prj.row_label;
      //t_ofs << "# " << tmp.str() << endl;
      t_ofs << "# " << title+"_"+prj.axis_row.label << endl;
      for(size_t i=0; i<prj.axis_row.size(); i++)
        t_ofs << prj.axis_row[i] << endl;
      t_ofs.close();
    }

    // 2. store col-Axis
    {
      string t_name=basename+"_"+prj.axis_col.label+".txt";
      t_ofs.open(t_name.c_str());
      if(t_ofs.fail()) {
        cerr<<"Could not open output file "<< t_name<<endl;
        exit(1);
      }
      t_ofs.setf(ios::showpoint|ios::scientific);
      //ostringstream tmp; tmp << title <<"_"<<prj.col_label;
      //t_ofs << "# " << tmp.str() << endl;
      t_ofs << "# " << title+"_"+prj.axis_col.label << endl;
      for(size_t i=0; i<prj.axis_col.size(); i++)
        t_ofs << prj.axis_col[i] << endl;
      t_ofs.close();
    }

  }


  void save(double t, Projection &prj) {
    if(!t_ofs.is_open()) {
      openFile();   // lazy opening,just before first write
      write_axis(prj);
    }
    if(prj.title != title) {           // sanity check
      cerr<<"Error: writing two different projections to one Output file"<<endl;
      cerr<<"       first "<<title<<" and now "<<prj.title<<endl;
      exit(1);
    }

    t_ofs << t << endl;             // log the time
    // write the value
    size_t width=sizeof(prj.data[0]);
    for(size_t i=0; i<prj.data.size(); i++)
      ofs.write((char*)&prj.data[i], width);
  }


  ~Output() {
    ofs.close();
    t_ofs.close();
  }
};

/**
 * The next few lines provide a specializiation of Output for vector<complex<double> >
 * outputs in Camilos binary format
 * Note: This will perform badly, since the vector<complex<double> > will be copied!
 */
template<class T>
class Output<vector<T> >{
  string label;
  //ofstream bin_ofs;                 // the output file for the observable, possibly binary
  ofstream txt_ofs;            // the text output file for the time, "# keyword = value" comments and descriptive comments
  public:
  Output(string _label) : label(_label) {};
  void openFiles() {
    string basename=outputFolder()+label;
    /*
    string fullname=basename+".dat";
    // binary opening only means, that line ends are not converted.
    bin_ofs.open(fullname.c_str(), ios::out | ios::binary);
    if(bin_ofs.fail()) {
      cerr<<"Could not open output file "<< fullname<<endl;
      exit(1);
    }
    bin_ofs.setf(ios::showpoint|ios::scientific);
    */
    string txtname=basename+".txt";
    txt_ofs.open(txtname.c_str());
    if(txt_ofs.fail()) {
      cerr<<"Could not open output file "<< txtname<<endl;
      exit(1);
    }
    txt_ofs.setf(ios::showpoint|ios::scientific);
    ParameterMap pMap;
    txt_ofs << pMap;                           // write all parameters in the header.
    txt_ofs << "#  time" << endl;
  }
  ~Output() {
    //bin_ofs.close();
    txt_ofs.close();
  }
  inline void save(double time, vector<T> vect) {
    if(!txt_ofs.is_open()) openFiles();   // lazy opening,just before first write
    /*
    txt_ofs << time << endl;             // log the time
    // write the value
    size_t width=sizeof(T);
    for(size_t i=0; i<vect.size(); i++)
      bin_ofs.write((char*)&vect[i], width);
    */
    txt_ofs << time << "  ";             // log the time
    for(size_t i=0; i < vect.size() - 1; i++)
      txt_ofs << vect[i] << "  ";
    txt_ofs << vect.back() << endl;
  }
};

/**
 * The next few lines provide a specializiation of Output for vector<complex<double> >
 * outputs in Camilos binary format
 * Note: This will perform badly, since the vector<complex<double> > will be copied!
 */
template <>
class Output<vector<complex<double> > >{
  string label;
  ofstream bin_ofs;                 // the output file for the observable, possibly binary
  ofstream txt_ofs;            // the text output file for the time, "# keyword = value" comments and descriptive comments
  public:
  Output(string _label) : label(_label) {};
  void openFiles() {
    string basename=outputFolder()+label;
    string fullname=basename+".dat";
    // binary opening only means, that line ends are not converted.
    bin_ofs.open(fullname.c_str(), ios::out | ios::binary);
    if(bin_ofs.fail()) {
      cerr<<"Could not open output file "<< fullname<<endl;
      exit(1);
    }
    bin_ofs.setf(ios::showpoint|ios::scientific);
    string txtname=basename+".txt";
    txt_ofs.open(txtname.c_str());
    if(txt_ofs.fail()) {
      cerr<<"Could not open output file "<< txtname<<endl;
      exit(1);
    }
    txt_ofs.setf(ios::showpoint|ios::scientific);
    ParameterMap pMap;
    txt_ofs << pMap;                           // write all parameters in the header.
    txt_ofs << "#  time" << endl;
  }
  ~Output() {
    bin_ofs.close();
    txt_ofs.close();
  }
  inline void save(double time, vector<complex<double> > vect) {
    if(!txt_ofs.is_open()) openFiles();   // lazy opening,just before first write
    txt_ofs << time << endl;             // log the time
    // write the value
    size_t width=sizeof(complex<double> );
    for(size_t i=0; i<vect.size(); i++)
      bin_ofs.write((char*)&vect[i], width);
  }
};

class Sampler{
  size_t every;
  size_t count;
  vector<size_t> at;
public:
  size_t sample;
  /// Conversion operator to bool
  bool operator() () {                // the purpose of this: counter logic
    bool now=false;
    if(0 == at.size()) {            // ah, the simple Constructor was used
      count++;
      if(count == every) {
        count = 0;                  // reset counter;
        now= true;
        sample++;
      }
    } else {                        // oh, the array was given
      if(sample < at.size()){       // take care that we don't drop off the end
        if(at[sample] == count) {   // match!
          now= true;
          sample++;
	}
        count++;
      }
    }
    return now;
  }
  operator bool() { return (*this)(); }; // shorthand
  bool operator==(const bool& rhs) { return (*this)() == rhs; };
  bool operator!=(const bool& rhs) { return (*this)() != rhs; };
  bool operator||(const bool& rhs) { return (*this)() || rhs; };
  bool operator&&(const bool& rhs) { return (*this)() && rhs; };
  Sampler(size_t _every) : every(_every) {
    count=0; sample=0;
  }
  Sampler(size_t *_at, size_t number) : at(_at, _at + number) {
    count=0; sample=0;
  }
  void reinit(size_t _every) {    // to reuse the Sampler with different rate
    every=_every;
    count=0; sample=0;
  }
  void reinit(size_t *_at, size_t number) {  // or with different array
    at.resize(number);
    for(size_t i=0; i<number; i++) {
      at[i]=_at[i];
    };
    count=0; sample=0;
  }
  void reset() {    // to reset the Sampler
    count=0; sample=0;
  }
};

/**
 * the general observable, that calculates some value from two physical things,
 * e.g. Wavefunction and LaserField. It allows the second thing to be omitted
 * A specialization observing just one physical thing will be given below
 */
template<class R,class T1, class T2 = void>
class Observable{
  string label;
  R (*fptr)( T1 &, T2& );    /**<pointer to measuring function*/
  T1* phys_obj1_ptr;    /**<pointer to first physical thing*/
  T2* phys_obj2_ptr;    /**<pointer to second physical thing*/
  Output<R> o;    /**<The Output Device*/
  Sampler pnow, rnow;
public:
  Observable(string _label, unsigned int _every,  R (*_fptr)( T1 &, T2& ), T1& po1, T2& po2) :
    label(_label), fptr(_fptr), phys_obj1_ptr(&po1), phys_obj2_ptr(&po2), o(_label), pnow(_every), rnow(_every) {}
  R measure() {
     // calculate the observable
     return (*fptr)(*phys_obj1_ptr, *phys_obj2_ptr);
  }
  void record(double time) {
    if(rnow()) o.save(time, measure());
  }
  void print(double time) {
    if(pnow()) cout << label<< ": "<< time << "\t" << measure() << endl;
  }
};

/**
 * Specialization of the general Observable, measuring just one physical thing
 */
template<class R, class T>
class Observable<R, T, void>{
  string label;
  R last;
  R (* fptr)( T &);    /**<pointer to measuring function*/
  T* phys_obj_ptr;    /**<pointer to physical thing*/
  Output<R> o;    /**<The Output Device*/
  Sampler pnow, rnow;
public:
  Observable(string _label, size_t _every,  R (* _fptr)(T &), T& po ) :
    label(_label), fptr(_fptr), phys_obj_ptr(&po), o(_label), pnow(_every), rnow(_every) { }
  // constructor with array-sampler
  Observable(string _label, size_t *_at, size_t num, R (* _fptr)(T &), T& po ) :
    label(_label), fptr(_fptr), phys_obj_ptr(&po), o(_label), pnow(_at, num), rnow(_at, num) { }
  R measure() {
     // calculate the observable
     return (*fptr)(*phys_obj_ptr);
  }
  void record(double time) {
    if(rnow()) o.save(time, measure());
  }
  void print(double time) {
    if(pnow()) cout << label<< ": "<< time << "\t" << measure() << endl;
  }
};

/**
 * Specialization of the general Observable, measuring on two physical things
 * and outputting directly to a file (good for matlab output, see
 * MatPrj_Wavefunction_X1 in Cylindrical3D_mat.cpp for example)
 */
template<class T1, class T2>
class Observable<void, T1, T2>{
  string label;
  Sampler now;
  size_t sample;
  void (*fptr)( T1&, T2&, string, string );    /**<pointer to measuring function*/
  T1* phys_obj1_ptr;    /**<pointer to first physical thing*/
  T2* phys_obj2_ptr;    /**<pointer to second physical thing*/
public:
  Observable(string _label, size_t _every,  void (*_fptr)( T1 &, T2&, string, string ), T1& po1, T2& po2) :
    label(_label), now(_every), fptr(_fptr), phys_obj1_ptr(&po1), phys_obj2_ptr(&po2) {
    sample=0;
  }
  void record() {
    if(now()) {
      // call the output routine
      string filename=outputFolder()+label+".mat";  
      ostringstream s; s << label << "__" <<setw(6)<<setfill('0')<<sample++;
      (*fptr)(*phys_obj1_ptr, *phys_obj2_ptr, filename, s.str());
    }
  }
};

/**
 * Specialization of the general Observable, measuring just one physical thing
 * and outputting directly to a file (good for matlab output, see
 * MatPrj_Wavefunction_X1 in Cylindrical3D_mat.cpp for example)
 */
template<class T>
class Observable<void, T, void>{
  string label;
  Sampler now;
  size_t sample;
  void (*fptr)( T &, string, string );    /**<pointer to measuring function*/
  T* phys_obj_ptr;    /**<pointer to first physical thing*/
public:
  Observable(string _label, size_t _every,  void (*_fptr)( T &, string, string ), T& po) :
    label(_label), now(_every), fptr(_fptr), phys_obj_ptr(&po) {
    sample=0;
  }
  void record() {
    if(now()) {
      // call the output routine
      string filename=outputFolder()+label+".mat";   
      ostringstream s; s << label << "__" <<setw(6)<<setfill('0')<<sample++;
      (*fptr)(*phys_obj_ptr, filename, s.str());
    }
  }
};


#endif /* OBSERVABLE_H */
