// $Id: ParameterMap.h 512 2006-11-21 15:22:16Z goerke $
#include <cstdlib>
#ifndef PARAMETER_H
#define PARAMETER_H
#include <string>
#include <sstream>
#include <map>
#include <iostream>
#include <vector>
#include "ConfigFileParser.h"
#include <typeinfo>
#include <complex>

using namespace std;

namespace ParameterHelper{

/**
 * Exceptions for the FromStringTo<T> and FromStringToVector<T> conversion
 * functions
 */
  struct Error {};              /**< general error class */
  struct NullError : public Error {};
  struct ConversionError : public Error {
    string msg;
    ConversionError(const string & str = "") : msg(str) {}
  };
  struct SizeChangeWarning : public Error {
    string msg;
    SizeChangeWarning(const string & str = "") : msg(str) {}
  };

/* ----------------------- Utility functions ---------------------- */

/** 
 * \class FromStringTo
 * try to convert a string to a value of type T 
 * \remarks Example : int i = FromStringTo<int>(s)
 * @param str string to convert 
 * @return  value of type T corresponding to the conversion of s to this type
 */
template<typename T> 
inline T FromStringTo(const std::string & str) {
    istringstream is(str);
    T t;
    is >> t;

    bool failFlag=0;
    if(is.fail()) failFlag=1;

    // check if every character of the string has been read
    string buf;
    is >> buf;
    if(failFlag || buf.size() !=0) {
      string type=typeid(t).name();
      throw ConversionError(
        "\"" + str + "\" has a bad format for conversion to type " + type \
        +" in FromStringTo<" + type + ">");
    }

    return t;
}
template<> 
inline string FromStringTo<string>(const std::string & str) {
    return str;
}
/** specialization for bool
 * \li "false", "F", "no", "n", "0" are interpreted as false
 * \li "true", "T", "yes", "y", "1", "-1", or anything else are interpreted as true
*/
template<> 
inline bool FromStringTo<bool>(const std::string & str) {
  /*
  static const char whitespace[] = " \n\t\v\r\f";
  string::size_type i=str.find_first_not_of(whitespace);
  string::size_type f=str.find_last_not_of(whitespace) + 1U;
  string tmp=str.substr(i,f-i);  // trim
  */

  // capitalisation to faciliate comparisons
  string sub=str;
  for(string::iterator it = sub.begin(); it != sub.end(); it++) *it = toupper(*it);
  if( sub==string("FALSE") || sub==string("F") || sub==string("NO") || 
    sub==string("N") || sub==string("0") || sub==string("NONE") )
    return false;
  if( sub==string("TRUE") || sub==string("T") || sub==string("YES") || 
    sub==string("Y") || sub==string("1") )
    return true;
  throw ConversionError(
      "\"" + str + "\" in not a boolean expression");
}

/** 
 * \class FromStringToVector
 * try to convert a string to a vector containing type T 
 * \remarks Example : vector<int> i = FromStringToVector<int>(s)
 * @param s string to convert 
 * @return  vector container of type T corresponding to the conversion of s to this type
 */
template<typename T> 
class FromStringToVector : public vector<T> {
public:
  FromStringToVector(const std::string & str) {
  string::size_type loc = str.find( "(", 0 );
  if( loc == string::npos ) {
    throw ConversionError("opening bracket of vector missing");
  }
  loc+=1;
  //cout<<"loc="<<loc<<endl;
  string::size_type _end = str.find( ")", loc );
  if( _end == string::npos ) {
    throw ConversionError("closing bracket of vector missing");
  }
  //cout<<"end="<<_end<<endl;
  string::size_type last=loc;
  string::size_type index = str.find( ",", loc );
  try {
  while(index != string::npos){
    //cout<<"last="<<last<<endl;
    //cout<<"index="<<index<<endl;
      push_back(FromStringTo<T>(str.substr(last, index-last)));
    //cout<<back()<<endl;
    last=index+1;
    index = str.find( ",", last );
  }
  push_back(FromStringTo<T>(str.substr(last, _end-last)));
  } catch (ConversionError & e) {
    ostringstream s; s << this->size();
    throw ConversionError("[" + s.str() + "]: " + e.msg);
  }		
  }
};

/** 
 * \class ToStringFrom
 * convert any value v of type T to string 
 * \remarks Example: string s = ToStringFrom(14)
 * @param v value to convert
 * @return  a string corresponding to v
 */
template<typename T> 
inline string ToStringFrom(const T & v) {
    ostringstream s;
    s << v;
    return s.str();
}
}	/* </namespace ParameterHelper> */

/* ======================= ParameterMap =========================== */

class Parameter_ABC;  // forward declaration
/**
 * ParameterMap is the registry for all Parameter Objects
 * this is coded as a Wrapper object containing a static map.
 * Since the internal map is static, all ParameterMap instances refer to the
 * same one. i.e. it is guaranteed that only one map exists, but several
 * ParameterMap objects can be instanciated locally in the code to access it
 */
class ParameterMap {
  // Write configuration
    friend void printParameterMap( std::ostream& os, const ParameterMap & f, string _prefix);//="");
private:
  static map<string,Parameter_ABC *> _contentMap; /**< the mapping (key, value) */ 
  static string _delimiter;	  /**< separator between key and value          */
  static string _comment;	  /**< character(s) used to introduce a comment */

public:
  /*  wrapper functions to access methods of contained _contentMap */
  // Return a string value associated to key (eventually throw KeyNotFound)
  //const Parameter_ABC * operator[](const string & key) const;
  Parameter_ABC* operator[](const string & key);
  /// Redefinition of iterators types to faciliate navigation in _contents
  typedef map<string,Parameter_ABC *>::iterator       mapi;
  typedef map<string,Parameter_ABC *>::const_iterator mapci;
  mapci find(const string& str) const { return _contentMap.find(str); }
  mapi find(const string& str) { return _contentMap.find(str); }
  /// Return iterators to the beginning and the end of _contentMap
  mapci begin() const { return _contentMap.begin(); }
  mapci end()   const { return _contentMap.end();   }
  void add(string &label, Parameter_ABC* ptr) {
    if ( _contentMap.count(label) ) {
      cerr<<"Error: parameter with label \""<<label<<"\" doubly defined"<<endl;
      throw ParameterHelper::Error();
    } else {
      /// register in static _contentMap
      _contentMap[label]=ptr;
    }
  }
  void remove(string &label) {
    _contentMap.erase(label);
  }

  /// Update the configuration file with the content of the map
  void saveToFile(const string & filename);
  // late initialization of all map entries
  void init(ConfigFileParser& _cfg);
  /// Exception types when the file can't be created
  struct FileNotCreated {
    string filename;
    FileNotCreated(const string & f = "") : filename(f) {}
  };
  /// Exception types when a key is not found
  struct KeyNotFound {
    KeyNotFound(const string & k = "") {
      cerr<<"No Parameter with label \""<<k<<"\" found! Please make sure to declare it."<<endl;
    }
  };
};

/* ------------------------- Parameter_ABC ------------------------ */
/* -------------- Abstract Base Class for Parameter<T> ------------ */

template<typename T> class Parameter;  // forward declaration

/**
 * Abstract Base Class(ABC) for the Parameter\<T\> class,
 * necessary to be able to store the Parameters of different type
 * in the ParameterMap.
 *
 * @see Parameter\<T\>
 */
class Parameter_ABC { 
  ParameterMap pm;		/**<Access Object for static ParameterMap */

  //Parameter_ABC(const Parameter_ABC&);          /**<hidden, never used */
  //Parameter_ABC operator=(const Parameter_ABC&);/**<hidden, never used */
protected: // the derived classes can access these things
  string label, comment;
public:
  const type_info* myType;
  Parameter_ABC(const type_info* _t, string _label, string _comment);
  virtual ~Parameter_ABC();
  string getLabel() const { return label; };
  string getComment() const { return comment; };
  virtual void setValue(const string & ) = 0;
  virtual void init(ConfigFileParser& _cfg) = 0;
  virtual string getValueString() const = 0;
  class Incompatible_Type_Exception {
    public:
          Incompatible_Type_Exception(const string& l, const string& f, const string& t) {
            cerr<<"Cannot convert Parameter<"<<f<<"> (labeled \""<<l<<"\") to "<<t<<endl;
          }
  };
  template<class T> operator T() {
    if ((* myType) == typeid(T)) {
      return (*static_cast<Parameter<T>*>(this))();
    } else throw(Incompatible_Type_Exception(label, myType->name(), typeid(T).name()));
  }
};

/* ======================= Parameter<T> =========================== */

/**
 * \class Parameter
 */
template<typename T>
class Parameter : public Parameter_ABC{ 
  T _data;
public:
  Parameter(string _label, T _value, string _comment = "");
  Parameter(ConfigFileParser& _cfg, string _label, T _value, string _comment = "");
  virtual void init(ConfigFileParser& _cfg);
  /** 
   * Set the _data from a string
   * @param str  the string representation of the value
   * @return       nothing
   */
  virtual void setValue(const string & str );
  /** 
   * Get a string from the _data
   * @return       the string representation of the value
   */
  virtual string getValueString() const;
  /// \name Access operators
  //@{ 
  const T operator()() const { return _data; };
  T& operator()() { return _data; };
  operator T() const { return _data; };
  operator T&() { return _data; };
  //@} 
  /// \name Assignment operators
  //@{ 
  void operator()(T _value){ _data=_value; };
  void set(T _value){ _data=_value; };
  //@} 
};

/* =================== Parameter<vector<T> > ====================== */

/**
 * specialized Parameter\<vector\<T\> \> class
 */
template<typename T>
class Parameter<vector<T> > : public Parameter_ABC{ 
  vector<T> _data;
  static bool sizeChangeException;
  static bool sizeChangeMessage;
public:
  Parameter(string _label, vector<T> _value, string _comment = "");
  Parameter(ConfigFileParser& _cfg, string _label, vector<T> _value, string _comment = "");
  virtual void init(ConfigFileParser& _cfg);
  /** 
   * Set the _data from a string
   * @param str  the string representation of the vector of values
   * @return       nothing
   */
  virtual void setValue(const string & str );
  /** 
   * Get a string from the _data
   * @return       the string representation of the vector of values
   */
  virtual string getValueString() const;
  /// \name Access operators
  //@{ 
  const vector<T>& operator()() const { return _data; };
  T operator[](const int i) const { return _data[i]; };
  T& operator[](const int i) { return _data[i]; };
  //@} 
  /// \name Assignment operators
  //@{ 
  void operator()(vector<T> _value){ _data=_value; };
  void set(vector<T> _value){ _data=_value; };
  //@} 
  size_t size() const { return _data.size(); };
};
template<typename T> bool Parameter<vector<T> >::sizeChangeException=false;
template<typename T> bool Parameter<vector<T> >::sizeChangeMessage=true;

/* ================== output operators ============================ */

ofstream& operator<<( std::ofstream& , const ParameterMap &);
ostream& operator<<( std::ostream& , const ParameterMap &);

ostream& operator<<( ostream& , const Parameter_ABC& );
//template<typename T> ostream& operator<<( ostream& os, const Parameter<T>& rh ) { return os << rh(); };
    
/* ------------------ general conveniences ------------------------ */
/**
 * output stream operator for vector\<T\> datatypes
 */
template<typename T>
inline ostream& operator<<( ostream& os, const vector<T>& rh ) {
  os << "(";
  typename vector<T>::const_iterator it;
  for (it=rh.begin(); it != rh.end() - 1 ; ++it){
    os << *it << ", ";
  }
  os << *it << ")";
  return os;
};

/**
 * for "const char*"+string this helps
 */
inline string operator+(const char *lhs, const string rhs) { string t(lhs); return t+=rhs; };


/* ======================= Implementations ======================== */

/**
 * Parameter\<T\> class constructors
 */
template<typename T>
Parameter<T>::Parameter(string _label, T _value, string _comment) :
  Parameter_ABC(&typeid(T), _label, _comment), _data(_value) {
  };
template<typename T>
Parameter<T>::Parameter(ConfigFileParser& _cfg, string _label, T _value, string _comment) :
  Parameter_ABC(&typeid(T), _label, _comment), _data(_value) {
  ConfigFileParser::mapci it=_cfg.find(_label);
  if ( _cfg.end() != it ) {
    setValue(it->second);
  } else {
    cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/**
 * Parameter\<T\> init(ConfigFileParser) for late initialization
 */
template<typename T>
void Parameter<T>::init(ConfigFileParser& _cfg) {
  ConfigFileParser::mapci it=_cfg.find(label);
  if ( _cfg.end() != it ) {
    setValue(it->second);
  } else {
    // no message for late initialization
    //cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/** 
 * Set the _data from a string
 * @param str  the string representation of the _data
 * @return       nothing
 */
template<typename T>
inline void Parameter<T>::setValue(const string & str ) {
  try {
    _data = ParameterHelper::FromStringTo<T>( str );
  } catch (ParameterHelper::ConversionError & e) {
    cerr<<"Error: "<<e.msg<<endl;
    throw ParameterHelper::ConversionError();
  }		
};
/** 
 * Get a string from the _data
 * @return       the string representation of the _data
 */
template<typename T>
inline string Parameter<T>::getValueString() const {
  return ToStringFrom<T>(_data);
};

/**
 * Parameter\<vector\<T\> \> class constructors
 */
template<typename T>
Parameter<vector<T> >::Parameter(string _label, vector<T> _value, string _comment) :
  Parameter_ABC(&typeid(T), _label, _comment), _data(_value) {};
template<typename T>
Parameter<vector<T> >::Parameter(ConfigFileParser& _cfg, string _label, vector<T> _value, string _comment) :
  Parameter_ABC(&typeid(T), _label, _comment), _data(_value) {
  ConfigFileParser::mapci it=_cfg.find(_label);
  if ( _cfg.end() != it ) {
    sizeChangeException=true;
    try {
      setValue(it->second);
    } catch (ParameterHelper::SizeChangeWarning & e) {
      cerr<<"Warning: Reading config file "<<_cfg.getFilename()<<": "<<endl;
      cerr<<"\t"<<e.msg<<endl;
    }
    sizeChangeException=false;
  } else {
    cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/**
 * Parameter\<vector\<T\>\> init(ConfigFileParser) for late initialization
 */
template<typename T>
void Parameter<vector<T> >::init(ConfigFileParser& _cfg) {
  ConfigFileParser::mapci it=_cfg.find(label);
  if ( _cfg.end() != it ) {
    sizeChangeException=true;
    try {
      setValue(it->second);
    } catch (ParameterHelper::SizeChangeWarning & e) {
      cerr<<"Warning: Reading config file "<<_cfg.getFilename()<<": "<<endl;
      cerr<<"\t"<<e.msg<<endl;
    }
    sizeChangeException=false;
  } else {
    // no message for late initialization
    //cerr << "Warning: Parameter \"" << _label << "\" not found in config file " << _cfg.getFilename() <<endl;
  }
};

/** 
 * Set the _data from a string
 * @param str  the string representation of the vector of values
 * @return       nothing
 */
template<typename T>
inline void Parameter<vector<T> >::setValue(const string & str) {
  try {
    size_t vs=_data.size();
    _data=ParameterHelper::FromStringToVector<T>( str );
    if (_data.size() != vs) {
      ostringstream size1; size1 << vs;
      ostringstream size2; size2 << _data.size();
      string warn("size of vector \""+getLabel()+"\" changed from "+size1.str()+" to "+size2.str());
      if(sizeChangeException) {
	throw ParameterHelper::SizeChangeWarning(warn);
      } else if(sizeChangeMessage)
        cerr<<"Warning: "<<warn<<endl;
    }
  } catch (ParameterHelper::ConversionError & e) {
    cerr<<"Error in ParameterHelper::FromStringToVector: "<<getLabel()<<e.msg<<endl;
    throw /*ParameterHelper::ConversionError()*/;
  }		
};
/** 
 * Get a string from the _data
 * @return       the string representation of the vector of values
 */
template<typename T>
inline string Parameter<vector<T> >::getValueString() const {
  return ParameterHelper::ToStringFrom<vector<T> >(_data);
};

#endif /* PARAMETER_H */
