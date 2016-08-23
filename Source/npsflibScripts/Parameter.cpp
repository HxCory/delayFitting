#include <string>
#include <cstdlib>
#include <map>
#include <iostream>
#include "ParameterMap.h"
#include <iomanip>

using namespace std;

/* ======================= static data ============================ */

map<string,Parameter_ABC *> ParameterMap::_contentMap; /**< the static core of the ParameterMap Class*/ 
string ParameterMap::_delimiter="="; // defaults
string ParameterMap::_comment="#";   // defaults

/* ----------------------- Friends -------------------------------- */

inline void printParameterMap( std::ostream& os, const ParameterMap & f, string _prefix) {
  const int offset = 20;
  os << _prefix << f._comment << " --------------------- <Configuration> ---------------------" << endl;
  // print date of run and directory
  for (ParameterMap::mapci it = f.begin(); it != f.end(); it++) {
    Parameter_ABC &p = (* it->second);   // short
    os << _prefix << it->first << " " << f._delimiter << " " << left << setw(offset) << p.getValueString();
    string comment = p.getComment();
    if (comment.size() != 0)
      os << " " << f._comment << " " << comment;
    os << endl;

  }
  os << _prefix << f._comment << " --------------------- </Configuration> --------------------" << endl;
}

/* ================== output operators ============================ */

/** 
 * This operator is used to save (redirect) a ParameterMap to (out)file stream
 * \todo : add comment management
 */
ofstream& operator<<( std::ofstream& os, const ParameterMap & f) {
  printParameterMap(os, f, "# ");
  return os;
}

/** 
 * This operator is used to save (redirect) a ParameterMap to output stream 
 * \todo : add comment management
 */
ostream& operator<<( std::ostream& os, const ParameterMap & f) {
  printParameterMap(os, f);
  return os;
}

/**
 * output stream operator for Parameter_ABC (Abstract Base Class)
 */
ostream& operator<<( ostream& os, const Parameter_ABC& rh ) {
  //return os << rh.getLabel() << " = " << rh.getValueString();
  return os << rh.getLabel() << " = " << rh.getValueString() << "\t\t# " << rh.getComment();
};

/* ======================= Implementations ======================== */

Parameter_ABC::Parameter_ABC(const type_info* _t, string _label, string _comment) :
 label(_label), comment(_comment), myType(_t) {
    pm.add(label, this);
};
Parameter_ABC::~Parameter_ABC(){		// plain virtual destructor is important for derived classes. This way delete(Parameter_ABC* x) calls the destructor of the derived class that x belongs to. Yet this routine is called implicitely after the call to the destructor of the derived class.
  //cerr<<"~Parameter_ABC("<<label<<")"<<endl;
  pm.remove(label);
};

void ParameterMap::saveToFile(const string & filename)
{
    ofstream cfgFile(filename.c_str(), ios::out);
    if ( !cfgFile) throw FileNotCreated(filename);	
    cfgFile << (*this);
    cfgFile.close();	
}
void ParameterMap::init(ConfigFileParser& _cfg)
{
  for (mapci it = begin(); it != end(); it++) {
    Parameter_ABC &p = (* it->second);   // short
    p.init(_cfg);
  }
}

/** 
 * Faciliate access to a value associated to a key. Throw KeyNotFound when the key 
 * doesn't exist.
 * \note
 * Difference to map\<string,string\>::operator[] is the behaviour on KeyNotFound
 * @param key    the name of the key
 * @return       the value associated (may be I shall return an empty string if 
 *               the key does not exist). 
 */
/*
const Parameter_ABC* ParameterMap::operator[](const string & key) const {
    mapci it=find(key);
    if ( end()==it ) throw KeyNotFound(key); 
    return it->second; 

}
*/
Parameter_ABC* ParameterMap::operator[](const string & key) {
    mapi it=find(key);
    if ( end()==it ) throw KeyNotFound(key); 
    return it->second; 

}

