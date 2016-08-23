/**
 * @file   ConfigFileParser.h
 * @author Sebastien Varrette <Sebastien.Varrette@imag.fr>
 * @date   Fri Dec  9 2005  
 *
 * Copyright (c) 2005 Sebastien Varrette  (http://www-id.imag.fr/~svarrett/)
 *
 * @brief  This library provides a mechanism to parse and update a configuration 
 *         file. It follows my experiments with confloader of Nicolas Bernard 
 *         (http://www.lafraze.net/nbernard/) and is partly inspired by the work of
 *         Jeff Schiller(CodeDread) http://www.codedread.com/code/ConfigParser/docs/
 *         (Yet, this work appears to be only available for windows users) and Rick 
 *         Wagner (http://www-personal.engin.umich.edu/~wagnerr/ConfigFile.html).
 *         This last work was really interesting yet I decided to improve this code 
 *         for a better management of file writing and entry access. 
 *         So finally I wrote my own library using my own programming style :-)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 * 
 * Sebastien Varrette                               \n
 * <Sebastien.Varrette@imag.fr>                     \n
 * University of Luxembourg / ID-IMAG Laboratory    \n
 * 162-A avenue de la faïencerie                    \n
 * L-1511 Luxembourg, LUXEMBOURG                    \n
 */
 /********************************************************************************/
#include <typeinfo>
#ifndef __CONFIGFILEPARSER_H
#define __CONFIGFILEPARSER_H


#include <string>
#include <cstdlib>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

/// These functions are used to faciliate conversion to/from string from/to a type T
template<class T> T      FromStringTo(const string & s);
template<class T> string ToStringFrom(const T & v);

/**
 * This is the main class representing a configuration file and making it possible 
 * to read/write configuration entries. 
 */
class ConfigFileParser {
    // Write or read configuration
    friend ostream& operator<<(ostream& os, const ConfigFileParser & f);
    friend istream& operator>>(istream& is, ConfigFileParser & f);

private:
    string _filename;		  /**< name of the configuration file           */
    string _delimiter;		  /**< separator between key and value          */
    string _comment;		  /**< character(s) used to introduce a comment */
    map<string,string> _contentMap; /**< the mapping (key, value)               */

    static void _serialize(string & s);/**< Remove leading and trailing whitespace */
    /// return 'true' if s contains a comment (which is returned in cs)
    bool _containsComment(const string & s, string & cs) const ;
    /// extract comment from the string line s
    string _getCommentFrom(const string & s) const; 
    /// extract key, value and comment from the string line s
    bool _checkKeyPresence(const string & s, string & key, string & value, 
			   string & cs, bool & valueWithGuillemet) const;
    
public:
    /// Default constructor without a file
    ConfigFileParser() : _delimiter("="), _comment("#") {}
    /// Main constructor with a file name.
    ConfigFileParser(const string & filename, const string & delimiter = "=", 
		     const string & comment = "#");
    ~ConfigFileParser() {}

    /// Read the configuration file (in case you used the default constructor)
    void read(const string & filename);

    /*** Accessors ***/
    string getFilename()          const { return _filename;  }
    string getDelimiter()         const { return _delimiter; }
    string getCommentIntroducer() const { return _comment;   }

    /// Redefinition of iterators types to faciliate navigation in _contents
    typedef map<string,string>::iterator       mapi;
    typedef map<string,string>::const_iterator mapci;
    /// Find key
    mapci find(const string& key) const { return _contentMap.find(key); }
    mapi find(const string& key) { return _contentMap.find(key); }
    // Get value associated to a given key 
    template<class T> T getValue(const string & key )                        const;  
    template<class T> T getValue(const string & key, const T & defaultValue) const;
    template<class T> bool getValueInto(T & var, const string & key)         const;
    template<class T> bool getValueInto(vector<T> & var, const string & key) const;
    template<class T> bool getValueInto(T & var, const string & key, 
					const T& defaultValue) const;
    // Return a string value associated to key (eventually throw KeyNotFound)
    string operator[](const string & key) const;
    string & operator[](const string & key);
    /// Return iterators to the beginning and the end of _contentMap
    mapci begin() const { return _contentMap.begin(); }
    mapci end()   const { return _contentMap.end();   }

    /*** Mutators ***/
    void setFilename(const string & f)          { _filename = f;  }
    void setDelimiter(const string & d)         { _delimiter = d; }
    void setCommentIntroducer(const string & c) { _comment = c;   }
    /// Add a new keys/values
    template<class T> void add(string key, const T & value);
    /// remove a key
    void remove(const string& key);
    /// Set value associated to a given key
    template<class T> void setValue(const string & key, const T & value);
    /// Update the configuration file with the content of the map
    void saveCfgFile(const string & filename = "");

    /*** Print methods ***/
    void printMapContent(ostream & os = cout);         /**< Print map content */
    void printCurrentFileContent(ostream & os = cout); /**< output file content */

    struct Error {};              // general error class
    /// Exception types when the file is not found
    struct FileNotFound : public Error {
	string filename;
	FileNotFound(const string & f = "") : filename(f)
        {
          cerr<<"ConfigFileParser: File \""<<f<<"\" could not be found"<<endl;
        }
    };
    /// Exception types when the file can't be created
    struct FileNotCreated : public Error {
	string filename;
	FileNotCreated(const string & f = "") : filename(f)
        {
          cerr<<"ConfigFileParser: File \""<<f<<"\" could not be created"<<endl;
        }
    };
    /// Exception types when a key is not found
    struct KeyNotFound : public Error {
	string key;
	KeyNotFound(const string & k = "") : key(k) {}
    };
    /**
     * Exceptions for the FromStringTo<T> and FromsStringToVector<T> conversion
     * functions
     */
    struct ConversionError : public Error {
      string msg;
      ConversionError(const string & str = "") : msg(str) {}
    };
};

/********************************************************************************/
/********************************************************************************/
/********************************************************************************/
// Template functions have to be declared in the same file

/** 
 * FromStringTo convert a string to a value of type T 
 * Example : int i = FromStringTo<int>(s)
 * @param s string to convert 
 * @return  value of type T corresponding to the conversion of s to this type
 */
template<class T> 
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
      throw ConfigFileParser::ConversionError(
        "\"" + str + "\" has a bad format for conversion to type " + type \
        +" in FromStringTo<" + type + ">");
    }

    return t;
}
/**
 * specialization for string
 */
template<> inline string FromStringTo<string>(const string & str) { return str; }
/**
 * specialization for bool
 * "false", "F", "no", "n", "0" are interpreted as false
 * "true", "T", "yes", "y", "1", "-1", or anything else are interpreted as true
 */
template<> inline bool FromStringTo<bool>(const string & str) { 
    // capitalisation to faciliate comparisons
    string sup = str;
    for(string::iterator it = sup.begin(); it != sup.end(); it++) *it = toupper(*it);
    if( sup==string("FALSE") || sup==string("F") || sup==string("NO") || 
	sup==string("N") || sup==string("0") || sup==string("NONE") )
	return false;
    if( sup==string("TRUE") || sup==string("T") || sup==string("YES") || 
      sup==string("Y") || sup==string("1") )
      return true;
    throw ConfigFileParser::ConversionError(
        "\"" + str + "\" in not a boolean expression");
}

/** 
 * FromStringToVector try to convert a string to a vector containing type T 
 * Example : vector<int> i = FromStringToVector<int>(s)
 * @param s string to convert 
 * @return  vector container of type T corresponding to the conversion of s to this type
 */
template<class T> 
class FromStringToVector : public vector<T> {
public:
  FromStringToVector(const std::string & str) {
  string::size_type loc = str.find( "(", 0 );
  if( loc == string::npos ) {
    throw ConfigFileParser::ConversionError("opening bracket of vector missing");
  }
  loc+=1;
  //cout<<"loc="<<loc<<endl;
  string::size_type _end = str.find( ")", loc );
  if( _end == string::npos ) {
    throw ConfigFileParser::ConversionError("closing bracket of vector missing");
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
  } catch (ConfigFileParser::ConversionError & e) {
    ostringstream s; s << this->size();
    throw ConfigFileParser::ConversionError("[" + s.str() + "]: " + e.msg);
  }		
  }
};

/** 
 * ToStringFrom convert any value v of type T to string 
 * Exemple: string s = ToStringFrom(14)
 * @param v value to convert
 * @return  a string corresponding to v
 */
template<class T> 
string ToStringFrom(const T & v) {
    ostringstream s;
    s << v;
    return s.str();
}


/** 
 * Get the value corresponding to key. Throw KeyNotFound when the key doesn't exist
 * @param key    the name of the key 
 * @return       the value corresponding
 */
template<class T> 
T ConfigFileParser::getValue(const string & key ) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    try {
      return FromStringTo<T>( i->second );
    } catch (ConfigFileParser::ConversionError & e) {
      cerr<<"Error: "<<e.msg<<endl;
      throw ConfigFileParser::ConversionError();
    }
}

/** 
 * Get the value corresponding to key. 
 * @param key            the name of the key 
 * @param defaultValue   default value to return if the key is not found
 * @return               the value corresponding to the key, or defaultValue if the 
 *                       key doesn't exist.
 */
template<class T> 
T ConfigFileParser::getValue(const string & key, const T & defaultValue) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) return defaultValue;
    try {
      return FromStringTo<T>( i->second );
    } catch (ConfigFileParser::ConversionError & e) {
      cerr<<"Error: "<<e.msg<<endl;
      throw ConfigFileParser::ConversionError();
    }
}

/** 
 * Get the value corresponding to key and store it in var
 * @param var    the variable where the key value will be stored (if the key exists)
 * @param key    the name of the key 
 * @return       'true' if key is found, 'false' otherwise and in that case, var is 
 *               unchanged.
 */
template<class T> 
bool ConfigFileParser::getValueInto(T & var, const string & key) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) return false;
    try {
      var = FromStringTo<T>( i->second );
    } catch (ConfigFileParser::ConversionError & e) {
      cerr<<"Error: "<<e.msg<<endl;
      throw ConfigFileParser::ConversionError();
    }
    return true;
}

/** 
 * Get the value corresponding to key and store it in var
 * @param var    the variable where the key value will be stored (if the key exists)
 * @param key    the name of the key 
 * @return       'true' if key is found, 'false' otherwise and in that case, var is 
 *               unchanged.
 */
template<class T> 
bool ConfigFileParser::getValueInto(vector<T> & var, const string & key) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) return false;
    try {
      var = FromStringToVector<T>( i->second );
    } catch (ConfigFileParser::ConversionError & e) {
      cerr<<"Error: "<<e.msg<<endl;
      throw ConfigFileParser::ConversionError();
    }
    return true;
}

/** 
 * Get the value corresponding to key and store it in var
 * @param var    the variable where the key value will be stored (if the key exists)
 * @param key    the name of the key 
 * @param defaultValue   default value to return if the key is not found
 * @return       'true' if key is found, 'false' otherwise and in that case, var is 
 *               unchanged.
 */
template<class T> 
bool ConfigFileParser::getValueInto(T & var, const string & key, 
				    const T& defaultValue) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) {
      var = defaultValue;
      return false;
    } else {
      try {
        var =  FromStringTo<T>( i->second );
      } catch (ConfigFileParser::ConversionError & e) {
        cerr<<"Error: "<<e.msg<<endl;
        throw ConfigFileParser::ConversionError();
      }
      return true;
    }
}

/** 
 * Set the value associated to key (throw KeyNotFound when the key doesn't exist)
 * @param key      the name of the key 
 * @param value    the associated value
 */
template<class T> 
void ConfigFileParser::setValue(const string & key, const T & value) {
    mapi i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    i->second = ToStringFrom<T>(value);
}
template<> 
inline void ConfigFileParser::setValue<bool>(const string & key, const bool & value) {
    mapi i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    i->second = (value ? "yes":"no");
}



/** 
 * Add a new map (key, value) to _contentMap. If the key already exist, its value is
 * overwritted.
 * @param key      the name of the key 
 * @param value    the associated value
 */
template<class T>
void ConfigFileParser::add(string key, const T & value) {
    string v = ToStringFrom<T>(value);
    _serialize(key); 
    _serialize(v);
    _contentMap[key] = v;
}
template<>
inline void ConfigFileParser::add(string key, const bool & value) {
    string v = (value ? "yes":"no");
    _serialize(key); 
    _contentMap[key] = v;
}

#endif  // __CONFIGFILEPARSER_H



