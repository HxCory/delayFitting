/**
 * @file   ConfigFileParser.cpp
 * @author Sebastien Varrette <Sebastien.Varrette@imag.fr>
 * @date   Fri Dec  9 2005   
 *
 * Copyright (c) 2005 Sebastien Varrette  (http://www-id.imag.fr/~svarrett/)
 * 
 * @brief  see ConfigFileParser.h 
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


#include <cstdlib>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;

#include "ConfigFileParser.h"

/// This constant is used to parse a line for guillemets (I ignore the english word) 
static const char guillemet[] = { '\"', '\'' };

/** 
 * Default constructor for class ConfigFileParser, getting keys and values from a 
 * given file
 * @param filename     the name of the configuration file
 * @param delimiter    separator between key and value
 * @param comment      character(s) used to introduce a comment
 */
ConfigFileParser::ConfigFileParser(const string & filename, const string & delimiter,
				   const string & comment) :
    _filename(filename), _delimiter(delimiter), _comment(comment)
{
    ifstream cfgFile(_filename.c_str());
    if (!cfgFile) throw FileNotFound(_filename);
    cfgFile >> (*this); // here I use operator >> to load the configuration file
}

/** 
 * read keys and values from a given file
 * @param filename     the name of the configuration file
 */
void ConfigFileParser::read(const string & filename) {
    _filename=filename;
    ifstream cfgFile(_filename.c_str());
    if (!cfgFile) throw FileNotFound(_filename);
    cfgFile >> (*this); // here I use operator >> to load the configuration file
}

/** 
 * This static function serialize the string s by removing leading and trailing 
 * whitespace
 * @param s 
 */
void ConfigFileParser::_serialize(string & s) {
    static const char whitespace[] = " \n\t\v\r\f";
    s.erase(0, s.find_first_not_of(whitespace));
    s.erase(s.find_last_not_of(whitespace) + 1U);
}


/** 
 * Check if the string s contains a comment of not. The commented line is returned 
 * in cs. 
 * @param s    the string to check 
 * @param cs   the commented string
 * @return     'true' if s contains a comment 
 */
bool ConfigFileParser::_containsComment(const string & s, string & cs) const {
    return ((cs = _getCommentFrom(s)) != "");
}

/** 
 * extract a comment for a string s (return "" if s contains no comment)
 * @param s   the string to check
 * @return    the comment contained in s (or emptystring if no comment is contained
 */
string ConfigFileParser::_getCommentFrom(const string & s) const {
    string tmp = s;
    return tmp.erase(0, s.find(_comment));
}


/** 
 * check if the string s contains a key and an associated value 
 * @param s        the string to check
 * @param key      the key 
 * @param value    the associated value
 * @param cs       a comment that could be in the string s
 * @param valueWithGuillemet  'true' if the value is enclosed by " or '
 * @return         'true' if s contains a key and a value
 */
bool ConfigFileParser::_checkKeyPresence(const string &s, string &key, string &value,
				   string & cs, bool & valueWithGuillemet) const {
    cs = _getCommentFrom(s); // get eventually the comment contained in the line
    valueWithGuillemet = false;
    string tmp = s; 
    tmp = tmp.substr(0, s.find(_comment)); // extract the comment
    string::size_type delimPos = tmp.find(_delimiter);
    if (delimPos < string::npos) { 
	key = tmp.substr(0, delimPos); // extract the key 
	_serialize(key);
	// Extract value
	value = tmp.substr(delimPos+1, tmp.find(_comment)-delimPos-1);
	_serialize(value);
	// Now check for guillemets 
	if ((value.find("\"") < string::npos) || (value.find("'") < string::npos)) {
	    value.erase(0, value.find_first_not_of(guillemet));
	    value.erase(value.find_last_not_of(guillemet) + 1U);
	    valueWithGuillemet = true;
	}
	return true;
    } else {
	key = "";
	value = "";
	return false;
    }
}



/** 
 * This operator is used to save (redirect) a configuration file to os 
 * TODO : add comment management
 */
ostream& operator<<( std::ostream& os, const ConfigFileParser & f) {
    ifstream cfgFile(f._filename.c_str());
    if ( cfgFile.is_open() ) {
	string line;
	// This map wil be used to check if added key should be printed 
	map<string,string> tmpMap = f._contentMap;
	while ( !cfgFile.eof() ) {
	    getline (cfgFile,line);
	    string key="", value="", comment="";
	    bool hasGuillemet = false;
	    // Check if the line contains a key
	    bool checkKey = f._checkKeyPresence(line,key,value,comment,hasGuillemet);
	    if ( checkKey ) {
		// check if the value has changed
		try {
		    if ( f.getValue<string>(key) != value ) {
			string g = (hasGuillemet?"\"":"");
			os << key << "\t" << f._delimiter << " " << 
			    g << f.getValue<string>(key) << g;
			if (comment != "") os << "\t" << comment;
		    } else os << line;
		    // Here, we're sure the key exists in the map
		    tmpMap.erase( tmpMap.find(key) ); 
		    os << endl;
		} catch (ConfigFileParser::KeyNotFound & e) {
		    // This is a unknown key. Either it should be added to 
		    // _contentMap[] or it is a deleted key. In all case, 
		    // I shall ignore it.
		}		
	    } else {
		os << line << endl;
	    }
	}
        cfgFile.close();

	// Now proceed with keys that could have been added and that are not in the 
	// configuration file
	for (ConfigFileParser::mapci it=tmpMap.begin(); it!=tmpMap.end(); it++) {
	    os << it->first << "\t" << f._delimiter << " " << it->second << endl;
	}	

    } else {
	// Il the file doesn't exist, create a new one
	os << f._comment << "--------------------" << endl;
	os << f._comment << " Configuration File " << endl;
	os << f._comment << "--------------------" << endl;	
	for (ConfigFileParser::mapci it = f.begin(); it != f.end(); it++) {
	    os << it->first << " " << f._delimiter << " " << it->second << endl;
	}
    }
    return os;
}

/** 
 * This operator is used to load (redirect) a configuration file from is 
 */
istream& operator>>( std::istream& is, ConfigFileParser & f) {
    const string::size_type skip = f._delimiter.length();   // length of separator
    
    string nextline = "";  // might need to read ahead to see where value ends
	
    while( is || nextline.length() > 0 ) {
	// Read an entire line at a time
	string line;
	if( nextline.length() > 0 ) {
	    line = nextline;  // we read ahead; use it now
	    nextline = "";
	} else getline(is, line );
	
	// Ignore comments
	line = line.substr(0, line.find(f._comment) );
	
	// Parse the line if it contains a delimiter
	string::size_type delimPos = line.find(f._delimiter);
	if (delimPos < string::npos) { 
	    // Reminder: the special value string::npos represents the maximum 
	    // number of characters there can be in a string
	    
	    // Extract the key
	    string key = line.substr( 0, delimPos );
	    line.replace( 0, delimPos+skip, "" );
	    
	    // See if value continues on the next line
	    // Stop at blank line, next line with a key, or end of stream
	    bool terminate = false;
	    while( !terminate && is ) {
		getline( is, nextline );
		terminate = true;
		
		string nlcopy = nextline;
		ConfigFileParser::_serialize(nlcopy);
		if( nlcopy == "" ) continue;
				
		nextline = nextline.substr( 0, nextline.find(f._comment) );
		if( nextline.find(f._delimiter) != string::npos ) continue;
		nlcopy = nextline;
		ConfigFileParser::_serialize(nlcopy);
		if( nlcopy != "" ) line += "\n";
		line += nextline;
		terminate = false;
	    }
	    // Check for guillemets arround the value
	    line.erase(0, line.find_first_not_of(guillemet));
	    line.erase(line.find_last_not_of(guillemet) + 1U);
	    // Store key and value
	    f.add<string>(key, line);
	}
    }
    return is;    
}

/** 
 * Faciliate access to a value associated to a key. Throw KeyNotFound when the key 
 * doesn't exist. The second form enable the writting 
 * \cfg[i]
 * @param key    the name of the key
 * @return       the value associated (may be I shall return an empty string if 
 *               the key does not exist). 
 */
string ConfigFileParser::operator[](const string & key) const {
    mapci i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    return i->second; 

}
string & ConfigFileParser::operator[](const string & key) {
    mapi i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    return i->second; 

}

/** 
 * Print map content
 * @param os 
 */
void ConfigFileParser::printMapContent(ostream & os) {
    int offset = 20;
    os << setfill ('-') << setw (2*offset + 3) << "" << setfill (' ') << endl;
    os << left << setw(offset) << "     Key    " << "  " << "    Value" << endl;
    os << setfill ('-') << setw (2*offset + 3) << "" << setfill (' ') << endl;;
    for (mapci it = begin(); it != end(); it++) {
	os << left << setw(offset) << it->first << "| " << it->second << endl;
    }
}
/** 
 * Print directly the content of the configuration file.
 * @param os 
 */
void ConfigFileParser::printCurrentFileContent(ostream & os) {
    ifstream cfgFile(_filename.c_str());
    string line;
    if (! cfgFile ) throw FileNotFound(_filename); 
    while ( !cfgFile.eof() ) {
	getline (cfgFile, line);
	os << line << endl;
    }
}


/** 
 * Remove an existing key. Throw KeyNotFound if the key doesn't exist.
 * @param key     the key to remove
 */void ConfigFileParser::remove(const string& key) {
    mapi i=_contentMap.find(key);
    if ( end() == i ) throw KeyNotFound(key); 
    _contentMap.erase( i ); 
}

/** 
 * Update the configuration file named filename with the content of the map
 * @param filename    the name of the file to update
 */
void ConfigFileParser::saveCfgFile(const string & filename) { 
    string f = ((filename == "")? _filename : filename);
    ofstream cfgFile(f.c_str(), ios::out);
    if ( !cfgFile) throw FileNotCreated(f);	
    cfgFile << (*this);
    cfgFile.close();	
}
 
