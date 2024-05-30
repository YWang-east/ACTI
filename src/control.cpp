/*
 *  control.cpp
 *
 *  Created  by Christian Wolfsteiner
 *  Modified by Patrick Jenny
 *
 *  Copyright Christian Wolfsteiner. All rights reserved.
 *
 */

#include "control.h"
#include <string>
#include <fstream>
#include <iostream>

Control::Control() {};

/**
 * ctor: control object from filename
 */
Control::Control(string file, string path) {
  _file = path + file;                     // record path of control file
  int i = _file.find_last_of('/');
  if(i)
    path = _file.substr(0,i + 1);
  addEntry("ctrl_path", path);
  read();
}

/**
 * dtor: loop over all map entries and delete pointers
 */
Control::~Control() {
  for(map<string,ControlPair>::iterator itm = _ctrl.begin(); itm != _ctrl.end(); itm++) {
    //cout << endl << itm->first << "  " << itm->second._val << "  " << itm->second._ptr;
    if(itm->second._ptr) {
      //itm->second._ptr = NULL;
      //cout << " now" << itm->second._ptr;
    }
  }
}

/**
 * add a new entry to the map
 * pass by value on purpose (so string constants can be used directly)
 * a val == "default" leaves a field key unassigned
 */
void Control::addEntry(string key, string val, void* ptr) {
  if(val == "default") return;
  ControlPair cp;
  cp._val = val;
  cp._ptr = ptr;
  _ctrl[key] = cp;
}

/**
 * add a pointer to the map, create entries on the fly
 * pass by value on purpose (so string constants can be used directly)
 */
void Control::addPointer(string key, void* ptr, string val) {
  addEntry(key,val,ptr);
}

/**
 * read the control file
 */
void Control::read() {
  string key,val,dummy;
  ifstream file(_file.c_str());
  if(!file) cout << "!! error opening file " << _file << endl;
  while(file) {
    file >> key >> val;
    getline(file,dummy,'\n');
    addEntry(key,val + dummy);
  }
  file.close();
}

/**
 * see if a keyword is in the control file
 * @param key keyword as it appears in the control file
 */
bool Control::exists (string key) {
  return _ctrl.count(key);
}

/**
 * returns a token from an item
 */
string Control::getToken(string key, int tok) {
	char token[200],*tokptr;
	strcpy(token,_ctrl[key]._val.c_str());
	tokptr = token;
	if(tok) {
		tokptr = strtok(token," 	");
		for(int i = 0; i < tok - 1; i++)
			tokptr = strtok(NULL," 	");
		if(!tokptr) {
			printf("\n!! Control: tokenizing item %s, no token %d in %s",key.c_str(),tok,_ctrl[key]._val.c_str());
			return string("");
		}
	}
	return string(tokptr);
}

/**
 * return the control parameter as bool
 * if key is not in map, false will be returned
 * @param key keyword as it appears in the control file
 * @param tok if tok> 0: read the tokth token
 */
bool Control::getBool(string key, int tok) {
  if (!_ctrl.count(key)) return false;
  return getInt(key,tok);
}

/**
 * return the control parameter as integer
 * if key is not in map, -1 will be returned
 * @param key keyword as it appears in the control file
 * @param tok if tok> 0: read the tokth token
 */
int Control::getInt(string key, int tok) {
  if (!_ctrl.count(key)) return -1;
  return atoi(getToken(key,tok).c_str());
}

/**
 * return the control parameter as double
 * a unit following the number will be considered
 * if key is not in map, "-1" will be returned
 * @param key keyword as it appears in the control file
 * @param tok if tok> 0: read the tokth token
 */
double Control::getDouble(string key, int tok) {
  if (!_ctrl.count(key)) return -1;
  return readWithUnit(getToken(key,tok).c_str());
}

/**
 * return the control parameter as string
 * if key is not in map, "default" will be returned
 * @param key keyword as it appears in the control file
 * @param tok if tok> 0: read the tokth token
 */
string Control::getString(string key, int tok) {
  if (!_ctrl.count(key)) return string("default");
  return getToken(key,tok);
}

/**
 * reads a "number unit" string and converts it to a
 * single double considering the unit
 * returns -1 if nothing was converted
 */
double Control::readWithUnit(const char *unitstr) {
  double val = 0;
  char c[100];
  string unames [] = {"","s","m","h","d","w","y"};
  double uvals  [] = {1,1,60,3600,86400,86400*7,86400*365.25};
  int i = 6;
  sscanf(unitstr,"%le %s",&val,c);
  while (unames[i] != string(c) && i > 0) i--;
  val *= uvals[i];
  return val;
}

/**
 * read lines in a stream until the key string is found
 * (at the beginning of a line)
 * @param istr the stream to be read from
 * @param key  the key we are looking for
 */
void Control::findKey(istream& istr, string key) {
  char dummy[256];
  while(key != dummy && istr) {
    istr.getline(dummy,256,'\n');
    while(istr.peek() == '\n') istr.ignore(); 							// this took me an hour: get dies if a line starts with a \n
    istr.get(dummy,key.length() + 1,'\n');
  }
  if(istr) istr.get(dummy,256,' '); 												// eat until next space
}

/**
 * attempts to open a the file with name
 * and place it past the optional key
 * (at the beginning of a line)
 * @param file the stream to be opened for reading
 * @param name filepath and name
 * @param key key to be found
 */
bool Control::openFile(ifstream& file, string name, string key) {
  file.open(name.c_str());
  if(!file) {
    printf("\n!! error opening file %s \n",name.c_str());
    exit(0);
  }
  if(key != "")
    findKey(file,key);  
  return bool(file);
}
