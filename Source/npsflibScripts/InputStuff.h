#include "ConfigFileParser.h"
#include <cstdlib>
#include "ParameterMap.h"
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include "Observable.h"

extern ConfigFileParser cfg; // global

inline std::string create_dir(std::string name);

void touchdir(std::string name);

void parseInput(int argc, char *argv[]);

