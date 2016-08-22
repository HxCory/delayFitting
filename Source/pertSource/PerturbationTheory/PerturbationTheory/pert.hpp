//
//  pert.hpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#ifndef pert_h
#define pert_h

#include <iostream>
#include <vector>
#include <fstream>
#include <istream>
#include <complex>

class pert
{
public:
	pert();
	virtual ~pert();
    static void openFile(std::ifstream &filename, std::string name);
	static void closeFile(std::ifstream &filename);
    static void readInput(std::ifstream &input);
    std::vector<double> a1;


private:
	
};


#endif   /* pert_h */
