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
    static void readDipoleInput(std::ifstream &input);
    static void readOmegaInput(std::ifstream &input);
    static void populateAlphaOne(double real, double imaginary);
    std::vector< std::complex<double> > alphaOne;
    std::vector< std::complex<double> > alphaThree;




	std::vector<double> omega;
	std::vector<int> m;
	std::vector<double> E_m;
	double E_0;
	double alpha;	
	double T;
	std::string OutputFolder;
	std::string dataFolder;
    static std::ifstream omegaDip;
    std::ifstream dipOne;
    std::ifstream dipThree;
    static std::vector<double> reAlphaOne;
    static std::vector<double> imAlphaOne;
private:

};


#endif /* pert_h */
