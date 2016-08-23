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
#include "wavefunction.h"

class pert
{
public:
	pert();
	virtual ~pert();
    static void openFile(std::ifstream &filename, std::string name);
	static void closeFile(std::ifstream &filename);
    static void readInput(std::ifstream &input);
    std::vector<double> a1;
    std::complex<double> secondIntegral(int m, double tPrime, double t, double omega,
    				std::vector<double> field, std::vector<double> energy, double dt);
    std::complex<double> firstIntegral(int m, double t, double omega, std::vector<double> field,
    									std::vector<double> energy, double dt);
    
    //wavefunction operators
    static void Initialize(wavefunction &wf, int nPoint, double spatialStep, int symm = 0);
    
    //Dipole Operators
    std::complex<double> dipole(wavefunction &wf, wavefunction &wf2);
    std::complex<double> dipolePlaneWave(wavefunction &wf, double &k);
    

private:
	double groundEnergy;
};


#endif   /* pert_h */
