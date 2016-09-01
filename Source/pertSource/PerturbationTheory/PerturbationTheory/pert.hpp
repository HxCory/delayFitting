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
    static void openFile(std::ofstream &filename, std::string name);
	static void closeFile(std::ofstream &filename);
    std::vector<double> a1;
    std::complex<double> secondIntegral(int m, double tPrime, double &t,
    							double omega, std::vector<double> &field, 
    							std::vector<double> &energy, double dt);
    std::complex<double> firstIntegral(int m, double &t, double &omega, 
    								std::vector<double> &field,
    								std::vector<double> &energy, double &dt);
    
    void setEnergies(std::vector<double> &energies);
    
    //wavefunction operators
    static void Initialize(wavefunction &wf, int nPoint, double spatialStep,
    														int symm = 0);
    
    //Dipole Operators
    std::complex<double> dipole(wavefunction &wf, wavefunction &wf2);
    std::complex<double> dipolePlaneWave(wavefunction &wf, double k);
    vector<double> stateEnergy;
    double getMomentum(double &omega);

    //Field Operators
    void createCosineSquare(vector<vector<double>> &dummy, vector<double> &dummyTime,
    			double dt, double amp, double length, vector<double> &dummyOmega,
    			double phi);
    void takeEnergy(vector<double> &dummy);

    //d's
	double dt0;
	double domega0;
	double omegaMin;
	double omegaMax;

private:
	double groundEnergy;
    double firstEnergy;
    double secondEnergy;
    double thirdEnergy;
    double fourthEnergy;
    double fifthEnergy;

};


#endif   /* pert_h */
