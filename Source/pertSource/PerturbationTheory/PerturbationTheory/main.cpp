//
//  main.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#include "pert.hpp"
#include "wavefunction.h"
#include "ParameterMap.h"
#include "constant.h"
#include "ConfigFileParser.h"

using namespace std;

Parameter<string> wavefunctionFolder("wavefunctionFolder",
 "/Users/cgoldsmith/repos/delayFitting/Data/eigenstates/WF", "directory location");
Parameter<string> outputFolder("outputFolder", "/Users/cgoldsmith/repos/delayFitting/Data/pertOutput");
Parameter<int> nState("nState", 4, "number of states");
Parameter<int> nPoint("nPoint", 8200, "number grid points");
Parameter<double> spatialStep("spatialStep", 0.0732, "dx");


int main(int argc, const char* argv[]) {
    
    pert pObject;
	vector<double> omega;
    
    vector<int> m;
	vector<double> E_m;
    double E0 = sqrt(1e+14/3.5101e+16);
	double T = 3/0.02419;
	
    vector<double> fieldX;
    vector< complex<double> > dip1f;
    vector< complex<double> > dip3f;
    vector< complex<double> > alphaOne;
    vector< complex<double> > alphaThree;
    
    wavefunction wfG; wavefunction wf1; wavefunction wf3;
    char wfGround[50]; char wfFirst[50]; char wfThird[50];
    sprintf(wfGround, "wf1D_R2.640000_0"); sprintf(wfFirst, "wf1D_R2.640000_1"); sprintf(wfThird, "wf1D_R2.640000_3");
    
    pert::Initialize(wfG, nPoint(), spatialStep(), 0);
    pert::Initialize(wf1, nPoint(), spatialStep(), 0);
    pert::Initialize(wf3, nPoint(), spatialStep(), 0);
    wfG.load(wfGround);
    wf1.load(wfFirst);
    wf3.load(wfThird);
    
    complex<double> dip01 = pObject.dipole(wfG, wf1);
    complex<double> dip03 = pObject.dipole(wfG, wf3);
    cout<<dip01<<"\t"<<dip03<<endl;
    
    pObject.setEnergies(omega, 0.6, 1.4, 0.01);
    
    for (int i = 0; i < omega.size(); i++)
    {
        complex<double> elmtOne = pObject.dipolePlaneWave(wf1, omega[i]);
        complex<double> elmtThree = pObject.dipolePlaneWave(wf3, omega[i]);
        dip1f.push_back(elmtOne);
        dip3f.push_back(elmtThree);
        alphaOne.push_back(elmtOne * dip01);
        alphaThree.push_back(elmtThree * dip03);
    }
    pObject.takeEnergy(E_m);
    pObject.createCosineSquare(fieldX, 0.01, E0, T, 0.8313, -0.5*pi);
    
    
    /*Tests*/
    cout<<pObject.secondIntegral(1, 1, T, 0.8313, fieldX, E_m, 0.01)<<endl;
    cout << "Go Dawgs!\n";
    return 0;
}
