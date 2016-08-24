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
Parameter<string> outputFolder("outputFolder", "/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/");
Parameter<int> nState("nState", 4, "number of states");
Parameter<int> nPoint("nPoint", 8200, "number grid points");
Parameter<double> spatialStep("spatialStep", 0.0732, "dx");

int main(int argc, const char* argv[]) {
/*Output Stuff*/
    ofstream outputField;
    pert::openFile(outputField, outputFolder() + "field.txt");
    outputField.precision(15);
    
    ofstream outputFieldTime;
    pert::openFile(outputFieldTime, outputFolder() + "fieldTime.txt");
    outputFieldTime.precision(15);
    
    ofstream outputAlphaOne;
    pert::openFile(outputAlphaOne, outputFolder() + "alphaOne.txt");
    outputAlphaOne.precision(15);
    
   	ofstream outputAlphaThree;
    pert::openFile(outputAlphaThree, outputFolder() + "alphaThree.txt");
    outputAlphaThree.precision(15);
    
    ofstream outputOmega;
    pert::openFile(outputOmega, outputFolder() + "omega.txt");
    outputOmega.precision(15);
    
    ofstream outputImagCf;
    pert::openFile(outputImagCf, outputFolder() + "recf.txt");
    outputImagCf.precision(15);

    ofstream outputRealCf;
    pert::openFile(outputRealCf, outputFolder() + "imcf.txt");
    outputRealCf.precision(15);

/*Declarations*/
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
    vector<double> timer;
    vector< complex<double> > cf;
    vector<double> dummyField;
    complex<double> fac;

    wavefunction wfG; wavefunction wf1; wavefunction wf3;
    char wfGround[50]; char wfFirst[50]; char wfThird[50];

    /*Initialize Things*/
    sprintf(wfGround, "wf1D_R2.640000_0"); 
    sprintf(wfFirst, "wf1D_R2.640000_1"); 
    sprintf(wfThird, "wf1D_R2.640000_3");
    pert::Initialize(wfG, nPoint(), spatialStep(), 0);
    pert::Initialize(wf1, nPoint(), spatialStep(), 0);
    pert::Initialize(wf3, nPoint(), spatialStep(), 0);
    wfG.load(wfGround);
    wf1.load(wfFirst);
    wf3.load(wfThird);

/*Ops*/    
    complex<double> dip01 = pObject.dipole(wf1, wfG);
    complex<double> dip03 = pObject.dipole(wf3, wfG);    
    pObject.setEnergies(omega, 0.58, 1.4, 0.001);
    pObject.takeEnergy(E_m);

    for (int i = 0; i < omega.size(); i++)
    {
        complex<double> elmtOne = dip01 * pObject.dipolePlaneWave(wf1, 
        	pObject.getMomentum(omega[i]));
        complex<double> elmtThree = dip03 * pObject.dipolePlaneWave(wf3,
         pObject.getMomentum(omega[i]));
        alphaOne.push_back(elmtOne);
        alphaThree.push_back(elmtThree);
    }
    vector< vector<double> > fieldVector(omega.size(), vector<double> (timer.size(), 0.0));

    pObject.createCosineSquare(fieldVector, timer, pObject.dt0, E0, T, omega, -0.5*pi);
    
    for (int i = 0; i < timer.size(); i++)
    {
        for (int j = 0; j < omega.size(); j++)
        {
            outputField<<fieldVector[j][i]<<"\t";
        }
        outputField<<endl;
        outputFieldTime<<timer[i]<<endl;
    }

    
    for(int j = 0; j < omega.size(); j++)
    {
    	outputOmega<<omega[j]<<endl;
    	outputAlphaOne<<real(alphaOne[j])<<endl;
    	outputAlphaThree<<real(alphaThree[j])<<endl;
        fac = (alphaOne[j] * pObject.firstIntegral(1, T, omega[j], fieldVector[j],
                E_m, pObject.dt0)) + (alphaThree[j] * pObject.firstIntegral
                (3, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        outputRealCf<<real(fac)<<endl;
        outputImagCf<<imag(fac)<<endl;
        cout<<fac<<endl;
    }
    
/*Tests*/
    cout<<real(fac)<<"\t"<<imag(fac)<<endl;
    cout << "Go Dawgs!\n";
    return 0;
}
