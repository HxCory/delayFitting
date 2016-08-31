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
Parameter<string> outputFolder("outputFolder", "/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/Tests/");
Parameter<int> nState("nState", 4, "number of states");
Parameter<int> nPoint("nPoint", 8200, "number grid points");
Parameter<double> spatialStep("spatialStep", 0.0732, "dx");
Parameter<double> alphaTwoTest("alphaTwoTest", 0.00091, "test value for alphaTwo");

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
    
    ofstream outputAlphaTwo;
    pert::openFile(outputAlphaTwo, outputFolder() + "alphaTwo.txt");
    outputAlphaTwo.precision(15);

   	ofstream outputAlphaThree;
    pert::openFile(outputAlphaThree, outputFolder() + "alphaThree.txt");
    outputAlphaThree.precision(15);
    
    ofstream outputOmega;
    pert::openFile(outputOmega, outputFolder() + "omega.txt");
    outputOmega.precision(15);
    
    ofstream outputImagCfOne;
    pert::openFile(outputImagCfOne, outputFolder() + "imcfOne.txt");
    outputImagCfOne.precision(15);

    ofstream outputImagCfTwo;
    pert::openFile(outputImagCfTwo, outputFolder() + "imcfTwo.txt");
    outputImagCfTwo.precision(15);

    ofstream outputImagCfThree;
    pert::openFile(outputImagCfThree, outputFolder() + "imcfThree.txt");
    outputImagCfThree.precision(15);

    ofstream outputRealCfOne;
    pert::openFile(outputRealCfOne, outputFolder() + "recfOne.txt");
    outputRealCfOne.precision(15);

    ofstream outputRealCfTwo;
    pert::openFile(outputRealCfTwo, outputFolder() + "recfTwo.txt");
    outputRealCfTwo.precision(15);

    ofstream outputRealCfThree;
    pert::openFile(outputRealCfThree, outputFolder() + "recfThree.txt");
    outputRealCfThree.precision(15);

    ofstream outputRealCf;
    pert::openFile(outputRealCf, outputFolder() + "recf.txt");
    outputRealCf.precision(15);

    ofstream outputImagCf;
    pert::openFile(outputImagCf, outputFolder() + "imcf.txt");
    outputImagCf.precision(15);

/*Declarations*/
    pert pObject;
	vector<double> omega;
    vector<int> m;
	vector<double> E_m;
    double E0 = sqrt(1e+14/3.5101e+16);
	double T = 3/0.02419;
	
    vector<double> fieldX;
    vector< complex<double> > dip1f;
    vector< complex<double> > dip2f;
    vector< complex<double> > dip3f;
    vector< complex<double> > alphaOne;
    vector< complex<double> > alphaTwo;
    vector< complex<double> > alphaThree;
    vector<double> timer;
    vector< complex<double> > cf;
    vector<double> dummyField;
    complex<double> fac;
    complex<double> facOne;
    complex<double> facTwo;
    complex<double> facThree;
    
    wavefunction wfG, wf1, wf2, wf3;
    char wfGround[50], wfFirst[50], wfSecond[50], wfThird[50];

    /*Initialize Things*/
    sprintf(wfGround, "wf1D_R2.640000_0"); 
    sprintf(wfFirst, "wf1D_R2.640000_1"); 
    sprintf(wfSecond, "wf1D_R2.640000_2");
    sprintf(wfThird, "wf1D_R2.640000_3");
    pert::Initialize(wfG, nPoint(), spatialStep(), 0);
    pert::Initialize(wf1, nPoint(), spatialStep(), 0);
    pert::Initialize(wf2, nPoint(), spatialStep(), 0);
    pert::Initialize(wf3, nPoint(), spatialStep(), 0);
    wfG.load(wfGround);
    wf1.load(wfFirst);
    wf2.load(wfSecond);
    wf3.load(wfThird);

/*Ops*/    
    complex<double> dip01 = pObject.dipole(wf1, wfG);
    complex<double> dip02 = pObject.dipole(wf2, wfG);
    complex<double> dip03 = pObject.dipole(wf3, wfG);    
    pObject.setEnergies(omega);
    pObject.takeEnergy(E_m);
    cout<<dip01<<"\t"<<dip02<<"\t"<<dip03<<endl;

    for (int i = 0; i < omega.size(); i++)
    {
        complex<double> elmtOne = dip01 * pObject.dipolePlaneWave(wf1,
        	pObject.getMomentum(omega[i]));
        complex<double> elmtTwo = dip02 * pObject.dipolePlaneWave(wf2,
            pObject.getMomentum(omega[i]));
        complex<double> elmtThree = dip03 * pObject.dipolePlaneWave(wf3,
         	pObject.getMomentum(omega[i]));
        alphaOne.push_back(elmtOne);
        alphaTwo.push_back(alphaTwoTest());
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
        outputAlphaTwo<<real(alphaTwo[j])<<"\t"<<imag(alphaTwo[j])<<endl;
    	outputAlphaThree<<real(alphaThree[j])<<endl;
        fac = (alphaOne[j] * pObject.firstIntegral(1, T, omega[j], fieldVector[j],
               E_m, pObject.dt0)) + (alphaTwo[j] * pObject.firstIntegral
                (2, T, omega[j], fieldVector[j], E_m, pObject.dt0)) + (alphaThree[j] * pObject.firstIntegral
                (3, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facOne = (alphaOne[j] * pObject.firstIntegral
                  (1, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facTwo = (alphaTwo[j] * pObject.firstIntegral
                  (2, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facThree = (alphaThree[j] * pObject.firstIntegral
                  (3, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        
        outputRealCf<<real(fac)<<endl;
        outputImagCf<<imag(fac)<<endl;
        outputRealCfOne<<real(facOne)<<endl;
        outputImagCfOne<<imag(facOne)<<endl;
        outputRealCfTwo<<real(facTwo)<<endl;
        outputImagCfTwo<<imag(facTwo)<<endl;
        outputRealCfThree<<real(facThree)<<endl;
        outputImagCfThree<<imag(facThree)<<endl;
        cout<<fac<<"\t"<<facTwo<<endl;
    }
    
/*Tests*/
    cout<<real(facTwo)<<"\t"<<imag(facTwo)<<endl;
    cout << "Go Dawgs!\n";
    return 0;
}
