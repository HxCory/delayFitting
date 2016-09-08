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
#include <Python/Python.h>

using namespace std;

Parameter<string> wavefunctionFolder("wavefunctionFolder",
 "/Users/cgoldsmith/repos/delayFitting/Data/eigenstates/WF", "directory location");
Parameter<string> outputFolder("outputFolder", "/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/MoreStateTests/");
//Parameter<string> outputFolder("outputFolder", "/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/");
Parameter<int> nState("nState", 4, "number of states");
Parameter<int> nPoint("nPoint", 8200, "number grid points");
Parameter<double> spatialStep("spatialStep", 0.0732, "dx");
Parameter<double> dipTwoTest("dipTwoTest", 0.0004086301, "test value for dip02");
//Parameter<double> alphaTwoTest("alphaTwoTest", -0.00009408645, "test value for alphaTwo");


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
    vector< complex<double> > alphaFour;
    vector< complex<double> > alphaFive;
    vector<double> timer;
    vector< complex<double> > cf;
    vector<double> dummyField;
    complex<double> fac;
    complex<double> facOne;
    complex<double> facTwo;
    complex<double> facThree;
    complex<double> facFour;
    complex<double> facFive;
    
    wavefunction wfG, wf1, wf2, wf3, wf4, wf5;
    char wfGround[50], wfFirst[50], wfSecond[50], wfThird[50], wfFourth[50], wfFifth[50];

    /*Initialize Things*/
    sprintf(wfGround, "testEwf_R2.640000_0");
    sprintf(wfFirst, "testEwf_R2.640000_1");
    sprintf(wfSecond, "testEwf_R2.640000_2");
    sprintf(wfThird, "testEwf_R2.640000_3");
    sprintf(wfFourth, "testEwf_R2.640000_4");
    sprintf(wfFifth, "testEwf_R2.640000_5");
    pert::Initialize(wfG, nPoint(), spatialStep(), 0);
    pert::Initialize(wf1, nPoint(), spatialStep(), 0);
    pert::Initialize(wf2, nPoint(), spatialStep(), 0);
    pert::Initialize(wf3, nPoint(), spatialStep(), 0);
    pert::Initialize(wf4, nPoint(), spatialStep(), 0);
    pert::Initialize(wf5, nPoint(), spatialStep(), 0);

    wfG.load(wfGround);
    wf1.load(wfFirst);
    wf2.load(wfSecond);
    wf3.load(wfThird);
    wf4.load(wfFourth);
    wf5.load(wfFifth);

/*Ops*/    
    complex<double> dip01 = pObject.dipole(wf1, wfG);
    // complex<double> dip02 = pObject.dipole(wf2, wfG);
    auto dip02 = complex<double> (0.0, dipTwoTest());
    complex<double> dip03 = pObject.dipole(wf3, wfG);    
    complex<double> dip04 = pObject.dipole(wf4, wfG);
    complex<double> dip05 = pObject.dipole(wf5, wfG);

    pObject.setEnergies(omega);
    pObject.takeEnergy(E_m);
    cout<<dip01<<"\t"<<dip02<<"\t"<<dip03<<"\t"<<dip04<<"\t"<<dip05<<endl;

    for (int i = 0; i < omega.size(); i++)
    {
        complex<double> elmtOne = dip01 * pObject.dipolePlaneWave(wf1,
        	pObject.getMomentum(omega[i]));
        complex<double> elmtTwo = dip02 * pObject.dipolePlaneWave(wf2,
            pObject.getMomentum(omega[i]));
        complex<double> elmtThree = dip03 * pObject.dipolePlaneWave(wf3,
         	pObject.getMomentum(omega[i]));
        complex<double> elmtFour = dip04 * pObject.dipolePlaneWave(wf4,
            pObject.getMomentum(omega[i]));
        complex<double> elmtFive = dip05 * pObject.dipolePlaneWave(wf5,
            pObject.getMomentum(omega[i]));
        alphaOne.push_back(elmtOne);
        alphaTwo.push_back(elmtTwo);
        alphaThree.push_back(elmtThree);
        alphaFour.push_back(elmtFour);
        alphaFive.push_back(elmtFive);
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
        
        facOne = (alphaOne[j] * pObject.firstIntegral
                  (1, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facTwo = (alphaTwo[j] * pObject.firstIntegral
                  (2, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facThree = (alphaThree[j] * pObject.firstIntegral
                  (3, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facFour = (alphaFour[j] * pObject.firstIntegral
                    (4, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        facFive = (alphaFive[j] * pObject.firstIntegral
                    (5, T, omega[j], fieldVector[j], E_m, pObject.dt0));
        
        fac = facOne + facTwo + facThree;// + facFour + facFive;
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
    cout<<real(facFour)<<"\t"<<imag(facFour)<<endl;
//    cout<<real(facFive)<<"\t"<<imag(facFive)<<endl;
    cout << "Go Dawgs!\n";
    return 0;
}
