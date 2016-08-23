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
Parameter<int> nState("nState", 4, "number of states");
Parameter<int> nPoint("nPoint", 8200, "number grid points");
Parameter<double> spatialStep("spatialStep", 0.0732, "dx");


int main(int argc, const char* argv[]) {
	vector<double> omega; 	
	vector<int> m;
	vector<double> E_m;
    double E0 = sqrt(1e+14/3.5101e+16);
	double alpha;	
	double T;
	string OutputFolder;
	string dataFolder;
    
    std::cout<<E0<<std::endl;
    
    wavefunction wfG; wavefunction wf1; wavefunction wf3;
    char wfGround[50]; char wfFirst[50]; char wfThird[50];
    sprintf(wfGround, "wf1D_R2.640000_0"); sprintf(wfFirst, "wf1D_R2.640000_1"); sprintf(wfThird, "wf1D_R2.640000_3");
    Initialize(wfG, nPoint(), spatialStep(), 0);
    wfG.load(wfGround);
    wf1.load(wfFirst);
    wf3.load(wfThird);
    
    //    static std::ifstream omegaDip;
//	static std::ifstream dipOne;
//    static std::ifstream dipThree;
//    static std::vector<double> AlphaOne;

//	pert::openFile(dipOne, "/Users/cgoldsmith/Desktop/text_files_data/alpha1_3fs.txt");
//   	pert::readInput(dipOne);
//   	pert::closeFile(dipOne);
//
//   	pert::openFile(omegaDip, "/Users/cgoldsmith/Desktop/text_files_data/omega4dip.txt");
//   	pert::readInput(omegaDip);
//   	pert::closeFile(omegaDip);

    cout << "Go Dawgs!\n";
    return 0;
}
