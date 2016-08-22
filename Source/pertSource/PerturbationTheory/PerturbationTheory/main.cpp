//
//  main.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#include "pert.hpp"

using namespace std;

int main(int argc, const char* argv[]) {
	vector<double> omega; 	
	vector<int> m;
	vector<double> E_m;
	double E_0;
	double alpha;	
	double T;
	string OutputFolder;
	string dataFolder;
    static std::ifstream omegaDip;
	static std::ifstream dipOne;
    static std::ifstream dipThree;
    static std::vector<double> AlphaOne;

	pert::openFile(dipOne, "/Users/cgoldsmith/Desktop/text_files_data/alpha1_3fs.txt");
   	pert::readInput(dipOne);
   	pert::closeFile(dipOne);

   	pert::openFile(omegaDip, "/Users/cgoldsmith/Desktop/text_files_data/omega4dip.txt");
   	pert::readInput(omegaDip);
   	pert::closeFile(omegaDip);

    cout << "Go Dawgs!\n";
    return 0;
}
