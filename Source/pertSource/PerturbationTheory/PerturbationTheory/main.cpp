//
//  main.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

// #include <iostream>
// #include <vector>
// #include <fstream>
#include "pert.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    pert::openFile(pert::omegaDip, "/Users/cgoldsmith/Desktop/text_files_data/alpha1_3fs.txt");
    pert::readDipoleInput(pert::omegaDip);
    cout << "Go Dawgs!\n";
    return 0;
}
