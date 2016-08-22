//
//  pert.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#include "pert.hpp"
#include <vector>

template <typename T>
T StringToNumber ( const std::string &Text )
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

pert::pert()
{
//	pert::a1.push_back(0.0);
//    std::cout<<a1[0]<<std::endl;
}

pert::~pert()
{
}

void pert::readInput(std::ifstream &input)
{
	pert Pert;
    std::string line;

	while(getline(input, line))
	{
		double realDip;
		std::string keyword;
		std::stringstream linestream(line);
		
		if(linestream >> keyword)
		{
            realDip = StringToNumber<double>(keyword);
            std::cout<<realDip<<std::endl;
//            a1.push_back(realDip);
        }
	}
}

void pert::openFile(std::ifstream &filename, std::string name)
{
	filename.open(name.c_str());
}

void pert::closeFile(std::ifstream &filename)
{
	filename.close();
}

