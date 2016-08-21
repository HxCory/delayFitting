//
//  pert.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#include "pert.hpp"
#include <vector>

pert::pert()
{
    // pert::reAlphaOne.push_back(0.0);
    // pert::imAlphaOne.push_back(0.0);
}

pert::~pert()
{
}

void pert::readDipoleInput(std::ifstream &input)
{
	pert Pert;
    std::string line;

	while(getline(input, line))
	{
		std::stringstream linestream(line);
		std::string keyword;

		if(linestream >> keyword)
		{
			double realDip;
			double imagDip;
		
			if(linestream>>realDip>>imagDip)
			{
				Pert.populateAlphaOne(realDip, imagDip);
			}
		}
	}
}

void populateAlphaOne(double real, double imag)
{
    pert::reAlphaOne.push_back(real);
    pert::imAlphaOne.push_back(imag);
    std::sort(pert::reAlphaOne.begin(), pert::reAlphaOne.end());
    std::sort(pert::imAlphaOne.begin(), pert::imAlphaOne.end());

}

void pert::openFile(std::ifstream &filename, std::string name)
{
	filename.open(name.c_str());
}

void pert::closeFile(std::ifstream &filename)
{
	filename.close();
}