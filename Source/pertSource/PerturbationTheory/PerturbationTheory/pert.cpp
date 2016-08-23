//
//  pert.cpp
//  PerturbationTheory
//
//  Created by C. Goldsmith on 8/19/16.
//  Copyright © 2016 C. Goldsmith. All rights reserved.
//

#include "pert.hpp"
#include <vector>
#include <complex>

template <typename T>
T StringToNumber ( const std::string &Text )
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

pert::pert()
	:  groundEnergy(-1.1591)
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

std::complex<double> pert::secondIntegral(int m, double tPrime, double t, 
	double omega, std::vector<double> field, std::vector<double> energy, double dt)
{
	std::complex<double> s = (0.0, 0.0);
	int numPoints=field.size();
	std::vector<double> time(numPoints, 0.0);
	int Nhalf=int(0.5 * t / dt);
	for(int i=0; i<numPoints; i++)
	{
		time[i] = (-Nhalf+i) * dt;
		
		if(time[i]<=tPrime)
		{
			s += exp(std::complex<double>(0.0, (energy[m] - groundEnergy) * time[i])) *
				std::complex<double>(field[i], 0.0);
		}
	}

	return s * dt;
}

std::complex<double> pert::firstIntegral(int m, double t, double omega,
 				std::vector<double> field, std::vector<double> energy, double dt)
{
	std::complex<double> s = (0.0, 0.0);
	int numPoints = field.size();
	std::vector<double> time(numPoints, 0.0);
	int Nhalf = int(0.5 * t / dt);
	for (int i = 0; i < numPoints; ++i)
	{
		time[i] = (-Nhalf + i) * dt;
		s += exp(std::complex<double> (0.0, ((2 * omega) + groundEnergy - energy[m])
			* time[i])) * std::complex<double> (field[i], 0.0) * secondIntegral(
				m, time[i], t, omega, field, energy, dt);
	}
	return s * dt;
}
