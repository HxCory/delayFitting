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

/*-------------------------------------------------------------------------------------*/


std::complex<double> pert::secondIntegral(int m, double tPrime, double t, 
	double omega, std::vector<double> field, std::vector<double> energy, double dt)
{
	std::complex<double> s;
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
	std::complex<double> s;
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

/*-------------------------------------------------------------------------------------*/

void pert::Initialize(wavefunction &wf, int nPoint, double spatialStep, int symm)
{
    wf.verbose = false;
    
    /*Initialize symmetry*/
  
    wf.symmetry_x1 = symm;
    wf.n1 = nPoint;

    /*Allocate the wavefx*/
    wf.wave.resize((wf.n1 + 2) * (2) ,0.);

    /*Allocate grid*/
    abs(symm) == 1
        ? wf.x1.type = HalfAxis
        : wf.x1.type = FullAxis;
    
    wf.x1.Init("x", nPoint, spatialStep, wf.x1.type);
    wf.dp1 = 2 * pi / (wf.dx1 * wf.n1);
    wf.p1.Init("p_x", nPoint, wf.dp1, FFTWAxis);
    
    /*Initialize spatial grid helper constants*/
    wf.one_by_dx1 = 1. / wf.dx1;
    wf.one_by_dx1sqr = 1. / (wf.dx1 * wf.dx1);
}

/*-------------------------------------------------------------------------------------*/


std::complex<double> pert::dipole(wavefunction &wf, wavefunction &wf2)
{
    complex<double> mu (0.0, 0.0);
    for (int i = 1; i <= wf.n1; i++)
    {
        mu += -conj(wf.wave[wf.in2(1, i)]) * std::complex<double> (wf.x1[i], 0.0) * wf2.wave[wf2.in2(1, i)];
    }
    mu *= wf.dx1;
    return mu;
}

std::complex<double> pert::dipolePlaneWave(wavefunction &wf, double &k)
{
    complex<double> mu (0.0, 0.0);
    for (int i = 1; i <= wf.n1; i++)
    {
        mu += -(wf.wave[wf.in2(1,i)]) * complex<double> (wf.x1[i], 0.0) * exp(complex<double>(0.0, -k * wf.x1[i]));
    }
    mu *= wf.dx1;
    return mu;
}

