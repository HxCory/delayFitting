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
,  firstEnergy(-0.8502)
,  secondEnergy(-0.4905)
,  thirdEnergy(-0.3278)
,  stateEnergy(4)
,  dt0(0.01)
{
    stateEnergy[0] = groundEnergy;
    stateEnergy[1] = firstEnergy;
    stateEnergy[2] = secondEnergy;
    stateEnergy[3] = thirdEnergy;
}

pert::~pert()
{
}

void pert::openFile(std::ofstream &filename, std::string name)
{
	filename.open(name.c_str());
}

void pert::closeFile(std::ofstream &filename)
{
	filename.close();
}

/*-------------------------------------------------------------------------------------*/


std::complex<double> pert::secondIntegral(int m, double tPrime, double &t,
	double omega, std::vector<double> &field, std::vector<double> &energy, double dt)
{
	std::complex<double> s;
	int numPoints = field.size();
	std::vector<double> time(numPoints, 0.0);
	int Nhalf=int(0.5 * t / dt);
	for(int i=0; i<numPoints; i++)
	{
		time[i] = (-Nhalf+i) * dt;
		
		if(time[i] <= tPrime)
		{
			s += exp(std::complex<double>(0.0, (energy[m] - groundEnergy)
				 * time[i])) * std::complex<double>(field[i], 0.0);
		}
	}

	return s * dt;
}

std::complex<double> pert::firstIntegral(int m, double &t, double &omega,
 				std::vector<double> &field, std::vector<double> &energy, 
 				double &dt)
{
	std::complex<double> s;
	int numPoints = field.size();
	std::vector<double> time(numPoints, 0.0);
	int Nhalf = int(0.5 * t / dt);
	for (int i = 0; i < numPoints; ++i)
	{
		time[i] = (-Nhalf + i) * dt;
		s += exp(std::complex<double> (0.0, ((2 * omega) 
			+ groundEnergy - energy[m])
			* time[i])) * std::complex<double> (field[i], 0.0) 
			* secondIntegral(m, time[i], t, omega, field, energy, dt);
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

void pert::setEnergies(std::vector<double> &energies, double min,
						 double max, double interval)
{
    double point = min;
    while (point <= max) {
        energies.push_back(point);
        point += interval;
    }
}

/*-------------------------------------------------------------------------------------*/


std::complex<double> pert::dipole(wavefunction &wf, wavefunction &wf2)
{
    complex<double> mu (0.0, 0.0);
    for (int i = 1; i <= wf.n1; i++)
    {
        mu += -conj(wf.wave[wf.in2(1, i)]) 
        	* std::complex<double> (wf.x1[i], 0.0) * wf2.wave[wf2.in2(1, i)];
    }
    mu *= wf.dx1;
    return mu;
}

std::complex<double> pert::dipolePlaneWave(wavefunction &wf, double k)
{
    complex<double> mu (0.0, 0.0);
    for (int i = 1; i <= wf.n1; i++)
    {
        mu += -(wf.wave[wf.in2(1, i)]) * complex<double> (wf.x1[i], 0.0) 
        		* exp(complex<double>(0.0, -k * wf.x1[i]));
    }
    mu *= wf.dx1;
    return mu;
}

double pert::getMomentum(double &omega)
{
	return sqrt(2*(2 * omega + groundEnergy));
} 

/*-------------------------------------------------------------------------------------*/

void pert::createCosineSquare(vector<vector<double>> &dummy,
 vector<double> &dummyTime,	double dt, double amp, double duration,
  vector<double> &dummyOmega, double phi)
{
    int Nhalf = int(floor(0.5 * duration / dt));
    int N = 2 * Nhalf;
    vector< vector<double> > s(N, vector<double>(dummyOmega.size(), 0.0));
    
    for (int i = 0; i < N; i++)
    {
   		double time = (-Nhalf + i) * dt;
   		dummyTime.push_back(time);

    	for (int j = 0; j < dummyOmega.size(); ++j)
    	{
            double currentOmega = dummyOmega[j];
       	 	s[i][j] = amp * pow(cos(pi * time / duration), 2.0) 
        				* cos(currentOmega * time + phi);
    	}
        
    }
    dummy = s;
}

void pert::takeEnergy(vector<double> &dummy)
{
    dummy = stateEnergy;
}

