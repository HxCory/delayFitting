#include <fftw3.h>
#include <time.h>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "constant.h"
#include "Cartesian1D.h"
#include "ParameterMap.h"
#include "InputStuff.h"
#include "Projection.h"
#include "Laser.h"
#include "Observable.h"
#include "utils.h"

using namespace std;
using namespace Cartesian_1D;  

double PI=3.14159265358979;

Parameter<int> n_point("n_point", 1000, "number of points for R");
Parameter<double> spatial_step("spatial_step", 0.1, "spatial_steps for R");

Parameter<double> sigma("sigma", 1.0, "width of the pulse");
Parameter<double> gauss_center("gauss_center", 0.0, "Center of Gaussian for ground state");
Parameter<double> initial_momentum("initial_momentum", 0.0, "Initial momentum of wavepacket");
Parameter<string> jobname("jobname", "/p5data2/cogo4490/example/H_eigenstate", "The Name of the job");
Parameter<string> outputFolder("outputFolder", "/p5data2/cogo4490/example", "The basename of the output directory");
Parameter<string> wavefunctionFolder("wavefunctionFolder", "/p5data2/cogo4490/example/Wavefunction/", "The basename of the output directory");
Parameter<double> dt("dt", 0.05, "Time step");
Parameter<double> Z("Z", 1.0, "Z in Coulomb potential");
Parameter<double> softcore("softcore", 0.001, "soft core parameter");
Parameter<int> free_propagation_steps("free_propagation_steps", 10000, "free propagation steps");
Parameter<int> x_symmetry("x_symmetry", 1, "symmetry in x direction");
Parameter<int> Nstate("Nstate", 3, "number of states we want to calculate");
Parameter<double> R_in("R_in", 2.0, "BO-fixed nuclear distance");

void openFile( ofstream &filename, string name )
{
  filename.open(name.c_str());
  filename.setf(ios::showpoint|ios::scientific);
}

void closeFile (ofstream &filename)
{
  filename.close();
}

int main( int argc, char *argv[] )
{

  cout<<"****************************************************************************************"<<endl;
  cout.precision(15); 

  time_t second_1 = time(NULL);
  parseInput(argc, argv);

  wavefunction psi;
  Initialize(psi, n_point(), spatial_step(), x_symmetry());
  Guess_Function_Gauss(psi, sigma(), gauss_center(), initial_momentum());
  Normalize(psi);

  Hamiltonian H(psi);

  wavefunction wf[Nstate()];
  for(int i=0; i<Nstate(); i++)
    {
      Initialize(wf[i], n_point(), spatial_step(), x_symmetry());
      Guess_Function_Gauss(wf[i], sigma(), gauss_center(), initial_momentum());
      Normalize(wf[i]);
    }

  double R=R_in();

  ofstream output_wf;
  openFile(output_wf, outputFolder()+"/wf.dat");
  output_wf.precision(15);

  ofstream output_c;
  openFile(output_c, outputFolder()+"/cout.dat");
  output_c.precision(15);

  Potential_H2_plus_e1D_revised pot(psi, R); //including 1/R in Ham
  Observable<double, wavefunction, ABVparam> energy( "Energy", 1, Obs_Energy, psi, pot );
  Observable<double, wavefunction> distance( "distance", 1, Obs_Expectation_Value_X, psi);

  double E_new=1.0, E_old=0.0, R_old=10.0, R_new=0.0;
  char wfname[50];
  char wfInd[50];

  for(int i=0; i<Nstate(); i++)
  {
      output_c<<"start to get the "<<i<<"th state:"<<endl;
      output_c<<"R"<<"\t"<<R_in()<<endl;
      E_new=1.0, E_old=0.0, R_old=10.0, R_new=0.0;
      //change threshold back to lower value later
      for(int counter=0; fabs(E_new-E_old)>=1e-10||fabs(R_old-R_new)>=5e-10; counter++)
      {
          H.X(-I*dt(), psi, pot);
          Normalize(psi);
          for(int j=0; j<i; j++)
            ProjectOUT_Diff_Sizes(psi, wf[j]);
          Normalize(psi);

          if(counter%2000==0)
          {
              E_old=E_new;
              E_new=energy.measure();
              R_old=R_new;
              R_new=distance.measure();



              output_c<<E_new<<"\t"<<R_new<<endl;
          }
      }
      Normalize(psi);

      //sprintf(wfname, "wf%f", R);
      sprintf(wfname, "wf1D_R%f_%d", R, i);
      PlaceWaveFunction(wf[i], psi);
      psi.save(wfname);

      ofstream outputInd;      
      sprintf(wfInd, "wf1D%d.txt", i);
      openFile(outputInd, outputFolder() + wfInd);

      for(int j=1; j<=wf[i].n1; j++)
        output_wf<<wf[i].x1[j]<<"\t"<<real(wf[i].wave[wf[i].in2(1, j)])<<"\t"<<imag(wf[i].wave[wf[i].in2(1, j)])<<endl;
      for(int j=1; j<=wf[i].n1; j++)
        outputInd<<wf[i].x1[j]<<"\t"<<real(wf[i].wave[wf[i].in2(1, j)])<<"\t"<<imag(wf[i].wave[wf[i].in2(1, j)])<<endl;
      closeFile(outputInd);
  }

  for(int k=0; k<Nstate(); k++)
  {
    output_c<<"start to test the "<<k<<"th state:"<<endl;
    PlaceWaveFunction(psi, wf[k]);
    for(int counter=0; counter<free_propagation_steps(); counter++)
    {
      H.X(dt(), psi, pot);

      Mask_Function(psi, 0.01, 0.01, 1.0/6.0);

      if(counter%2000==0)
      {
       E_new=energy.measure();
R_new=distance.measure();
        output_c<<E_new<<"\t"<<R_new<<"\t"<<norm(Project(psi, wf[k]))<<endl;
      }
    }

    for(int i=1; i<=wf[k].n1; i++)
      output_wf<<wf[k].x1[i]<<"\t"<<real(wf[k].wave[wf[k].in2(1, i)])<<"\t"<<imag(wf[k].wave[wf[k].in2(1, i)])<<endl;
  }

  output_c<<"****************************************************************************************"<<endl;
  output_c<<difftime(time(NULL),second_1)/60<<" minutes elapsed in running this program."<< endl;
}

