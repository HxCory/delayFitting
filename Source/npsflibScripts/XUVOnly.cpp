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

Parameter<int> n_point_small("n_point_small", 1000, "number of points for R"); //smaller grid for ground state.
Parameter<int> n_point("n_point", 2000, "number of points for R"); //grid size
Parameter<double> spatial_step("spatial_step", 0.1, "spatial_steps for R");
Parameter<string> jobname("jobname", "/p5data2/cogo4490/example/XUV_propagation", "The Name of the job");
Parameter<string> outputFolder("outputFolder", "/p5data2/cogo4490/example/", "The basename of the output directory");
Parameter<string> wavefunctionFolder("wavefunctionFolder", "/p5data2/cogo4490/example/Wavefunction/", "The basename of the output directory");
Parameter<string> InitialState("InitialState", "/p5data2/cogo4490/example/Wavefunction/", "The basename of the output directory");
Parameter<double> R_in("R_in", 2.64, "internuclear distance");
Parameter<double> dt("dt", 0.05, "Time step");
Parameter<double> Z("Z", 1.0, "Z in Yukawa potential");
Parameter<double> softcore("softcore", 2, "soft core parameter");
Parameter<double> omega_au_XUV("omega_au_XUV", 2.0, "photon energy of the XUV pulse");
Parameter<double> intensity_Wcm2_XUV("intensity_Wcm2_XUV", 2.9e14, "intensity of the XUV pulse (W/cm2)");
Parameter<double> phi_XUV("phi_XUV", 0.0, "phase for the XUV pulse (in uinits of pi)");
Parameter<double> tau_XUV_fs("tau_XUV_fs", 3.0, "pulse duration of the XUV pulse in fs");
Parameter<int> XUV_position("XUV_position", 1, "center position of the XUV pulse with respect to the IR");
Parameter<double> scan_step("scan_step", 2.0, "center position of the XUV pulse with respect to the IR");
Parameter<int> free_propagation_steps("free_propagation_steps", 30000, "free propagation steps");
Parameter<int> savestep("savestep", 3000, "save how many step");
Parameter<int> FTsavestep ("FTsavestep", 2, "how often save FFT");
Parameter<int> tau ("tau", 20, "delay between XUV and streaking pulses in au");
Parameter<double> cut_radius("cut_radius", 45.0, "where to cut bound wavefx");
Parameter<int> Nstate("Nstate", 1, "when start to cut the wavefunction");

void Open_File( ofstream &filename, string name ) 
{
  filename.open(name.c_str());
  filename.setf(ios::showpoint|ios::scientific);
}

int main( int argc, char *argv[] )
{

  cout<<"****************************************************************************************"<<endl;
  cout.precision(15);

  time_t second_1 = time(NULL);
  parseInput(argc, argv);
 
  ofstream output_field; 
  Open_File(output_field, outputFolder()+"/field.dat");
  output_field.precision(15);
 
  ofstream output_vp;
  Open_File(output_vp, outputFolder()+"vp.dat");
  output_vp.precision(15);
  
  ofstream output_c;
  Open_File(output_c, outputFolder()+"cout.dat");
  output_c.precision(15);

  ofstream output_p;
  Open_File(output_p, outputFolder()+"inputpy.dat");
  output_p.precision(15);

  ofstream output_field_tot;
  Open_File(output_field_tot, outputFolder()+"/field_t.dat");
  output_field_tot.precision(15);

  int Nhalf=int(floor(0.5 * tau_XUV/dt()));
  int N=Nhalf*2+1;
  double E_XUV=sqrt(intensity_Wcm2_XUV()/3.5101e+16);
  double cycle_num_XUV=(tau_XUV_fs()*omega_au_XUV())/(2*PI*0.02419);
  double XUV_center=XUV_position()*scan_step();
  double tau_XUV=cycle_num_XUV*(2*PI)/omega_au_XUV();
  vector<double> XUV(N, 0.0);

  for(int i=0; i<N; i++)
  {
    double time=(-Nhalf+i)*dt();
    if((time-XUV_center+0.5*tau_XUV)>=0&&(time-XUV_center+0.5*tau_XUV)<=tau_XUV)
      XUV[i]=E_XUV*pow(cos(PI*(time-XUV_center)/tau_XUV),2.0)*cos(omega_au_XUV*(time-XUV_center)+phi_XUV()*PI);
    output_field_tot<<time<<"\t"<<IR[i]<<"\t"<<XUV[i]<<"\t"<<IR[i]+XUV[i]<<endl;
  }
  output_c<<R_in()<<"\t"<<XUV_position()<<"\t"<<tau_XUV<<"\t"<<omega_au_XUV<<endl;


  vector<double> myField_total(N+free_propagation_steps(), 0.0);
  for(int i=0; i<N; i++)
    myField_total[i]=IR[i]+XUV[i];
  for(int i=N; i<N+free_propagation_steps(); i++)
    myField_total[i]=0.0;

  for(int i=0; i<N; i++){
    output_field<<i*dt()*0.02419<<"\t"<<myField_total[i]<<endl;
    output_vp<<i*dt()*0.02419<<"\t"<<A_IR[i]<<endl;
  }
  int Nomega=N;
  complex<double> sumfft=0.0;
  vector<double> omega(Nomega, 0.0);
  vector< complex<double> > Pomega(Nomega,0.0);

  for(int i=0; i<N; i++){
    omega[i]=(i+1)*(2*PI/(Nomega*dt()));
  }

  for(int j=0; j<Nomega; j++){
    sumfft=0.0;
    for(int i=0; i<N; i++){
      double time=(i*dt());
      sumfft+=XUV[i]*exp(I*time*omega[j])*dt();
    }
    Pomega[j] = sumfft;
    output_field_tot<<omega[j]<<"\t"<<real(Pomega[j])<<"\t"<<imag(Pomega[j])<<"\t"<<norm(Pomega[j])<<endl;
  }

  double R=R_in();
  double XUVpos=XUV_position();
  double omega_au=omega_au_XUV();
  wavefunction wf_small;
  Initialize(wf_small, n_point_small(), spatial_step(), 0);


  wavefunction wf[Nstate()]; 
  char wfname[50];
  char wfprop[50];
  output_c<<"still going?"<<endl;
  for(int i=0; i<Nstate(); i++)
  {
    Initialize(wf[i], n_point(), spatial_step(), 0);
    sprintf(wfname, "testEwf_R%f_%d", R, i);
    sprintf(wfprop, "wfprop_%2.4FauX%f%d", omega_au, XUVpos, i);
    output_c<<"loadwfx"<<endl;
    wf_small.load(wfname);  
  }
  wavefunction psi;
  Initialize(psi, n_point(), spatial_step(), 0);
  output_c<<"place wfx"<<endl;
  PlaceWaveFunction(psi, wf_small);
  Hamiltonian H(psi);

  double t_start=0.5*tau_IR+XUV_center-0.5*tau_XUV;
  int index_start=((int)floor((t_start/dt())/200))*200;
  vector<double> potential(n_point(), 0.0);
  
  wavefunction psiFT;
      Initialize(psiFT, n_point(), spatial_step(), 0);
  
  ofstream output_wf;   
  Open_File(output_wf, outputFolder()+"/wf.dat");
  output_wf.precision(15);

  Potential_H2_plus_e1D_revised pot(psi, R);  

  double prob_loss=0.0, t_new=0.0;
  for(int time_index=index_start+1; time_index<myField_total.size(); time_index++)
  {
    H.X_ECS_E_R(dt(), psi, pot, myField_total[time_index], 0.02);
    
   
    if(time_index%savestep()==0)
    {     
      for(int j=1; j<=psi.n1; j++)
      {
  output_wf<<"\t"<<psi.x1[j]<<"\t"<<norm(psi.wave[psi.in2(1, j)])<<"\n";    
      }
       output_wf<<endl;
      
      t_new=(time_index+1)*dt();
      prob_loss=1.0-pow(Obs_Norm(psi), 2.0);
      output_c<<t_new*0.02419<<"\t";
      for(int i=0; i<Nstate(); i++)
      {
  double prob=norm(Project(psi, wf[i]));
  output_c<<prob<<"\t";     
      }
      output_c<<prob_loss<<"\t"<<endl;
    }
  }
  psi.save(wfprop);

  output_c<<"****************************************************************************************"<<endl;
  output_c<<difftime(time(NULL),second_1)/60<<" minutes elapsed in running this program."<< endl;
  return 0;
}
