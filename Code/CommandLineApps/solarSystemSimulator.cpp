#include <nbsimMyFunctions.h>
#include <nbsimExceptionMacro.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <vector>
#include <iomanip>
#include <chrono>
#include <ctime>
#include<cstdlib>
#include <thread>
// Example, header-only library, included in project for simplicity's sake.
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
// reading the SolarSystemData
#include "nbsimSolarSystemData.ipp"
// 
#include "nbsimBasicTypes.h"

namespace nbsim
{
  // calculate kinetic energy
  double kinetic_energy(std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9])
  {
    double E_k=0.0;
    // #pragma omp parallel for reduction(+: E_k)
    for (int i=0; i< bodys_ptr->use_count();i++)
    {
      {
        E_k += (bodys_ptr[i]->getMu()*(bodys_ptr[i]->getVelocity().squaredNorm()))/2;
      }
    }
    return E_k;
  };

  // calculate potential energy
  double potential_energy(std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9])
  {
    double E_p=0.0;
    for (int i=0; i< bodys_ptr->use_count();i++)
    {
      // #pragma omp parallel for reduction(+: E_p)
      for(int ii=0; ii<bodys_ptr[0]->attractors_ptr.size();ii++)
      {
        if (bodys_ptr[i]!=bodys_ptr[ii])
        {
          E_p += (-bodys_ptr[i]->getMu()*bodys_ptr[ii]->getMu()/(bodys_ptr[i]->getPosition()-bodys_ptr[ii]->getPosition()).norm())/2;
        }
      }
    }
    return E_p;
  };

  void helpmessage()
  {
    std::cout<<" -- Help message"<<"\n"
             <<"    This is the solar system simulator.\n"
             <<"    You can choose the timestep(unit year), and numbers of timestep(unit year) to see the beginning and the end of the bodies in the solar system.\n"
             <<"    You can input two arguments in the terminal: ./bin/solarSystemSimulator (timestep) (number of timestep, must be double)\n"
             <<"    Then the program print the the beginning and the end of the bodies, their energy, running time and the maximum number of processor can be used.\n"
             <<"    For example: ./bin/solarSystemSimulator 0.000274 1.0\n\n"
             <<std::endl;
  };
  
  void errormessage()
  {
    std::cout<<" -- Error message"<<"\n"
             <<"    You can input:\n"
             <<"    ./bin/solarSystemSimulator -help\n"
             <<"    or\n"
             <<"    ./bin/solarSystemSimulator -h\n"
             <<"    to get help message."<<std::endl;
  };
}

// main function
int main(int argc, char** argv)
{
  // terminal input 0 arg, like: ./bin/solarSystemSimulator
  if (argc==1)
  {
    nbsim::helpmessage();
  }
  // terminal input 1 arg
  else if (argc==2)
  {
    std::string help =argv[1];
    if (help =="-help"||help =="-h")
    {
      nbsim::helpmessage();
    }
    else
    {
      std::cout<<"Error: Input arguments are not enough.\n\n";
      nbsim::errormessage();
    }
  }
  // terminal input 4 arg, like: ./bin/solarSystemSimulator 0.000274 1 16 1
  else if (argc>=4)
  {
    std::cout<<"Error: Input arguments are more then needed.\n\n";
    nbsim::errormessage();
  }
  // terminal input 2 or 3 arg
  else if (argc==3)
  {
  // start of the clock
  std::clock_t clock_start = std::clock();
  auto chrono_start = std::chrono::high_resolution_clock::now();

  double mu;
  std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9];
  // initialize the bodies
  for (int i=0; i<9; i++) 
  {
    std::shared_ptr<nbsim::MassiveParticle> body_ptr_i(new nbsim::MassiveParticle(nbsim::solarSystemData[i].name,nbsim::solarSystemData[i].position, nbsim::solarSystemData[i].velocity,nbsim::solarSystemData[i].mu));
    bodys_ptr[i]=body_ptr_i;
	}
  #ifdef DEBUG
  // cout the beginning of body positions
  std::cout<<"\nStart of position:\n\n";
  for (int i=0;i<9;i++)
  {
		std::cout<<bodys_ptr[i]->name<<": "<<bodys_ptr[i]->getPosition().transpose()<<std::endl;
  }
  std::cout<<std::endl;
  #endif
  // add the attractors
  for (int i=0; i<9; i++)
  {
		for (int ii=0; ii<9; ii++)
    {
      if (bodys_ptr[i]->name!=bodys_ptr[ii]->name)
      {
        bodys_ptr[i]->addAttractor(bodys_ptr[ii]);
      }
		}
	}  
  // calculate the initial energy
  double K_energy_init=nbsim::kinetic_energy(bodys_ptr);
  double P_energy_init=nbsim::potential_energy(bodys_ptr);
  double Total_energy_init= K_energy_init+P_energy_init;
  // calculate the acceleration and update the position and velocity
  double step_size=atof(argv[1]);
  int step_num= atoi(argv[2])*3650;
  // #pragma omp parallel
  for (int i=0; i<step_num; i++)
  {
    // #pragma omp for 
		for (int ii=0;ii<9;ii++)
    {
			bodys_ptr[ii]->calculateAcceleration();
		}
    // #pragma omp for nowait
		for (int iii=0;iii<9;iii++)
    {
			bodys_ptr[iii]->integrateTimestep(step_size);
		}	
	}
  #ifdef DEBUG
  // cout the end of body positions
  std::cout<<"\nEnd of position:\n\n";
  for (int i=0;i<9;i++)
  {
		std::cout<<bodys_ptr[i]->name<<": "<<bodys_ptr[i]->getPosition().transpose()<<std::endl;
  }
  std::cout<<std::endl;
  #endif
  // calcualte the end of total energy
  double K_energy_end=nbsim::kinetic_energy(bodys_ptr);
  double P_energy_end=nbsim::potential_energy(bodys_ptr);
  double Total_energy_end= K_energy_end+P_energy_end;
  std::cout<< "System initial energy are :"
  << " kinetic energy is "<<K_energy_init
  << ", potential energy is "<< P_energy_init
  << ", total_energy is "<<Total_energy_init<<"\n"<< std::endl;
  std::cout<< "System final energy are :"
  << " kinetic energy is "<< K_energy_end
  << ", potential energy is "<< P_energy_end
  << ", total_energy is "<< Total_energy_end<<"\n"<< std::endl;
  // end of the clock
  std::clock_t clock_end = std::clock();
  auto chrono_end = std::chrono::high_resolution_clock::now();
  std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
            << 1000.0 * (clock_end - clock_start) / CLOCKS_PER_SEC << " ms\n"
            << "Wall clock time has passed: "
            << std::chrono::duration<double, std::milli>(chrono_end-chrono_start).count()
            << " ms\n"<<std::endl;
  // get the maximum of the processor number
  int nProcessors=omp_get_max_threads();
  std::cout<<"Processor number: "<<nProcessors<<std::endl;
  }
  return 0;
}