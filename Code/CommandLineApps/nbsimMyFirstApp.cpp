/*=============================================================================

  PHAS0100ASSIGNMENT2: PHAS0100 Assignment 2 Gravitational N-body Simulation

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include <nbsimMyFunctions.h>
#include <nbsimExceptionMacro.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <list>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <thread>
#include <iostream>
#include <random>

// Example, header-only library, included in project for simplicity's sake.
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
// reading the SolarSystemData
#include "nbsimSolarSystemData.ipp"
// 
#include "nbsimBasicTypes.h"

/**
 * \brief Demo file to check that includes and library linkage is correct.
 */

namespace nbsim
{
  
  std::vector<std::shared_ptr<nbsim::MassiveParticle>> simulate_body()
  {
    // use random number
    std::default_random_engine e;
    std::uniform_real_distribution<double> rand_num(1e-10, 50.0);
    std::uniform_real_distribution<double> rand_angle(-3.1416, 3.1416);
    
    int tot_num=2000;
    double centre_Mu=rand_num(e);
    Eigen::Vector3d centre_position(0.0,0.0,0.0),centre_velocity(0.0,0.0,0.0);
    // initialize the first body
    std::vector<std::shared_ptr<nbsim::MassiveParticle>> sim_particles;
    std::shared_ptr<nbsim::MassiveParticle> first_ptr(new nbsim::MassiveParticle("-",centre_position, centre_velocity,centre_Mu));
    sim_particles.push_back(first_ptr);
    // initialize other bodies
    Eigen::Vector3d ini_position, ini_velocity;
    double angle=rand_angle(e);
    double r;
    double mu;
    for (int i=1; i<tot_num; i++) 
    {
      r=rand_num(e);
      ini_position<< r*sin(angle), r*cos(angle), 0.0;
      ini_velocity<< -sin(angle)/sqrt(r), cos(angle)/sqrt(r), 0.0;
      mu=rand_num(e);
      std::shared_ptr<nbsim::MassiveParticle> body_ptr_i(new nbsim::MassiveParticle("-",ini_position, ini_velocity,mu));
      sim_particles.push_back(body_ptr_i);
    }
    return sim_particles;
  };
  double kinetic_energy(std::vector<std::shared_ptr<nbsim::MassiveParticle>> bodies)
  {
    double E_k=0.0;
    #pragma omp parallel for
    for (int i=0; i< bodies.size();i++)
    {
      #pragma omp critical
      {
        E_k += (bodies[i]->getMu()*(bodies[i]->getVelocity().squaredNorm()))/2;
      }
    }
    return E_k;
  };

  // calculate potential energy
  double potential_energy(std::vector<std::shared_ptr<nbsim::MassiveParticle>> bodies)
  {
    double E_p=0.0;
    for (int i=0; i< bodies.size();i++)
    {
      for(int ii=0; ii<bodies[0]->attractors_ptr.size();ii++)
      {
        if (bodies[i]!=bodies[ii])
        {
          E_p += (-bodies[i]->getMu()*bodies[ii]->getMu()/(bodies[i]->getPosition()-bodies[ii]->getPosition()).norm())/2;
        }
      }
    }
    return E_p;
  };
}

int main(int argc, char** argv)
{
  // start of the clock
  std::clock_t clock_start = std::clock();
  auto chrono_start = std::chrono::high_resolution_clock::now();
  // initialize
  std::vector<std::shared_ptr<nbsim::MassiveParticle>> bodies;
  bodies = nbsim::simulate_body();
  #ifdef DEBUG_ON
  for (int i=0;i<bodies.size();i++)
  {
    std::cout<<bodies[i]->getPosition().transpose()<<" // "
    <<bodies[i]->getVelocity().transpose()<<"\n"<<std::endl;
  }
  #endif
  // add the attractors
  for (int i=0; i<bodies.size(); i++)
  {
		for (int ii=0; ii<bodies.size(); ii++)
    {
      if (i!=ii)
      {
        bodies[i]->addAttractor(bodies[ii]);
        // bodies[i]->addAttractor(bodies[ii]);
      }
		}
	}  
  // calculate the initial energy
  double K_energy_init=nbsim::kinetic_energy(bodies);
  double P_energy_init=nbsim::potential_energy(bodies);
  double Total_energy_init= K_energy_init+P_energy_init;
  // calculate the acceleration and update the position and velocity
  double step_size= 0.000274;
  omp_set_num_threads(8);
  #pragma omp parallel
  for (int i=0; i<10; i++)
  {
    #pragma omp for 
		for (int ii=0;ii<2000;ii++)
    {
			bodies[ii]->calculateAcceleration();
		}
    #pragma omp for nowait
		for (int iii=0;iii<2000;iii++)
    {
			bodies[iii]->integrateTimestep(step_size);
		}	
	}
  // cout the end of body positions
  #pragma omp parallel for
  for (int i=0;i<bodies.size();i++)
  {
		std::cout<<bodies[i]->name<<"\n end of position:"<<bodies[i]->getPosition().transpose()<<std::endl;
  }
  std::cout<<std::endl;

  // calcualte the end of total energy
  double K_energy_end=nbsim::kinetic_energy(bodies);
  double P_energy_end=nbsim::potential_energy(bodies);
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

  return 0;
}