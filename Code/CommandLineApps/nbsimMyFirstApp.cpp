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

}

int main(int argc, char** argv)
{
  std::vector<std::shared_ptr<nbsim::MassiveParticle>> bodies;
  bodies = nbsim::simulate_body();
  for (int i=0;i<10;i++)
  {
    std::cout<<bodies[i]->getPosition().transpose()<<" // "
    <<bodies[i]->getVelocity().transpose()<<"\n"<<std::endl;
  }
  return 0;
}