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
int main(int argc, char** argv)
{
  // 读取 data
  // int ibody = 0;
  // std::cout << nbsim::solarSystemData.at(ibody).name << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).mu << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).position << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).velocity << std::endl;

  Eigen::Vector3d a(0,0,0);
  nbsim::particle k1;
  // a << nbsim::solarSystemData.at(0).position;
  // std::cout << a << std::endl;
  k1.position << a;
  std::cout << k1.getPosition() << std::endl;



  return 0;
}

int other()
{
/*   double mu;
                    std::vector<std::shared_ptr<nbsim::MassiveParticle>> bodieslist;
                    Eigen::Vector3d pos, vel;
                    for (auto body : nbsim::solarSystemData)
                    {
                        mu = body.mu;
                        pos = body.position;
                        vel = body.velocity;

                        std::shared_ptr<nbsim::MassiveParticle> ptr_body(new nbsim::MassiveParticle(pos, vel, mu));
                        bodieslist.push_back(ptr_body);
                    } */
  return 0;
}