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

double kinetic_energy(std::vector <nbsim::MassiveParticle> body_vect)
{
  double E_k=0.0;
  for (int i=0; i< body_vect.size();i++)
  {
    E_k += body_vect[i].Mu*(body_vect[i].getVelocity().squaredNorm());
  }
  return E_k/2;
};

double potential_energy(std::vector <nbsim::MassiveParticle> body_vec)
{
  double E_p=0.0;
  for (int i=0; i< body_vec.size();i++)
  {
    for(int ii=0; ii<body_vec[0].mass_particle_vec.size();ii++)
    {
      E_p += -body_vec[i].getMu()*body_vec[i].mass_particle_vec[ii].getMu()/(body_vec[i].getPosition()-body_vec[i].mass_particle_vec[ii].getPosition()).norm();
    }
  }
  return E_p/2;
};

int main(int argc, char** argv)
{
  std::string planet_name[9];
  Eigen::Vector3d init_position, init_velocity;
  double mu;
  std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9];

  for (int i=0; i<9; i++) {
    planet_name[i]=nbsim::solarSystemData[i].name;
    init_position=nbsim::solarSystemData[i].position;
    init_velocity=nbsim::solarSystemData[i].velocity;
    mu=nbsim::solarSystemData[i].mu;
    std::shared_ptr<nbsim::MassiveParticle> body_ptr_i(new nbsim::MassiveParticle(planet_name[i],init_position, init_velocity,mu));
    bodys_ptr[i]=body_ptr_i;
	}
  // check the body's information
/*   for(int i=0; i<9; i++)
  {
    std::cout<<bodys_ptr[i]->name<<"\n"<<bodys_ptr[i]->getPosition()<<"\n"<<bodys_ptr[i]->getVelocity()<<"\n"<<bodys_ptr[i]->getMu()<<"\n"<<std::endl;
  } */
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
  // calculate the acceleration and update the position and velocity
  double step_size= 0.000274;
  for (int i=0; i<3650; i++)
  {
		for (int ii=0;ii<9;ii++)
    {
			bodys_ptr[ii]->calculateAcceleration();
		}
		for (int iii=0;iii<9;iii++)
    {
			bodys_ptr[iii]->integrateTimestep(step_size);
		}	
	}
  // cout the end of body positions
  for (int i=0;i<9;i++)
  {
		std::cout<<planet_name[i]<<"\n end of position:"<<bodys_ptr[i]->getPosition().transpose()<<std::endl;
  }
  // delete the memory
/*   for (int i=0;i<9;i++)
  {
    delete bodys_ptr[i];
  } */

  return 0;
}