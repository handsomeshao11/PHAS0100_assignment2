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

double kinetic_energy(std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9])
{
  double E_k=0.0;
  for (int i=0; i< bodys_ptr->use_count();i++)
  {
    E_k += (bodys_ptr[i]->getMu()*(bodys_ptr[i]->getVelocity().squaredNorm()))/2;
  }
  return E_k;
};

double potential_energy(std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9])
{
  double E_p=0.0;
  for (int i=0; i< bodys_ptr->use_count();i++)
  {
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

int main(int argc, char** argv)
{
  std::string bodys_name[9];
  Eigen::Vector3d init_position, init_velocity;
  double mu;
  std::shared_ptr<nbsim::MassiveParticle> bodys_ptr[9];
  // initialize the bodies
  for (int i=0; i<9; i++) {
    bodys_name[i]=nbsim::solarSystemData[i].name;
    init_position=nbsim::solarSystemData[i].position;
    init_velocity=nbsim::solarSystemData[i].velocity;
    mu=nbsim::solarSystemData[i].mu;
    std::shared_ptr<nbsim::MassiveParticle> body_ptr_i(new nbsim::MassiveParticle(bodys_name[i],init_position, init_velocity,mu));
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
		std::cout<<bodys_name[i]<<"\n end of position:"<<bodys_ptr[i]->getPosition().transpose()<<std::endl;
  }
  // delete the memory
/*   for (int i=0;i<9;i++)
  {
    delete bodys_ptr[i];
  } */
  // std::cout<<bodys_ptr->use_count()<<std::endl;

  // // calcualte total energy
  double K_energy_end=kinetic_energy(bodys_ptr);
  double P_energy_end=potential_energy(bodys_ptr);
  double Total_energy_end= K_energy_end+P_energy_end;
  std::cout<< K_energy_end<<"\n"<< std::endl;
  std::cout<< P_energy_end<<"\n"<< std::endl;
  std::cout<< Total_energy_end<<"\n"<< std::endl;

  return 0;
}