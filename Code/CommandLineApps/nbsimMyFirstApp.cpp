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
  // 读取 data
  // int ibody = 0;
  // std::cout << nbsim::solarSystemData.at(ibody).name << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).mu << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).position << std::endl;
  // std::cout << nbsim::solarSystemData.at(ibody).velocity << std::endl;

  // Eigen::Vector3d pos(1,0,0);
  // Eigen::Vector3d vel(1,0,0);
  // const Eigen::Vector3d acc(2,0,0);

  // nbsim::Particle k1;
  // // a << nbsim::solarSystemData.at(0).position;
  // // std::cout << a << std::endl;
  // k1.position << pos;
  // k1.velocity << vel;
  // double k=1.0;
  // // std::cout << k1.getPosition() << std::endl;
  // k1.integrateTimestep(acc,k);
  // std::cout << k1.getPosition() << std::endl;
  // std::cout << k1.getVelocity() << std::endl;

  // std::string kk=typeid(k1);
/*   std::cout << typeid(k1).name() << std::endl;
  std::cout <<acc << std::endl; */

  // test massiveparticle
/*   nbsim::MassiveParticle x1;
  nbsim::MassiveParticle x2;

  x1.Mu=1;
  x2.Mu=1;

  x1.position<<1.,.0,.0;
  x2.position<<-1.,.0,.0;

  x1.velocity<<.0,0.5,.0;
  x2.velocity<<.0,-0.5,.0;

  x1.addAttractor(x2);
  x2.addAttractor(x1);

  x1.calculateAcceleration();
  x2.calculateAcceleration();
  std::cout <<"x1 :acc "<<x1.acc << std::endl;
  x1.integrateTimestep(1.0);
  x2.integrateTimestep(1.0);
  std::cout <<"x1 :position "<<x1.position << std::endl;
  std::cout <<"x1 :vel "<<x1.velocity << std::endl;

  std::cout<<sqrt((x1.position-x2.position).dot(x1.position-x2.position))<<std::endl; */

    // initialize the bodies
/*     nbsim::MassiveParticle Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune;
    std::vector <nbsim::MassiveParticle> body_vec;
    body_vec.push_back(Sun);
    body_vec.push_back(Mercury);
    body_vec.push_back(Venus);
    body_vec.push_back(Earth);
    body_vec.push_back(Mars);
    body_vec.push_back(Jupiter);
    body_vec.push_back(Saturn);
    body_vec.push_back(Uranus);
    body_vec.push_back(Neptune);
    for(int i=0; i<9;i++)
    {
        body_vec[i].name=nbsim::solarSystemData.at(i).name;
        body_vec[i].Mu=nbsim::solarSystemData.at(i).mu;
        body_vec[i].position=nbsim::solarSystemData.at(i).position;
        body_vec[i].velocity=nbsim::solarSystemData.at(i).velocity;
    };
 */
    nbsim::MassiveParticle Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune;
    std::vector <nbsim::MassiveParticle> body_vec;
    body_vec.push_back(Sun);
    body_vec.push_back(Mercury);
    body_vec.push_back(Venus);
    body_vec.push_back(Earth);
    body_vec.push_back(Mars);
    body_vec.push_back(Jupiter);
    body_vec.push_back(Saturn);
    body_vec.push_back(Uranus);
    body_vec.push_back(Neptune);
    for(int i=0; i<9;i++)
    {
        body_vec[i].name=nbsim::solarSystemData.at(i).name;
        body_vec[i].Mu=nbsim::solarSystemData.at(i).mu;
        body_vec[i].position=nbsim::solarSystemData.at(i).position;
        body_vec[i].velocity=nbsim::solarSystemData.at(i).velocity;
    };
    for(int iii=0; iii<9;iii++)
    {
        for (int ii=0; ii<9;ii++)
        {
            if(body_vec[iii].name!=body_vec[ii].name)
            {
                body_vec[iii].addAttractor(body_vec[ii]);
            }
        }
    };
    // 测试所有的body——vec的数量
    // for(int i=0; i<9;i++)
    // {
    //   std::cout<< body_vec[i].name<<": \n"<<std::endl;
    //   std::cout<< body_vec[i].Mu<<std::endl;  
    //   std::cout<< body_vec[i].position<<std::endl;  
    //   std::cout<< body_vec[i].velocity<<std::endl;  
    //   for (int ii=0; ii<body_vec[ii].mass_particle_vec.size();ii++)
    //   {
    //     std::cout<< body_vec[i].mass_particle_vec[ii].name<<std::endl;       

    //   } 
    //   std::cout<<"\n"<<std::endl;
    // };

    // calcualte initial total energy
    double K_energy_initial=kinetic_energy(body_vec);
    double P_energy_initial=potential_energy(body_vec);
    double Total_energy_initial= K_energy_initial+P_energy_initial;
    std::cout<< K_energy_initial<<"\n"<< std::endl;
    std::cout<< P_energy_initial<<"\n"<< std::endl;
    std::cout<< Total_energy_initial<<"\n"<<std::endl;

    double timestep= 0.000274;
    for (int i=0;i<3650;i++)
    {
        for (int ii=0;ii<body_vec.size();ii++)
        {
            body_vec[ii].calculateAcceleration();
            // std::cout<<body_vec[ii].name<<" acc are: \n";
            // std::cout<<body_vec[ii].acc<<std::endl;

        };

        for (int iii=0;iii<body_vec.size();iii++)
        {
            body_vec[iii].integrateTimestep(timestep);
            // std::cout<<"positions are: \n";
            // std::cout<<body_vec[iii].position<<std::endl;

        };
        // std::cout<<body_vec[0].position<<std::endl;
    };

    std::cout<<"Final positions are: \n";
    for (int i=0;i<body_vec.size();i++)
    {
        std::cout<<body_vec[i].name<<": \n";
        std::cout<<body_vec[i].position<<"\n\n";
    }

    // calcualte total energy
    double K_energy_end=kinetic_energy(body_vec);
    double P_energy_end=potential_energy(body_vec);
    double Total_energy_end= K_energy_end+P_energy_end;
    std::cout<< K_energy_end<<"\n"<< std::endl;
    std::cout<< P_energy_end<<"\n"<< std::endl;
    std::cout<< Total_energy_end<<"\n"<< std::endl;

  return 0;
}