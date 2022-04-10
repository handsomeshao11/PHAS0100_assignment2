/*=============================================================================

  PHAS0100ASSIGNMENT2: PHAS0100 Assignment 2 Gravitational N-body Simulation

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "nbsimBasicTypes.h"

namespace nbsim {
  void acc_not_zero(Eigen::Vector3d acceleration)
  {
    Eigen::Vector3d zero(0,0,0);
    if (zero.isApprox(acceleration))
    throw std::logic_error("Input acceleration is 0.");
  }

  void acc_is_const(Eigen::Vector3d acceleration)
  {
    if (typeid(acceleration)!=typeid(const Eigen::Vector3d))
    throw std::logic_error("Input acceleration type is not const.");// ！！！ debug

  }
  // class Particle
  Eigen::Vector3d Particle::getPosition()
  {
    return position;
  };

  Eigen::Vector3d Particle::getVelocity()
  {
    return velocity;
  };

  void Particle::integrateTimestep(Eigen::Vector3d acceleration, double timestep)
  {
    // check acceleration
    acc_not_zero(acceleration);
    acc_is_const(acceleration);
    // update position
    Eigen::Vector3d step=timestep*velocity;
    position+=step;

    // update velocity
    acceleration*=timestep;
    velocity+=acceleration;
  };

  // class MassiveParticle
  double MassiveParticle::getMu()
  {
    return Mu;
  };

  void MassiveParticle::addAttractor(MassiveParticle Mass_instance)
  {
    mass_particle_vec.push_back(Mass_instance);
  };

  void MassiveParticle::removeAttractor(MassiveParticle Mass_instance)
  {
    for( int i=0; i<mass_particle_vec.size();i++)
    {
      if (mass_particle_vec[i].name == Mass_instance.name)
      {
        mass_particle_vec.erase(mass_particle_vec.begin()+i);
      };
    };

    // for( auto instance:mass_particle_vec)
    // {
    //   if (instance.name == Mass_instance.name)
    //   {
    //     mass_particle_vec.erase(instance);
    //   };
    // };
  };

  void MassiveParticle::calculateAcceleration()
  {
    acc << .0,.0,.0;
    Eigen::Vector3d ri(.0,.0,.0);
    for( auto body_instance:mass_particle_vec)
    {
      ri = position -body_instance.position;
      acc += -body_instance.Mu/(ri.dot(ri))*ri;
    };
  };

  void MassiveParticle::integrateTimestep(double timestep)
  {
    // check acceleration
    acc_not_zero(acc);
    acc_is_const(acc);
    // update position
    Eigen::Vector3d step=timestep*velocity;
    position+=step;
    // update velocity
    acc*=timestep;
    velocity+=acc;
  };



} // end namespace
