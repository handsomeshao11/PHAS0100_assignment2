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

  // class
  Eigen::Vector3d particle::getPosition()
  {
    return position;
  };

  Eigen::Vector3d particle::getVelocity()
  {
    return velocity;
  };
  void particle::integrateTimestep(Eigen::Vector3d acceleration, double timestep)
  {
    // check acceleration is not 0
    acc_not_zero(acceleration);
    // update position
    velocity*=timestep;
    position+=velocity;

    velocity/=timestep;
    // update velocity
    acceleration*=timestep;
    velocity+=acceleration;
  };




} // end namespace
