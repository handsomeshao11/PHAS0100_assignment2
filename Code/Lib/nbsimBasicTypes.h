/*=============================================================================

  PHAS0100ASSIGNMENT2: PHAS0100 Assignment 2 Gravitational N-body Simulation

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#ifndef nbsimBasicTypes_h
#define nbsimBasicTypes_h

// /**
// * \defgroup internal internal
// * \brief Internal stuff, not for end-users.
// */

// /**
// * \defgroup types types
// * \brief Package-wide data types.
// */

// /**
// * \defgroup utilities utilities
// * \brief Groups of c-style functions.
// */

// /**
// * \defgroup applications applications
// * \brief Small, end-user applications, most likely command line.
// */

// /**
// * \file nbsimBasicTypes.h
// * \brief Defines types and typedefs used in this library.
// * \ingroup types
// */

// //! Single namespace for all code in this package
#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
#include <stdlib.h>
#include <vector>
#include <set>
#include <memory>
namespace nbsim
{
  void acc_not_zero(Eigen::Vector3d acceleration); // check acceleration is not 0
  void acc_is_const(Eigen::Vector3d acceleration); // check acceleration is const

  

  class Particle
  {
    public:
    Particle(){};
    // ~Particle();
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    Eigen::Vector3d acc;
    Eigen::Vector3d getPosition(); // return position
    Eigen::Vector3d getVelocity(); // return velocity
    void integrateTimestep(Eigen::Vector3d acceleration, double timestep);

  };

  class MassiveParticle : public Particle 
  {
    private:
    double G = 6.6743e-11; // units m3 kg-1 s-2 

    public:
    MassiveParticle(){};
    MassiveParticle(std::string init_name,Eigen::Vector3d init_position,Eigen::Vector3d init_velocity,double init_Mu);
    // ~MassiveParticle();
    std::string name;
    double Mu;
    double getMu(); 
    std::vector <MassiveParticle> mass_particle_vec; // A vector to store MassiveParticle instance
    std::set<std::shared_ptr<MassiveParticle>> attractors_ptr;
    // void addAttractor(MassiveParticle Mass_instance);
    // void removeAttractor(MassiveParticle Mass_instance);
    void addAttractor(std::shared_ptr<MassiveParticle> attractor_ptr);
    void removeAttractor(std::shared_ptr<MassiveParticle> attractor_ptr);
    void calculateAcceleration();
    void integrateTimestep(double timestep);
  };


} // end namespace

#endif
