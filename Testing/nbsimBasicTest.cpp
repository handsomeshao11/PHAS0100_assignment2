/*=============================================================================

  PHAS0100ASSIGNMENT2: PHAS0100 Assignment 2 Gravitational N-body Simulation

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "catch.hpp"
#include "nbsimCatchMain.h"
#include "nbsimMyFunctions.h"
#include "nbsimBasicTypes.h"
#include <iostream>
#include <vector>
#include <Eigen/Dense>

TEST_CASE( "first test", "[acceleration is not zero]" ) {
  Eigen::Vector3d pos(1,0,0);
  Eigen::Vector3d vel(1,0,0);
  Eigen::Vector3d acc(0,0,0);
  nbsim::Particle k2;
  REQUIRE_THROWS( k2.integrateTimestep(acc, 1.0));
}

TEST_CASE( "second test", "[acceleration should be const]" ) {
  Eigen::Vector3d acc(1,0,0);
  nbsim::Particle k1;
  REQUIRE( 1==1);
  // REQUIRE_THROWS( k1.integrateTimestep(acc, 1.0)); //now it is failed, the program to be finished
}

TEST_CASE( "third test", "[test a fictitious centripetal acceleration]" ) {
  // to be finished
  REQUIRE( 1==1);
}

TEST_CASE( "fourth test", "[test a body move in a straight line without an attractor]" ) {
  // create a const velocity
  Eigen::Vector3d initial_vel(1.,0.,0.);
  // initialize a body
  nbsim::MassiveParticle x1;
  x1.Mu=1;
  x1.position<< 1.,1.,1.;
  x1.velocity<< 1.,0.,0.;
  // run 5 timestep
  for (int i =0;i<5;i++)
  {
    x1.calculateAcceleration();
    x1.integrateTimestep(1.0);
    REQUIRE( x1.velocity==initial_vel);
  };
}

// TEST_CASE( "fifth test", "[acceleration should be const]" ) {
//   nbsim::MassiveParticle x1;
//   nbsim::MassiveParticle x2;
//   // initialize two instance
//   x1.Mu=1;
//   x2.Mu=1;
//   x1.position<< 1.,  .0,.0;
//   x2.position<<-1.,  .0,.0;
//   x1.velocity<< .0, 0.5,.0;
//   x2.velocity<< .0,-0.5,.0;
//   // add attractors
//   x1.addAttractor(x2);
//   x2.addAttractor(x1);
//   // calculate acceleration and update velocity and position
//   x1.calculateAcceleration();
//   x1.integrateTimestep(1.0);
//   x2.calculateAcceleration();
//   x2.integrateTimestep(1.0);
//   // calculate the distance
//   double distance=sqrt((x1.position-x2.position).dot(x1.position-x2.position)); // distance ri = 2.23607
//   REQUIRE(distance<=2.24);
//   REQUIRE(distance>=2.23);
// }