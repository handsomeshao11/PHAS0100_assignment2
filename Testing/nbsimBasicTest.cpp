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
#include <math.h>
#define _USE_MATH_DEFINES

TEST_CASE( "first test", "[acceleration is not zero]" ) {
  Eigen::Vector3d pos(1,0,0);
  Eigen::Vector3d vel(1,0,0);
  Eigen::Vector3d acc(0,0,0);
  nbsim::Particle k2;
  REQUIRE_THROWS( k2.integrateTimestep(acc, 1.0));
}

TEST_CASE( "second test", "[test move with const acceleration]" ) {
  Eigen::Vector3d acc(1,0,0);
  Eigen::Vector3d pos(0,1,0);
  Eigen::Vector3d vel(0,0,0);
  nbsim::Particle k1;
  k1.position=pos;
  k1.velocity=vel;
  k1.integrateTimestep(acc,1.0);
  pos+=1.0*vel;
  vel+=1.0*acc;
  REQUIRE(k1.getPosition().isApprox(pos,0.01));
  REQUIRE(k1.getVelocity().isApprox(vel,0.01));
  k1.integrateTimestep(acc,1.0);
  pos+=1.0*vel;
  vel+=1.0*acc;
  REQUIRE(k1.getPosition().isApprox(pos,0.01));
  REQUIRE(k1.getVelocity().isApprox(vel,0.01));
}

TEST_CASE( "third test", "[test move with a fictitious centripetal acceleration]" ) {
  Eigen::Vector3d pos(1,0,0),vel(0,1,0);
  nbsim::Particle k3;
  k3.position=pos;
  k3.velocity=vel;
  double time_interval=0.001, end_time=2*M_PI;
  Eigen::Vector3d acc;
  for(double i=0;i<end_time;i+=time_interval)
  {
    acc=-k3.getPosition();
    k3.integrateTimestep(acc,time_interval);
  }
  REQUIRE((k3.getVelocity().isApprox(vel,0.01) && k3.getPosition().isApprox(pos,0.01))==1);
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

TEST_CASE( "fifth test", "[two bodies remain at a distance]" ) {
  Eigen::Vector3d x1(1.,0.,0.),v1(0.,0.5,0.),x2(-1.,0.,0.),v2(0.,-0.5,0.);
  std::shared_ptr<nbsim::MassiveParticle> first(new nbsim::MassiveParticle("-",x1,v1,1.0));
  std::shared_ptr<nbsim::MassiveParticle> second(new nbsim::MassiveParticle("-",x2,v2,1.0));
  first->addAttractor(second);
  second->addAttractor(first);
  for (int i=0;i<4;i++)
  {
    first->calculateAcceleration();
    second->calculateAcceleration();
    first->integrateTimestep(0.5);
    second->integrateTimestep(0.5);
    REQUIRE((first->getPosition()-second->getPosition()).norm()<2.3);
    REQUIRE((first->getPosition()-second->getPosition()).norm()>1.7);
  }
}