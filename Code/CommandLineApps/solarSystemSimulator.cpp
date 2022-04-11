#include <nbsimMyFunctions.h>
#include <nbsimExceptionMacro.h>
#include <iostream>
#include <stdlib.h>
#include <vector>

// Example, header-only library, included in project for simplicity's sake.
#include <Eigen/Dense>
#include <unsupported/Eigen/NonLinearOptimization>
// reading the SolarSystemData
#include "nbsimSolarSystemData.ipp"
// 
#include "nbsimBasicTypes.h"

int main(int argc, char** argv)
{
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
/*     for(int i=0; i<9;i++)
    {
    std::cout<< body_vec[i].mass_particle_vec.size()<<std::endl;        
    } */

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

    // output the position
    std::cout<<"Final positions are: \n";
    for (int i=0;i<body_vec.size();i++)
    {
        std::cout<<body_vec[i].name<<": \n";
        std::cout<<body_vec[i].getPosition()<<"\n\n";
    }
    return 0;
}