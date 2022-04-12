#include <nbsimMyFunctions.h>
#include <nbsimExceptionMacro.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
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

// double kinetic_energy(std::vector <nbsim::MassiveParticle> body_vect)
// {
//   double E_k=0.0;
//   for (int i=0; i< body_vect.size();i++)
//   {
//     E_k += body_vect[i].Mu*(body_vect[i].getVelocity().squaredNorm());
//   }
//   return E_k/2;
// }

// double potential_energy(std::vector <nbsim::MassiveParticle> body_vec)
// {
//   double E_p=0.0;
//   for (int i=0; i< body_vec.size();i++)
//   {
//     for(int ii=0; ii<body_vec[0].mass_particle_vec.size();ii++)
//     {
//       E_p += -body_vec[i].getMu()*body_vec[i].mass_particle_vec[ii].getMu()/(body_vec[i].getPosition()-body_vec[i].mass_particle_vec[ii].getPosition()).norm();
//     }
//   }
//   return E_p/2;
// }

int main(int argc, char** argv)
{
    // std::clock_t c_start=std::clock();
    // auto t_start = std::chrono::high_resolution_clock::now();
    // // initialize the class
    // nbsim::MassiveParticle Sun,Mercury,Venus,Earth,Mars,Jupiter,Saturn,Uranus,Neptune;
    // std::vector <nbsim::MassiveParticle> body_vec;
    // body_vec.push_back(Sun);
    // body_vec.push_back(Mercury);
    // body_vec.push_back(Venus);
    // body_vec.push_back(Earth);
    // body_vec.push_back(Mars);
    // body_vec.push_back(Jupiter);
    // body_vec.push_back(Saturn);
    // body_vec.push_back(Uranus);
    // body_vec.push_back(Neptune);
    // for(int i=0; i<9;i++)
    // {
    //     body_vec[i].name=nbsim::solarSystemData.at(i).name;
    //     body_vec[i].Mu=nbsim::solarSystemData.at(i).mu;
    //     body_vec[i].position=nbsim::solarSystemData.at(i).position;
    //     body_vec[i].velocity=nbsim::solarSystemData.at(i).velocity;
    // };
    // for(int iii=0; iii<9;iii++)
    // {
    //     for (int ii=0; ii<9;ii++)
    //     {
    //         if(body_vec[iii].name!=body_vec[ii].name)
    //         {
    //             body_vec[iii].addAttractor(body_vec[ii]);
    //         }
    //     }
    // };
    // // 测试所有的body——vec的数量
    // // for(int i=0; i<9;i++)
    // // {
    // //   std::cout<< body_vec[i].name<<": \n"<<std::endl;
    // //   std::cout<< body_vec[i].Mu<<std::endl;  
    // //   std::cout<< body_vec[i].position<<std::endl;  
    // //   std::cout<< body_vec[i].velocity<<std::endl;  
    // //   for (int ii=0; ii<body_vec[ii].mass_particle_vec.size();ii++)
    // //   {
    // //     std::cout<< body_vec[i].mass_particle_vec[ii].name<<std::endl;       

    // //   } 
    // //   std::cout<<"\n"<<std::endl;
    // // };

    // // calcualte initial total energy
    // double K_energy_initial=kinetic_energy(body_vec);
    // double P_energy_initial=potential_energy(body_vec);
    // double Total_energy_initial= K_energy_initial+P_energy_initial;
    // std::cout<< K_energy_initial<<"\n"<< std::endl;
    // std::cout<< P_energy_initial<<"\n"<< std::endl;
    // std::cout<< Total_energy_initial<<"\n"<<std::endl;

    // double timestep= 0.000274;
    // for (int i=0;i<3650;i++)
    // {
    //     for (int ii=0;ii<body_vec.size();ii++)
    //     {
    //         body_vec[ii].calculateAcceleration();
    //         // std::cout<<body_vec[ii].name<<" acc are: \n";
    //         // std::cout<<body_vec[ii].acc<<std::endl;

    //     };

    //     for (int iii=0;iii<body_vec.size();iii++)
    //     {
    //         body_vec[iii].integrateTimestep(timestep);
    //         // std::cout<<"positions are: \n";
    //         // std::cout<<body_vec[iii].position<<std::endl;

    //     };
    //     // std::cout<<body_vec[0].position<<std::endl;
    // };
    // // end of the clock
    // std::clock_t c_end = std::clock();
    // auto t_end = std::chrono::high_resolution_clock::now();
    // std::cout << std::fixed << std::setprecision(2) << "CPU time used: "
    //           << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
    //           << "Wall clock time passed: "
    //           << std::chrono::duration<double, std::milli>(t_end-t_start).count()
    //           << " ms\n\n";

    // std::cout<<"Final positions are: \n";
    // for (int i=0;i<body_vec.size();i++)
    // {
    //     std::cout<<body_vec[i].name<<": \n";
    //     std::cout<<body_vec[i].position<<"\n\n";
    //     // std::cout<<body_vec[i].position.transpose()<<"\n\n";
    // }

    // // calcualte total energy
    // double K_energy_end=kinetic_energy(body_vec);
    // double P_energy_end=potential_energy(body_vec);
    // double Total_energy_end= K_energy_end+P_energy_end;
    // std::cout<< K_energy_end<<"\n"<< std::endl;
    // std::cout<< P_energy_end<<"\n"<< std::endl;
    // std::cout<< Total_energy_end<<"\n"<< std::endl;
    return 0;
}
