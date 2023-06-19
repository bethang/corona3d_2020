/*
 * vtally.cpp
 *
 *  Created on: Jun 12, 2023
 *      Author: bethan
 */
#include "Atmosphere.hpp"
#include "Vtally.hpp"
#include <iostream>

Vtally::Vtally(vector<shared_ptr<Particle>> my_parts, int num_parts, int num_EDFs, vector<int> EDF_alts, double vtally_x, double vtally_dx, double vtally_psi, double vtally_w) {

  num_vel_bins = 90; // set number of velocity bins here
  vel_min = -45.0;   // minimum velocity for velocity distributions (km/s)
  vel_max = 45.0;    // maximum velocity for velocity distributions (km/s)
  dx = 1;            // if particle's x coordinate is between chosen altitude and chosen altitude + dx, add it to the count (km)
  no_LOS_angle_bins = 12;   // no of LOS angles over which vel dist is averaged
  w = 100e5;             // width of line of sight, in y-z plane (cm)
  
  vtally_matrix.resize(num_vel_bins);
  for (int i=0; i<num_vel_bins; i++)
    {
      vtally_matrix[i].resize(num_EDFs);
      for (int j=0; j<num_EDFs; j++)
      {
        vtally_matrix[i][j] = 0.0;
      }
    }

  //std::cout << "num_parts " << num_parts << "\n";
  //  print_things(vtally_x, my_parts);
  
  
  
}

	  Vtally::~Vtally() {
  }

// can probably remove dt
void Vtally::update_vtally(int i, vector<shared_ptr<Particle>> my_parts, int num_EDFs, vector<int> EDF_alts){
  for (int j=0; j<num_EDFs; j++)
    {
      if (is_inside(i, my_parts, EDF_alts[j]))
     	{
	  // 1) calculate average LOS velocity
	    	  LOS_velocity = calculate_LOS_velocity(i, my_parts);
	  //  std::cout << LOS_velocity << "\n";
	  // 2) work out which velocity bin it should go into
       	  vel_bin = choose_vel_bin();
	  //  vel_bin = 5;
	  //  std::cout << vel_bin << "\n";
	  // 3) add count to velocity bin
	  //  vtally_matrix[i][j] += 1;
		      //   std:: cout << j << "\t" << vtally_matrix[1][j] << "\n";
		      //std::cout << EDF_alts[j] << "\t" << is_inside(i, my_parts, EDF_alts[j]) << "\n";
       	}
    }
}

double Vtally::calculate_LOS_velocity(int i, vector<shared_ptr<Particle>> my_parts){
  // std::cout << w << "\t" << no_LOS_angle_bins << "\n";
  double x = my_parts[i]->get_x(); //(cm)
  double y = my_parts[i]->get_y(); //(cm)
  double z = my_parts[i]->get_z(); //(cm)
  double vx = my_parts[i]->get_vx(); //(cm)
  double vy = my_parts[i]->get_vy(); //(cm)
  double vz = my_parts[i]->get_vz(); //(cm)
  angle_spacing = constants::pi/no_LOS_angle_bins;
    for (int l=0; l<no_LOS_angle_bins; l++)
       {
	 angle = l*angle_spacing;
	 // the below is for if it lies exactly along the line along the middle of the line of sight. Adjust to take into account w
	 if (angle == 0.0){
	   if (z >= -1*w/2 and z<= w/2){
	     std::cout << "in angle = 0; x is " << x << ", y is " << y << ", z is " << z << "\n";
	     
	   }
	 }
	 else if (angle == constants::pi/2){
	   if (y>= -1*w/2 and y<= w/2){
	   std::cout << "in angle = pi/2; x is " << x << ", y is " << y << ", z is " << z << "\n";
	   }
         }
	 else {
	   if (y >= (-1*z)/tan(angle) - w/(2*sin(angle)) and y <= (-1*z)/tan(angle) + w/(2*sin(angle))){
	     std::cout << "in angle = " << angle << "; x is " << x << ", y is " << y << ", z is " << z << "\n";
	 }
	 }
	 
	 //    angle_list.push_back(l*angle_spacing);
  // for each angle, find out if the particle is in the line of sight
      // for (theta = 0; theta < 360; theta++)
       }
  return 50.3;
}

int Vtally::choose_vel_bin(){
  //  std::cout << LOS_velocity << "\t" << num_vel_bins << "\n";
  //  std::cout << vel_min << "\t" << vel_max << "\n";
  return 5;
}

void Vtally::print_things(double vtally_x, vector<shared_ptr<Particle>> my_parts){
  std::cout << "vtally_x " << vtally_x << "\n";
  std::cout << "parts " << my_parts[0] << "\n";
}

// find out whether particle i has x co-ordinate at chosen altitude (between chosen altitude and 1 km greater)
bool Vtally::is_inside(int i, vector<shared_ptr<Particle>> my_parts, double alt){  //(int i)
 double x = my_parts[i]->get_x(); // in cm
 
 if (x*1e-5 >= alt and x*1e-5 < alt+dx)
   {
     return true;
   }
 else
   {
     return false;
   }
}
  
