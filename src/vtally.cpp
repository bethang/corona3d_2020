/*
 * vtally.cpp
 *
 *  Created on: Jun 12, 2023
 *      Author: bethan
 */
#include "Atmosphere.hpp"
#include "Vtally.hpp"
#include <iostream>

Vtally::Vtally(int num_parts, int num_EDFs, vector<int> EDF_alts, double vtally_x, double vtally_dx, double vtally_psi, double vtally_w, double dt){//, double rate) {

  std::cout << "num_parts\t" << num_parts;
  num_vel_bins = 36; // set number of velocity bins here
  vel_min = -45.0;   // minimum velocity for velocity distributions (km/s)
  vel_max = 45.0;    // maximum velocity for velocity distributions (km/s)
  dx = 1;            // if particle's x coordinate is between chosen altitude and chosen altitude + dx, add it to the count (km)
  no_LOS_angle_bins = 12;   // no of LOS angles over which vel dist is averaged
  angle_spacing = constants::pi/no_LOS_angle_bins;
  w = 10000e5;             // width of line of sight, in y-z plane (cm)

  vtally_matrix.resize(num_vel_bins);
  for (int i=0; i<num_vel_bins; i++)
    {
      vtally_matrix[i].resize(num_EDFs);
      for (int j=0; j<num_EDFs; j++)
      {
        vtally_matrix[i][j] = 0.0;
      }
    }
  cos_sin_angle.resize(2); // array for cos(angle) and sin(angle) values
  for (int i=0; i<2; i++){
    cos_sin_angle[i].resize(no_LOS_angle_bins);
  }
  for (int j=0; j<no_LOS_angle_bins; j++){
    cos_sin_angle[0][j] = cos(j*angle_spacing);
    cos_sin_angle[1][j] = sin(j*angle_spacing);
   }
  
  
}

	  Vtally::~Vtally() {
  }

// this is not currently used
void Vtally::update_vtally(int i, const vector<shared_ptr<Particle>> &my_parts, int num_EDFs, vector<int> EDF_alts, double dt, double rate, int num_parts){
  for (int j=0; j<num_EDFs; j++)
    {
      if (is_inside(i, my_parts, EDF_alts[j]))
     	{
	  // 1) calculate average LOS velocity
	  calculate_LOS_velocity(i, my_parts, j, dt, rate, num_parts);
       	}
    }
}

void Vtally::calculate_LOS_velocity(int i, vector<shared_ptr<Particle>> my_parts, int j, double dt, double rate, int num_parts){ //, int l, double angle){
  double x = my_parts[i]->get_x(); //(cm)
  double y = my_parts[i]->get_y(); //(cm)
  double z = my_parts[i]->get_z(); //(cm)
  double vx = my_parts[i]->get_vx(); //(cm)
  double vy = my_parts[i]->get_vy(); //(cm)
  double vz = my_parts[i]->get_vz(); //(cm)

  weighting_factor = dt*rate/(no_LOS_angle_bins*num_parts*(dx*1e5)*w); // contribution of each count to column density

  for (int l=0; l<no_LOS_angle_bins; l++){
  	 angle = l*angle_spacing;
	 if (angle == 0.0){
	   if (z >= -1*w/2 and z<= w/2){
	     //    std::cout << "in angle = 0; x is " << x << ", y is " << y << ", z is " << z << ", vy is " << vy << "\n";
	          vel_bin = choose_vel_bin(vy);
	          vtally_matrix[vel_bin][j] += weighting_factor;
	          vel_bin = choose_vel_bin(-1*vy);
		  vtally_matrix[vel_bin][j] += weighting_factor;
	     
	   }
	 }
	 else if (angle == constants::pi/2){
	   if (y>= -1*w/2 and y<= w/2){
	     //     std::cout << "in angle = pi/2; x is " << x << ", y is " << y << ", z is " << z << "\n";
	     vel_bin = choose_vel_bin(vz);
	     vtally_matrix[vel_bin][j] += weighting_factor;
	     vel_bin = choose_vel_bin(-1*vz);
	     vtally_matrix[vel_bin][j] += weighting_factor;
	     
	   }
         }
	 else {
	   if (y >= (-1*z)/tan(angle) - w/(2*sin(angle)) and y <= (-1*z)/tan(angle) + w/(2*sin(angle))){
	     vel_LOS = vy*cos_sin_angle[0][l] - vz*cos_sin_angle[1][l]; // this is equal to the component of vy along the LOS + the component of vz along the LOS = vy*cos(angle) = vz*sin(angle)
	     vel_bin = choose_vel_bin(vel_LOS);
	     vtally_matrix[vel_bin][j] += weighting_factor;
	     vel_bin = choose_vel_bin(-1*vel_LOS);
	     vtally_matrix[vel_bin][j] += weighting_factor;
	 }
	 }
	 
	 //    angle_list.push_back(l*angle_spacing);
  // for each angle, find out if the particle is in the line of sight
      // for (theta = 0; theta < 360; theta++)
	 //      }
  }
  } // NOTE THAT, AS IS, IT WILL SOMETIMES/OFTEN RETURN NOTHING!!

int Vtally::choose_vel_bin(const double &LOS_velocity) {
  bin = std::floor(((LOS_velocity*1e-5)-vel_min)/2.5);
  return bin;
}

void Vtally::print_things(double vtally_x, vector<shared_ptr<Particle>> my_parts){
  std::cout << "vtally_x " << vtally_x << "\n";
  std::cout << "parts " << my_parts[0] << "\n";
}

// find out whether particle i has x co-ordinate at chosen altitude (between chosen altitude and 1 km greater)
bool Vtally::is_inside(int i, vector<shared_ptr<Particle>> my_parts, double alt){  //(int i)
 const double x = my_parts[i]->get_x(); // in cm
 
 if (x*1e-5 >= alt and x*1e-5 < alt+dx)
   {
     return true;
   }
 else
   {
     return false;
   }
}

void Vtally::record_vtallies(vector<int> stats_EDF_alts, int i, string output_dir){
  ofstream vtally_out;
  vtally_out.open(output_dir + "vtally_" + to_string(stats_EDF_alts[i]) + "km.csv");
  vtally_out << "#velocity_bin_no,min_velocity_in_bin_[km/s],column_density_[cm-2]\n";
  for (int j=0; j < num_vel_bins; j++){
    vtally_out << j << "," << vel_min + j*((vel_max-vel_min)/num_vel_bins) << "," << vtally_matrix[j][i] <<"\n"; // here, j is the vel_bin count and i is the altitude count
  }
}

  
