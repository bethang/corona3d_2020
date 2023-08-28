/*
 * vtally.cpp
 *
 *  Created on: Jun 12, 2023
 *      Author: bethan
 */
#include "Atmosphere.hpp"
#include "Vtally.hpp"
#include <iostream>

Vtally::Vtally(const int &num_EDFs, const vector<int> &EDF_alts, const double &vtally_x, const double &vtally_dx,  const double &vtally_w, const double &dt, const double &rate, const int &num_parts, const double &radius){ // 

  num_vel_bins = 60; // set number of velocity bins here
  vel_min = -45.0;   // minimum velocity for velocity distributions (km/s)
  vel_max = 45.0;    // maximum velocity for velocity distributions (km/s)
  dx = 100;            // if particle's x coordinate is between chosen altitude and chosen altitude + dx, add it to the count (km)
  no_LOS_angle_bins = 12;   // no of LOS angles over which vel dist is averaged
  angle_spacing = constants::pi/no_LOS_angle_bins;
  w = 100e5;             // width of line of sight, in y-z plane (cm)
  weighting_factor = dt*rate/(no_LOS_angle_bins*num_parts*(dx*1e5)*w);
  radius_km = radius/1e5;

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
void Vtally::update_vtally(const shared_ptr<Particle> &p, const int &num_EDFs, const vector<int> &EDF_alts){
  for (int j=0; j<num_EDFs; j++)
    {
      if (is_inside(p, EDF_alts[j]))
     	{
	  // calculate average LOS velocity
	  calculate_LOS_velocity(p, j);
       	}
    }
}

void Vtally::calculate_LOS_velocity(const shared_ptr<Particle> &p, const int &j){
  const double y = p->get_y(); //(cm)
  const double z = p->get_z(); //(cm)
  const double vy = p->get_vy(); //(cm)
  const double vz = p->get_vz(); //(cm)
  
  for (int l=0; l<no_LOS_angle_bins; l++){ // make this more efficient - probably don't need this loop... maybe? Probably a more efficient way
  	 angle = l*angle_spacing;
	 if (angle == 0.0){
	   if (z >= -1*w/2 and z<= w/2){
		  vtally_matrix[choose_vel_bin(vy)][j] += weighting_factor;
		  vtally_matrix[choose_vel_bin(-1*vy)][j] += weighting_factor;
	     
	   }
	 }
	 else if (angle == constants::pi/2){
	   if (y>= -1*w/2 and y<= w/2){
	     vtally_matrix[choose_vel_bin(vz)][j] += weighting_factor;
	     vtally_matrix[choose_vel_bin(-1*vz)][j] += weighting_factor;
	     
	   }
         }
	 else {
	   if (y >= (-1*z)/tan(angle) - w/(2*sin(angle)) and y <= (-1*z)/tan(angle) + w/(2*sin(angle))){
	     vel_LOS = vy*cos_sin_angle[0][l] - vz*cos_sin_angle[1][l]; // this is equal to the component of vy along the LOS + the component of vz along the LOS = vy*cos(angle) = vz*sin(angle)
	     vtally_matrix[choose_vel_bin(vel_LOS)][j] += weighting_factor;
	     vtally_matrix[choose_vel_bin(-1*vel_LOS)][j] += weighting_factor;
	 }
	 }
	 
  }
  }

int Vtally::choose_vel_bin(const double &LOS_velocity) {
  bin = std::floor(((LOS_velocity*1e-5)-vel_min)/((vel_max-vel_min)/num_vel_bins));
  return bin;
}


// find out whether particle i has x co-ordinate at chosen altitude (between chosen altitude and 1 km greater)
bool Vtally::is_inside(const shared_ptr<Particle> &p, const double &alt){  //(int i)
 const double alt_p = p->get_x()*1e-5 - radius_km; // in km

 
 if (alt_p >= alt and alt_p < alt+dx)
   {
     return true;
   }
 else
   {
     return false;
   }
}

void Vtally::record_vtallies(const vector<int> &stats_EDF_alts, const int &i, const string &output_dir){
  ofstream vtally_out;
  vtally_out.open(output_dir + "vtally_" + to_string(stats_EDF_alts[i]) + "km.csv");
  vtally_out << "#velocity_bin_no,min_velocity_in_bin_[km/s],column_density_[cm-2]\n";
  for (int j=0; j < num_vel_bins; j++){
    vtally_out << j << "," << vel_min + j*((vel_max-vel_min)/num_vel_bins) << "," << vtally_matrix[j][i] <<"\n"; // here, j is the vel_bin count and i is the altitude count
  }
}

  
