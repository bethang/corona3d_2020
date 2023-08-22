/*
 * vtally.hpp
 *
 *  Created on: Jun 19, 2020
 *      Author: rodney
 */

#ifndef VTALLY_HPP_
#define VTALLY_HPP_

#include "Atmosphere.hpp"
#include "Particle.hpp"

class Vtally
{
private:
  
  
public:
  Vtally(vector<shared_ptr<Particle>> my_parts, int num_parts, int num_EDFs,vector<int> EDF_alts, double vtally_x, double vtally_dx, double vtally_psi, double vtally_w, double dt); //, double rate);
  virtual ~Vtally();
       int num_vel_bins;
       double LOS_velocity_spare;
       int vel_bin;
       double vel_min;
       double vel_max;
       double dx;
       int no_LOS_angle_bins;
       double angle_spacing;
       double angle;
       double w;
       int bin;
       double vel_LOS;
       double weighting_factor;
  
  // some of the definitions from a different file, for reference
  //    void tallying(double vtally_x, double vtally_dx, double vtally_psi, double vtally_w);
  //     vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
  //	static const double mass;
  //	static const string name;
  //	double get_mass() const;
  //	string get_name() const;

  void update_vtally(int i, const vector<shared_ptr<Particle>> &my_parts, int num_EDFS, vector<int> EDF_alts, double dt, double rate, int num_parts);
  void calculate_LOS_velocity(int i, vector<shared_ptr<Particle>> my_parts, int j, double dt, double rate, int num_parts);
  int choose_vel_bin(const double &vel);
  void print_things(double vtally_x, vector<shared_ptr<Particle>> my_parts);
  bool is_inside(int i, vector<shared_ptr<Particle>> my_parts, double alt);
  void record_vtallies(vector<int> stats_EDF_alts, int i, string output_dir);

private:
        vector<vector<double>> vtally_matrix; // rows = no. of velocity bins; columns = no. of altitudes (from corona3d.cfg); contains counts and is added to at each timestep.
        vector<vector<double>> cos_sin_angle; // rows: [0] = cos(angle); [1] = sin(angle); columns: one for each angle

};

#endif /* VTALLY_HPP_ */
