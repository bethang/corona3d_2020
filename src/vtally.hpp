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
  
  
public:
  Vtally(const int &num_EDFs, const vector<int> &EDF_alts, const double &vtally_x, const double &vtally_dx,  const double &vtally_w, const double &dt, const double &rate, const int &num_parts, const double &radius); //, double rate);
  virtual ~Vtally();

  
  // some of the definitions from a different file, for reference
  //    void tallying(double vtally_x, double vtally_dx, double vtally_psi, double vtally_w);
  //     vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
  //	static const double mass;
  //	static const string name;
  //	double get_mass() const;
  //	string get_name() const;

  void update_vtally(const shared_ptr<Particle> &p, const int &num_EDFS, const vector<int> &EDF_alts);
  void calculate_LOS_velocity(const shared_ptr<Particle> &p, const int &j);
  int choose_vel_bin(const double &vel);
  bool is_inside(const shared_ptr<Particle> &p, const double &alt);
  void record_vtallies(const vector<int> &stats_EDF_alts, const int &i, const string &output_dir);

private:
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
       double radius_km;
       vector<vector<double>> vtally_matrix; // rows = no. of velocity bins; columns = no. of altitudes (from corona3d.cfg); contains counts and is added to at each timestep.
       vector<vector<double>> cos_sin_angle; // rows: [0] = cos(angle); [1] = sin(angle); columns: one for each angle

};

#endif /* VTALLY_HPP_ */
