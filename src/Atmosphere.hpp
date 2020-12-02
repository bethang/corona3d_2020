/*
 * Atmosphere.hpp
 *
 *  Created on: Jun 8, 2020
 *      Author: rodney
 */

#ifndef ATMOSPHERE_HPP_
#define ATMOSPHERE_HPP_

#include <vector>
#include <iomanip>
#include "Background_Species.hpp"
#include "Distribution_Hot_H.hpp"
#include "Distribution_Hot_O.hpp"
#include "Distribution_Import.hpp"
#include "Distribution_MB.hpp"
#include "Common_Functions.hpp"
using namespace std;

class Atmosphere {
public:
	Atmosphere(int n, int num_to_trace, Planet p, vector<shared_ptr<Particle>> parts, shared_ptr<Distribution> dist, Background_Species bg, int pos_out_freq, string pos_out_dir);
	virtual ~Atmosphere();

	void init_shell(double bottom_r, double top_r, int num_bins, double bin_width, string output_dir);
	void output_positions(string datapath);
	void output_altitude_distro(double bin_width, string datapath);
	void output_velocity_distro(double bin_width, string datapath);
	void run_simulation(double dt, int num_steps);

private:
	int num_parts;                      // number of particles initially spawned
	int num_traced;                     // number of particles to output trace data on
	int active_parts;                   // number of active particles
	Planet my_planet;                   // contains planet mass and radius
	vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
	shared_ptr<Distribution> my_dist;              // distribution class to initialize particles
	Background_Species bg_species;      // background species used for collisions
	vector<int> traced_parts;           // indices of randomly selected trace particles
	int output_pos_freq;                // number of timesteps between outputting all active particle positions
	string output_pos_dir;              // directory to output active particles positions to

	bool shell_active;
	double shell_bottom;
	double shell_top;
	int shell_enter_top;
	int shell_enter_bottom;
	int shell_exit_top;
	int shell_exit_bottom;
	int shell_numvelbins;
	int shell_velbinwidth;
	vector<int> shell_velbins;
	string shell_output_dir;

	vector<int> stats_alt_bins;
	vector<int> stats_dens_counts;
	void update_stats();
	void output_stats();

	// output test particle trace data for selected particles
	void output_collision_data();
	void output_trace_data();

	// output accumulated shell data
	void output_shell_data();

	// update shell numbers (velocity distro, etc.)
	void update_shell_data();
};

#endif /* ATMOSPHERE_HPP_ */
