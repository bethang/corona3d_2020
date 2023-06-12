/*
 * vtally.hpp
 *
 *  Created on: Jun 19, 2020
 *      Author: rodney
 */

#ifndef VTALLY_HPP_
#define VTALLY_HPP_

//#include "Atmosphere.hpp"
#include "Particle.hpp"

class Vtally { // public Atmosphere {
public:
  Vtally(vector<shared_ptr<Particle>> my_parts,double vtally_x, double vtally_dx, double vtally_psi, double vtally_w);
  virtual ~Vtally();
  //    void tallying(double vtally_x, double vtally_dx, double vtally_psi, double vtally_w);
        vector<shared_ptr<Particle>> my_parts;         // array of particles to be tracked
  //	static const double mass;
  //	static const string name;
  //	double get_mass() const;
  //	string get_name() const;
        double vtally_x;
        double vtally_dx;
        double vtally_psi;
        double vtally_w;

private:
        double vtally_arr; //????
  //    bool is_inside(particle) const;
  //    void tally(particle);
};

#endif /* VTALLY_HPP_ */
