/*
 * Distribution_MB.hpp
 *
 *  Created on: Jul 21, 2020
 *      Author: rodney
 */

#ifndef DISTRIBUTION_MB_HPP_
#define DISTRIBUTION_MB_HPP_

#include "Distribution.hpp"

class Distribution_MB: public Distribution {
public:
	Distribution_MB(Planet my_p, double ref_h, double ref_T);
	virtual ~Distribution_MB();
	void init(Particle* p);
	void init_vonly(Particle* p, double v_avg);

};

#endif /* DISTRIBUTION_MB_HPP_ */
