/*
 * phys.cpp
 *
 *  Created on: May 29, 2020
 *      Author: rodney
 */

#include "phys.h"
#include "constants.h"
#include "Particle.h"
#include "Planet.h"
#include <cstdlib>
#include <cmath>

void initPlanet(Planet plnt)
{

}

void initParticles(Planet plnt, Particle parts[], int n)
{
	double scaleHeight = 20e3; // should calc/assign in planet class eventually
	for (int i=0; i<n; i++)
	{
		double r = plnt.radius + 160e3 - log(getRand())*scaleHeight;
		double phi = twopi*(getRand());
		double u = 2.0*getRand() - 1;
		parts[i].radius = r;
		parts[i].inverseRadius = 1.0/r;
		parts[i].position[0] = r*sqrt(1-(u*u))*cos(phi);
		parts[i].position[1] = r*sqrt(1-(u*u))*sin(phi);
		parts[i].position[2] = r*u;
	}
}

// Iterate Equation of Motion using Velocity-Verlet Algorithm
void stepParticles(Particle parts[], double dt, int n, double k_g)
{
	double a[3] = {0.0, 0.0, 0.0};
	for (int i=0; i<n; i++)
	{
		if (parts[i].active == false)
		{
			continue;
		}
		else
		{
			// calculate acceleration at current position
			double invRcube = pow(parts[i].inverseRadius, 3);
			a[0] = k_g*parts[i].position[0]*invRcube;
			a[1] = k_g*parts[i].position[1]*invRcube;
			a[2] = k_g*parts[i].position[2]*invRcube;

			// calculate next position and update particle
			parts[i].position[0] = parts[i].position[0] + parts[i].velocity[0]*dt + 0.5*a[0]*dt*dt;
			parts[i].position[1] = parts[i].position[1] + parts[i].velocity[1]*dt + 0.5*a[1]*dt*dt;
			parts[i].position[2] = parts[i].position[2] + parts[i].velocity[2]*dt + 0.5*a[2]*dt*dt;
			parts[i].radius = sqrt(pow(parts[i].position[0], 2) + pow(parts[i].position[1], 2) + pow(parts[i].position[2], 2));
			parts[i].inverseRadius = 1 / parts[i].radius;

			// calculate acceleration at next position
			invRcube = pow(parts[i].inverseRadius, 3);
			a[0] = k_g*parts[i].position[0]*invRcube;
			a[1] = k_g*parts[i].position[1]*invRcube;
			a[2] = k_g*parts[i].position[2]*invRcube;

			// calculate next velocity using acceleration at next position and update particle
			parts[i].velocity[0] = parts[i].velocity[0] + 0.5*a[0]*dt;
			parts[i].velocity[1] = parts[i].velocity[1] + 0.5*a[1]*dt;
			parts[i].velocity[2] = parts[i].velocity[2] + 0.5*a[2]*dt;
		}
	}
}

double getRand()
{
	return ((double)rand() / RAND_MAX);
}