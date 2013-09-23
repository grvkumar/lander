// Mars lander simulator
// Version 1.7
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, March 2013

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

//includes
#include "lander.h"
#include <stdio.h>
#include <iostream>
#include <cmath>

#define VERLET //use the verlet integrator, alternatively define EULER

// variables omega and dist_geostat have been created to store the angular velocity and the initial distance for the geostationary orbit
double omega = 2*PI/MARS_DAY;
double dist_geostat = pow(GRAVITY*MARS_MASS/(omega*omega), (double) 1/3);

using namespace std;

void autopilot (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
	//define the control constants upfront
	const double K_H = 0.065;
	const double K_P = 1.0;
	const double DELTA = 0.3;

	double h = position.abs() - MARS_RADIUS; //height of the lander from surface
	double e = -1.0*(0.5 + K_H*h + velocity*position.norm()); //the error term

	double p_out = K_P*e;

	if(p_out <= -1.0*DELTA)
	{
		throttle = 0.0;
	}
	else if (p_out < 1.0 - DELTA)
	{
		throttle = p_out + DELTA;
	}
	else
	{
		throttle = 1.0;
	}	

	if( e > 1.0 && safe_to_deploy_parachute())
	{
		parachute_status = DEPLOYED;
	}
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{

	// List of forces to consider for accelerations
	vector3d grav_attraction;
	vector3d thrust;
	vector3d atm_drag;
	vector3d chute_drag;
	vector3d net_acc;

    // Mass of lander
    double mass_lander = UNLOADED_LANDER_MASS - (1-fuel)*FUEL_CAPACITY*FUEL_DENSITY;

    // Calculate acceleration due to gravitional attraction
    grav_attraction = -1.0*((GRAVITY*MARS_MASS)/position.abs2())*position.norm();

    // Calculate acceleration due to thrust
    thrust = thrust_wrt_world()/mass_lander;

    // Calculate the acceleration due to atmospheric drag
    atm_drag = -0.25*((atmospheric_density(position)*DRAG_COEF_LANDER*LANDER_SIZE*LANDER_SIZE*PI*velocity.abs2())*velocity.norm())/mass_lander;


    // Calculate the acceleration due to drag due to parachutes
    if(parachute_status == DEPLOYED)
    {
        chute_drag = -10.0*((atmospheric_density(position)*DRAG_COEF_CHUTE*LANDER_SIZE*LANDER_SIZE*velocity.abs2())*velocity.norm())/mass_lander;
    }
    else
    {
        chute_drag = vector3d(0.0,0.0,0.0);
    }

    // Net acceleration

    net_acc = grav_attraction + thrust + atm_drag + chute_drag;

#ifdef VERLET
	// Update position and velocity using the Euler algorithm
    static vector3d prev_position(position), prev_velocity(velocity), prev_acc(net_acc);
    vector3d temp = position;
    vector3d temp_v = velocity;
    vector3d temp_a = net_acc;

    if(simulation_time == 0.0)
    {
    	position = position + velocity*delta_t + 0.5*net_acc*delta_t*delta_t;
    	velocity = velocity + net_acc*delta_t;
    }

    else
    {
    	position = 2*position - prev_position + delta_t*delta_t*net_acc;
    	velocity = 2*velocity - prev_velocity + delta_t*(net_acc - prev_acc);
    }

    prev_position = temp;
    prev_velocity = temp_v;
    prev_acc = temp_a;
#endif

#ifdef EULER
    position = position + velocity*delta_t + 0.5*net_acc*delta_t*delta_t;
    velocity = velocity + net_acc*delta_t;
#endif




  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled) autopilot();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "geostationary orbit";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;


    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;


    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;


    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;


    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;


    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;


    break;

  case 6:
	// a geostationary orbit

	position = vector3d(dist_geostat, 0.0, 0.0);
    velocity = vector3d(0.0, dist_geostat*omega, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.5;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
