/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 10;
  
  std::vector<double> particle_weights(num_particles,1.0);
  weights = particle_weights;
  
  // Initialize particles
    // Random Number generator
  default_random_engine gen;
  std::normal_distribution<double> distribution_x(x,std[0]);
  std::normal_distribution<double> distribution_y(y,std[1]);
  std::normal_distribution<double> distribution_theta(theta,std[2]);
  
  for (int i = 0;i < num_particles; i++){
    Particle p;
    p.x = distribution_x(gen);
    p.y = distribution_y(gen);
    p.theta = distribution_theta(gen);
    // Add to particles vector
    particles.push_back(p);
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  
  default_random_engine gen;
  std::normal_distribution<double> distribution_x(0,std_pos[0]);
  std::normal_distribution<double> distribution_y(0,std_pos[1]);
  std::normal_distribution<double> distribution_theta(0,std_pos[2]);
  
  for (vector<Particle>::iterator p_it = particles.begin(); p_it != particles.end(); p_it++) {
    // Apply motion model
    double x1,y1,theta1;
    double x0 = p_it->x;
    double y0 = p_it->y;
    double theta0 = p_it->theta;
    if (yaw_rate == 0){
      x1 = x0 + delta_t*velocity*cos(theta0);
      y1 = y0 + delta_t*velocity*sin(theta0);
    } else {
      x1 = x0 + (velocity/yaw_rate)*(sin(theta0 + yaw_rate*delta_t) - sin(theta0));
      y1 = y0 + (velocity/yaw_rate)*(cos(theta0) - cos(theta0 + yaw_rate*delta_t));
    }
    theta1 = theta0 + yaw_rate*delta_t;
    
    p_it->x = x1 + distribution_x(gen);
    p_it->y = y1 + distribution_y(gen);
    p_it->theta = theta1  + distribution_theta(gen);
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  
  // For each observation, find the map landmark with the shortest distance
  for (std::vector<LandmarkObs>::iterator obs_it = observations.begin();
       obs_it!=observations.end(); obs_it++){
    double d_min = INFINITY; // Min dist
    double d = 0;
    
    for (std::vector<LandmarkObs>::iterator pred_it = predicted.begin();
         pred_it!=predicted.end(); pred_it++){
      d = dist(obs_it->x, obs_it->y, pred_it->x, pred_it->y);
      if (d < d_min){
        // New minimum distance, update observation association to current landmark
        obs_it->id = pred_it->id;
        d_min = d;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  
  // Get map landmarks and convert to vector of LandmarkObs
  std::vector<LandmarkObs> predicted;
  LandmarkObs landmark_obs;
  for (int ii = 0; ii < map_landmarks.landmark_list.size();ii++){
    Map::single_landmark_s landmark = map_landmarks.landmark_list[ii];
    landmark_obs.x = landmark.x_f;
    landmark_obs.y = landmark.y_f;
    landmark_obs.id = landmark.id_i;
    predicted.push_back(landmark_obs);
  }
  
  // Total Particle Weights
  double W = 0;
  for (std::vector<Particle>::iterator p_it = particles.begin(); p_it != particles.end();p_it++){
    // Particle data
    double p_x = p_it-> x;
    double p_y = p_it->y;
    double p_theta = p_it->theta;

    // Transform Observations to MAP
    std::vector<LandmarkObs> observations_map = observations;
    double o_x, o_y; // Observation X,Y in PARTICLE
    for (std::vector<LandmarkObs>::iterator p_it = observations_map.begin();
         p_it != observations_map.end(); p_it++){
      // Get Observation X,Y
      o_x = p_it->x; // P
      o_y = p_it->y; // P
      // Transform observations from PARTICLE to MAP (inplace)
      p_it->x = cos(p_theta)*o_x - sin(p_theta)*o_y + p_x; // MAP
      p_it->y = sin(p_theta)*o_x + cos(p_theta)*o_y + p_y; // MAP
      p_it->id = NULL;
    }
    
    // Set associations between prediction and observed
    dataAssociation(predicted, observations_map);
    std::vector<int> associations;
    std::vector<double> sense_x, sense_y;
    for (int ii = 0; ii < observations_map.size();ii++){
      associations.push_back(observations_map[ii].id);
      sense_x.push_back(observations_map[ii].x);
      sense_y.push_back(observations_map[ii].y);
    }
    SetAssociations(*p_it, associations, sense_x, sense_y);
    
    // Calculate Weights
    double p_weight = 1.0;
    for (int ii = 0; ii < observations_map.size();ii++){
      LandmarkObs obs = observations_map[ii];
      LandmarkObs lnd = predicted[obs.id - 1];
      double dx = obs.x - lnd.x;
      double dy = obs.y - lnd.y;
      double P_obs = 1.0/sqrt(2*M_PI*std_landmark[0]*std_landmark[1])*
        exp(-pow(dx,2)/(2*pow(std_landmark[0],2)))*
        exp(-pow(dy,2)/(2*pow(std_landmark[1],2)));
      p_weight *= P_obs;
    }
    p_it->weight = p_weight;
    // Update Total Weight
    W += p_weight;
  }
  
  // Normalize Weights
  for (int ii =0 ; ii < num_particles; ii++){
    particles[ii].weight = particles[ii].weight / W;
    weights[ii] = particles[ii].weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::default_random_engine generator;
  std::discrete_distribution<double> distribution(weights.begin(),weights.end());
  std::vector<Particle> new_particles;
  for (int ii = 0; ii < num_particles; ii++){
    int p_idx = distribution(generator);
    new_particles.push_back(particles[p_idx]);
  }
  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
