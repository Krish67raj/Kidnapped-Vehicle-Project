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
	
	if(!is_initialized){
	num_particles = 100;  // can put anything between 30 - 1000, as I have tested
	
	default_random_engine gen;
	weights.resize(num_particles);
	
	normal_distribution<double> noise_x(x, std[0]);
	normal_distribution<double> noise_y(y, std[1]);
	normal_distribution<double> noise_theta(theta, std[2]);
	
	for (int i = 0; i<num_particles ; i++){
		Particle  P;
		P.id = i;
		P.x = noise_x(gen);
		P.y = noise_y(gen);
		P.theta = noise_theta(gen);
		P.weight = 1.0;

		particles.push_back(P);
	}	
	is_initialized = true;	
	}	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	
	normal_distribution<double> noise_x(0, std_pos[0]);
	normal_distribution<double> noise_y(0, std_pos[1]);
	normal_distribution<double> noise_theta(0, std_pos[2]);
	
	for (int i = 0; i<num_particles ; i++){		 
		if (fabs(yaw_rate) > 0.0001){			// the value can be put between 0.001 to 0.00001 as I have tested
			particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + (yaw_rate*delta_t)) - sin(particles[i].theta)); 
			particles[i].y += (velocity/yaw_rate)*(-cos(particles[i].theta + (yaw_rate*delta_t)) + cos(particles[i].theta));
			particles[i].theta += yaw_rate*delta_t;
		}else{
			particles[i].x += velocity * cos(particles[i].theta) * delta_t;
			particles[i].y += velocity * sin(particles[i].theta) * delta_t;
		}		
		particles[i].x += noise_x(gen);
		particles[i].y += noise_y(gen);
		particles[i].theta += noise_theta(gen);		
	}	
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
		
	for (int i = 0; i < observations.size() ; i++){
		double lowest_dist = numeric_limits<double>::max();
			
		LandmarkObs Ob = observations[i];		
		
		for (int j = 0; j < predicted.size() ; j++){
			LandmarkObs Pr = predicted[j];
			double current_dist = dist(Ob.x,Ob.y,Pr.x,Pr.y);
			
			if(current_dist < lowest_dist){
				lowest_dist = current_dist;
				observations[i].id = j;
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
	
	double sig_x = std_landmark[0];
	double sig_y = std_landmark[1];
	
	double gauss_norm = 1/(2 * M_PI * sig_x * sig_y);
	double weight_norm = 0.0;
	
	for (int i = 0; i<num_particles; i++){
		
		// creating predicted measurement vector from map_landmarks
		vector<LandmarkObs> Pred;		
		for (int j = 0; j<map_landmarks.landmark_list.size() ; j++){			
			if (dist(map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f,particles[i].x,particles[i].y) < sensor_range){
				double pred_X = map_landmarks.landmark_list[j].x_f;
				double pred_Y = map_landmarks.landmark_list[j].y_f;
				int pred_id = map_landmarks.landmark_list[j].id_i;
				
				Pred.push_back(LandmarkObs{pred_id,pred_X,pred_Y});
			}	
		}				
		//Transforming the observations from car coordinate to map coordinates
		vector<LandmarkObs> New_obs;
		
		for (int j = 0; j<observations.size() ; j++){
			double o_x = observations[j].x;
			double o_y = observations[j].y;
			 
			double new_x = particles[i].x + (cos(particles[i].theta)*o_x) - (sin(particles[i].theta)*o_y);
			double new_y = particles[i].y + (sin(particles[i].theta)*o_x) + (cos(particles[i].theta)*o_y);
			
			New_obs.push_back(LandmarkObs{observations[j].id,new_x,new_y});			
		}
		// data association
		dataAssociation(Pred,New_obs);
		
		std::vector<int> associations;
		std::vector<double> sense_x;
		std::vector<double> sense_y;
		for(int j = 0; j < New_obs.size(); j++){
			sense_x.push_back(New_obs[j].x);
			sense_y.push_back(New_obs[j].y);
			associations.push_back(Pred[New_obs[j].id].id);
		}
		SetAssociations(particles[i],associations,sense_x,sense_y);
		
		// calculating weights 
		double prob = 1.0;
		for (int j = 0; j<observations.size() ; j++){
			
			int oid = New_obs[j].id;			
			double exponent = (((New_obs[j].x-Pred[oid].x)*(New_obs[j].x-Pred[oid].x))/(2*sig_x*sig_x)) + (((New_obs[j].y-Pred[oid].y)*(New_obs[j].y-Pred[oid].y))/(2*sig_y*sig_y));
			prob *= gauss_norm * exp(-exponent);	
		
			New_obs.clear();
			Pred.clear();				
		 }
		particles[i].weight = prob;
		weights[i] = particles[i].weight;		
		weight_norm += particles[i].weight;
	}	
	
	// normalizing weights
	for (int i = 0;i<particles.size(); i++){
		particles[i].weight /= weight_norm;
		weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	std:vector<Particle> resampled_particles;
	default_random_engine gen;
	double min_val = numeric_limits<double>::min();
	for(int i =0; i< num_particles; i++){
		if(particles[i].weight > min_val){
			min_val = particles[i].weight; 
		}
	}
	//creating uniform distributions for particles and weights
	uniform_int_distribution<int> dist_particles(0, num_particles-1);
	uniform_real_distribution<double> dist_weights(0, min_val);

	int index = dist_particles(gen);
	double beta = 0.0;
	for(int i = 0; i < num_particles; i++){
		beta += dist_weights(gen) * 2 * min_val;
		while(particles[index].weight<beta){
			beta -= particles[index].weight;
			index = (index+1)%num_particles;
		}
		resampled_particles.push_back(particles[index]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
	return particle;
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
