// MolecularDynamicsMod2.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include "NR/nr3.h"
#include "functiondefs.h"
using namespace std;

int main()
{
	// Initializing scalar variables
	const int N = 100; // number of particles
	const double tau = .001; // timestep
	const double T = .01; //temperature
	const int dim = 2; // number of dimensions
	const int G = 2*N-4; //degrees of freedom
	const int ITER = 20000; //number of iterations
	const double L = 11.3; // box size
	const int OutputFreq = 100; //Number of steps to take before outputting coordinates to text file
	const double Rthresh = 3; //radius threshold for calculating forces
	const double Rmin = 0; // minimum radius for interaction
	const double radmin = .9; //potential energy cutoff at which the energy is made constant (to eliminate explosive accelerations)
	double percentage = 0; //percentage done in calculation
	time_t starttime = clock();
	const double PI = std::atan(1.0) * 4;
	const double epsilon = 0; //radius correction
	double eta = 0.; //Nose-Hoover correction
	double KE = 0; // total kinetic energy
	double M_s = G*T; // Nose thermal inertia parameter

	//Initializing arrays
	double x[N][dim];  //particle positions
	double v[N][dim];  //particle velocities
	double g[N][dim] = { { 0. } };  //particle accelerations
	double g_old[N][dim] = { { 0. } }; //timestep k-1 acceleration
	double dist[N]; // interaction distance for a given particle (index of array is particle number)
	double deltax[N];  //interatomic distances for x component for a given particle
	double deltay[N];  // same for y component
	double Temp[ITER/OutputFreq]; //instantaneous temperature array
	int Tempindex = 0;
	bool energy = true; //whether or not to calculate energy
	double Etot[ITER / OutputFreq] = { 0. };  //Total energy of the system
	double Vtot[ITER / OutputFreq] = { 0. }; //Total potential energy
	double KEtot[ITER / OutputFreq] = { 0. }; //Total KE
	double time[ITER / OutputFreq] = { 0. };  //time elapsed in simulation
	int Eindex = 0;

	//pair correlation variables
	bool paircorr = true; // whether or not to calculate pair correlation
	const int pcn = 10000; // number of steps in pair correlation calculation
	double pair[pcn]= { 0. };  // pair correlation for pcn values of r before averaging over all particles
	double pairR[pcn] = { 0. };	//radii for the pair correlation function (from 0 to 2)
	double dr = .001;	//increment for pair correlation
	for (int rind = 1; rind < pcn+1; rind++){
		pairR[rind-1] = double(rind*dr);
	}

	//creating a 3D matrix in which to store the positions
	const int OutNum = floor(ITER / OutputFreq)+2;
	vector<vector<vector<double> > > xstore(N, vector<vector<double> >(dim, vector <double>(OutNum, 0.0)));
	int Outindex = 1; //index for timestep dimension of matrix

	//Intitializing positions in a simple square lattice
	int nxy = int(ceil(pow(N, 1.0 / dim))); // number of atoms in each direction
	double latconst = L / nxy;
	for (int i = 0; i < N; i++){
		x[i][0] = (i%nxy) * latconst;		// x-position
		x[i][1] = (floor(i / nxy)+1) * latconst; //y-position
		xstore[i][0][0] = x[i][0];	//storing initial positions
		xstore[i][1][0] = x[i][1];
	}

	//Initializing mean square displacement variables
	double MSD[N] = { 0. };
	double MSDavg[ITER / OutputFreq] = { 0. };
	double xinit[N][dim];
	int MSDindex = 0;

	//Initializing velocities using Maxwell distribution (m/kT =1)
	for (int i = 0; i < N; i++){
		for (int j = 0; j < dim; j++){
			v[i][j] = gasdev();
			xinit[i][j] = x[i][j]; //saving inital positions for mean square displacement
		}
	}

	//Initial velocity rescaling
	double vSqdSum = 0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < dim; j++){
			vSqdSum += v[i][j] * v[i][j];
		}
	}
	double lambda = sqrt(G* T / vSqdSum);
	for (int i = 0; i < N; i++){
		for (int j = 0; j < dim; j++){
			v[i][j] *= lambda;
		}
	}

	//Velocity Verlet
	for (int k = 0; k < ITER; k++){ //iterations in time step
		for (int i = 0; i < N; i++){ //iterations in particle number i
			for (int q = 0; q < N; q++){ // iterations in particle number q=/=i (for interaction distances)
				deltax[q] = x[i][0] - x[q][0] - (round((x[i][0] - x[q][0]) / L))*L;
				deltay[q] = x[i][1] - x[q][1] - (round((x[i][1] - x[q][1]) / L))*L;
				dist[q] = sqrt(pow(deltax[q], 2.0) + pow(deltay[q], 2.0)); //interaction distances between particle i and particles q
			}	

			//building list of neighbors within the threshold
			vector <double> NeighborRadii = threshdist(dist, Rthresh, Rmin, N);
			vector <int> NeighborIndicies = threshindex(dist, Rthresh, Rmin, N);
			int RadNum = NeighborRadii.size(), IndNum = NeighborIndicies.size();
			if (RadNum != IndNum){
				cout << "Error!  Sizes of NeighborRadii and NeighborIndicies are not equal.";
			}
			//computing Nose-Hoover correction and calculating the instantaneous kinetic energy
			double sum = 0;
			for (int j = 0; j < dim; j++){
					sum += v[i][j] * v[i][j];
			}
			KE = 0.5 * sum;
			eta = (1/M_s)*tau*(2*KE - G*T);
			//eta = 0;

			for (int j = 0; j < dim; j++){
				for (int vi = 0; vi < RadNum; vi++){
					if (NeighborRadii[vi] <= radmin){ //Setting an interaction radius minimum where the potential levels off
						NeighborRadii[vi] = radmin;
					}
					//initializing acceleration
					g_old[i][j] = g_old[i][j] + ((x[i][j] - x[NeighborIndicies[vi]][j])*(pow(NeighborRadii[vi]+epsilon, -14.0) - 0.5*pow(NeighborRadii[vi]+epsilon, -8.0))) - eta*v[i][j];
				}
			}
			//applying the velocity verlet algorithm
			for (int j = 0; j < dim; j++){
				x[i][j] = x[i][j] + tau*v[i][j] + 0.5 * pow(tau, 2.0) * g_old[i][j];  //particle position for k+1
				x[i][j] = x[i][j] - L*floor(x[i][j] / L);  //Brings particle back into the center box
				for (int vi = 0; vi < RadNum; vi++){
					g[i][j] = g[i][j] + (((x[i][j] - x[NeighborIndicies[vi]][j])*(pow(NeighborRadii[vi] + epsilon, -14.0) - 0.5*pow(NeighborRadii[vi] + epsilon, -8.0)))) - eta*v[i][j]; //acceleration components
				}
				v[i][j] = v[i][j] + 0.5* tau * (g[i][j] + g_old[i][j]);
				g[i][j] = g_old[i][j];
			}
			//calculating the pair correlation function for particle i on the last iteration
			if (k==ITER-1 && paircorr == true){
				for (int s = 0; s < pcn; s++){
					int m = threshNint(dist, (pairR[s] + dr), pairR[s], N);
					pair[s] += double(m) / (pairR[s] * dr);
				}
			}
			//Adding the energy of the present iteration to KE and V arrays
			if (k % OutputFreq == 0 && energy == true){
				for (int vi = 0; vi < RadNum; vi++){
					Vtot[Eindex] = Vtot[Eindex] + 4*((pow(NeighborRadii[vi] + epsilon, -12.0) - pow(NeighborRadii[vi] + epsilon, -6.0)));
					KEtot[Eindex] = KE;
					time[Eindex] = tau*k;
				}
			}
		}
		if (k % OutputFreq == 0){					
			//computing calculation progress
			percentage = (100.0*k / ITER);
			double timeelapsed = difftime(clock(), starttime);
			double tottime = timeelapsed*ITER / k;
			double timeremain = (tottime - timeelapsed)/1000;
			cout << "Percentage Completed: " << percentage << "\t Seconds remaining: " << timeremain << endl;
			
			//rescaling velocities
			double vSqdSum = 0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < dim; j++)
					vSqdSum += v[i][j] * v[i][j];
			double lambda = sqrt(G * T / vSqdSum);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < dim; j++)
					v[i][j] *= lambda;
			
			//storing position coordinates
			for (int i = 0; i < N; i++){
				for (int j = 0; j < dim; j++){
					xstore[i][j][Outindex] = x[i][j];
				}
			}
			//Calculating Mean Square Displacement
			for (int i = 0; i < N; i++){
				for (int j = 0; j < dim; j++){
					MSD[i] += pow(x[i][j] - xinit[i][j] - (round((x[i][j] - xinit[i][j]) / L))*L, 2.0);
				}
			}
			for (int i = 0; i < N; i++){
				MSDavg[MSDindex] += (1 / double(N))*MSD[i];
				//MSD[i] = 0.;
			}

			//calculating the instantaneous temperature
			Temp[Tempindex] = 2*KE / G;		
			
			//Updating indicies
			Tempindex += 1;
			Eindex += 1;
			Outindex += 1;
			MSDindex += 1;
 		}


	}

	//Calculating the Total Energy
	if (energy == true){
		for (int k = 0; k < ITER / OutputFreq; k++){
			Etot[k] = KEtot[k] + Vtot[k];
		}
	}
	//Printing Total Energy to a text file
	ofstream fout("total_energy.txt");
	std::stringstream Ebuffer;
	for (int s = 0; s < ITER / OutputFreq; s++){
		Ebuffer << time[s] << '\t' << Etot[s] << '\n';
	}
	fout << Ebuffer.str();
	fout.close();

	//Printing MSD to a text file
	fout.open("msd.txt");
	std::stringstream msdbuffer;
	for (int s = 0; s < ITER / OutputFreq; s++){
		msdbuffer << time[s] << '\t' << MSDavg[s] << '\n';
	}
	fout << msdbuffer.str();
	fout.close();

	//Printing position coordinates to a text file
	fout.open("molecular_dynamics.txt");
	std::stringstream mdbuffer;
	for (int k = 0; k < OutNum; k++) {
		for (int i = 0; i < N; i++){
			for (int j = 0; j < dim; j++){
				mdbuffer << xstore[i][j][k] << "\t";
			}
			mdbuffer << '\n';
		}
		mdbuffer << '\n';
	}
	fout << mdbuffer.str();
	fout.close();


	//Taking the average pair correlation over all particles
	for (int s = 0; s < pcn; s++){
		pair[s] = pair[s] * (1 / double(N));
	}
	//Printing Pair Correlation Function to a text file
	fout.open("paircorr.txt");
	std::stringstream pcbuffer;
	for (int s = 0; s < pcn; s++){
		pcbuffer << pairR[s] << '\t' << pair[s] << '\n';
	}
	fout << pcbuffer.str();
	fout.close();

	print1Ddoub(Temp, ITER / OutputFreq, "Inst. Temperature:   ");
	keep_window_open();
}