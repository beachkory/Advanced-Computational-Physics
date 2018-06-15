// FractalsMod4.cpp : Defines the entry point for the console application.
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

int main()
{
	const int N = 512; // box size
	const int ITER = 2000;  // number of iterations
	time_t starttime = clock();
	vector <vector<bool>> mat;
	mat.resize(N, vector<bool>(N, false)); //N x N matrix containing coordinates of particles
	int nucl = 5; // number of nucleation perimeters
	int nucsites =int(pow(2*nucl+1, 2)); // initalizing variable for number of nucleation sites
	vector<int> neigharr(16); //preallocates vector for neighboring coordinates
	cout << nucsites << '\n';

	for (int i = -1 * nucl; i < nucl + 1; i++){
		for (int j = -1 * nucl; j < nucl + 1; j++){
			mat[N / 2 + i][N / 2 + j] = true;// initialization of first particle position
		}
	}


//	printarray(mat, N, N);
	int x = 0; //x-coordinate of moving particle
	int y = 0; //y-coordinate of moving particle
	int b = 5; //border parameter
	int bfac = 1; //loop-varying border parameter factor

	//prepopulating a list of random numbers
	int q = 0, r = 0; //random number generation variable
	srand(time(NULL)); //random number seed
	int enterrandx[10000];
	int enterrandy[10000];
	int qrands[10000];
	int rrands[10000];
	for (int i = 0; i < 10000; i++){
		q = 0, r = 0;
		while (q == 0){
			q = (rand() % 3) - 1;
		}
		qrands[i] = q;
		while (r == 0){
			r = (rand() % 3) - 1;
		}
		rrands[i] = r;
		enterrandx[i] = rand() % N;
		enterrandy[i] = rand() % N;
	}
	int randindex = 0;
	int xentrandindex = 0;
	int yentrandindex = 1;

	//Introduces a new particle into the system for every iteration
	for (int i = 0; i < ITER; i++){
		x = enterrandx[(xentrandindex % 10000)] * (xentrandindex % 2); //Initializing new particle
		y = enterrandy[(yentrandindex % 10001)] * (yentrandindex % 2); //Initializing new particle
		xentrandindex += 1;
		yentrandindex += 1;
		cout << x << ' ' << y << endl;
		cout << i << endl;
		while (!(neighboristrue(x, y, N, mat, neigharr))){

			mat[x][y] = false;

			//bfac = 1;
			//if (!(x >= N / 2 - N / 10 && x <= N / 2 + N / 10 && y >= N / 2 - N / 10 && y <= N / 2 + N / 10) && (x >= b + 1 && y >= b + 1 && x <= N - b + 1 && y <= N - b + 1)){
			//	bfac *= b;
			//}
			x += qrands[bfac*randindex%10000];
			x -= (N)*int(floor(float(x) / N));
			y += rrands[bfac*randindex%10000];
			y -= (N)*int(floor(float(y) / N));
			randindex += 1;

			mat[x][y] = true;
		}
		//cout << difftime(clock(), starttime);
		//printarray(mat, N, N);
}
	vector<int> xout(ITER + nucsites);
	vector<int> yout(ITER + nucsites);
	int index = 0;

	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			if (mat[i][j] == true){
				xout[index] = i;
				yout[index] = j;
				index += 1;
			}
		}
	}
	ofstream fout("fractal.txt");
	std::stringstream buffer;
	for (int i = 0; i < ITER + nucsites; i++){
		buffer << xout[i] << '\t' << yout[i] << '\n';
	}
	fout << buffer.str();
	fout.close();



//	printarray(mat, N, N);
	keep_window_open();
}

