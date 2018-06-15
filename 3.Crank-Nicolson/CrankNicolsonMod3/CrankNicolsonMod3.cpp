// CrankNicolsonMod3.cpp : Defines the entry point for the console application.
//
//
#include "stdafx.h"
#include <iostream>
#include "NR/nr3.h"
#include <fstream>
#include <complex>
#include "NR/sort.h"
#include "NR/sparse.h"
#include "NR/tridag.h"
#include "functiondefs.h"

int main()
{
	//Tunable Parameters
	const int N = 10000; //number of domain elements
	double L = 20; //box size
	double dt = .0001; //timestep
	double k0 = 1; //wavenumber
	double x0 = L / 2; //inital position
	double sig0 = .1; //gaussian width
	const int ts = 500; //total number of timesteps
	double Vloc = L / 2; //location of potential barrier
	double Vo = 0.;
	//double Vo = 10000; // potential height
	double Vw = .5; //potential half width
	double V[N] = {0.}; //potential barrier vector

	//Other Stuff
	double pi = 3.14159265359;
	Complex li = { 0.0, 1.0 };
	Complex alpha = 1.*li;
	VecDoub x(N+2);
	VecComplex psi(N+2);  // wavepacket psi
	VecComplex psic(N+2); // complex conjugate of psi
	VecComplex psinew(N); //psi for time n+1
	VecComplex b(N);  // b vector (Ax=Bc=b where we are solving for x)
	VecDoub rho(N+2); //(psi*)(psi)
	int t[ts] = { 0. };
	for (int i = 0; i < ts; i++){
		t[i] = i;
	}
	double dx = L / double(N);
	Complex mu = 2*pow(dx, 2.0) / (dt*alpha); // crank nicolson mu
	double area[ts] = { 0. }; //area under the curve	
	
	//Implicit scheme variables
	VecComplex psiimp(N + 2);  // wavepacket psi
	VecComplex psicimp(N + 2); // complex conjugate of psi
	VecComplex psinewimp(N); //psi for time n+1
	Complex muimp = pow(dx, 2.0) / (dt*alpha); // implicit scheme mu
	VecDoub rhoimp(N + 2);
	double areaimp[ts] = { 0. }; //area under the curve	
	
	//Explicit scheme variables
	VecComplex psiexp(N + 2);
	VecComplex psicexp(N + 2);
	VecComplex psinewexp(N);
	VecDoub rhoexp(N + 2);
	double areaexp[ts] = { 0. };

	//setting boundary values to 0
	x[0] = 0.; x[N+1] = 0.; 
	psi[0] = { 0., 0. }; psi[N + 1] = { 0., 0. };
	psic[0] = { 0., 0. }; psic[N + 1] = { 0., 0. };
	rho[0] = 0.; rho[N + 1] = 0.;
	psiimp[0] = { 0., 0. }; psiimp[N + 1] = { 0., 0. };
	psicimp[0] = { 0., 0. }; psicimp[N + 1] = { 0., 0. };
	rhoimp[0] = 0.; rhoimp[N + 1] = 0.;
	psiexp[0] = { 0., 0. }; psiexp[N + 1] = { 0., 0. };
	psicexp[0] = { 0., 0. }; psicexp[N + 1] = { 0., 0. };
	rhoexp[0] = 0.; rhoexp[N + 1] = 0.;


	//defining the domain
	double q = 0.;
	for (int i = 1; i < N + 1; i++){
		x[i] = q;
		q += dx;
	}

	//Defining a square potential
	for (int i = 1; i < N + 1; i++){
		if (x[i] >= Vloc - Vw && x[i] <= Vloc + Vw){
			V[i - 1] = Vo;
		}
	}
//	print1Ddoub(V, N, "V:  ");

	//Tridiagonal matrix vectors for Crank Nicolson
	VecComplex ldiag(N), rdiag(N), middiag(N);
	for (int i = 0; i < N; i++){
		ldiag[i] = Complex(-1);
		rdiag[i] = Complex(-1);
	}
	ldiag[0] = { 0., 0. }; rdiag[N] = { 0., 0. };
	for (int i = 0; i < N; i++){
		middiag[i] = Complex(2) + mu - .5*alpha*mu*dt*V[i];
	}

	//Tridiagonal vectors for implicit scheme
	VecComplex ldiagimp(N), rdiagimp(N), middiagimp(N);
	for (int i = 0; i < N; i++){
		ldiagimp[i] = ldiag[i] / (muimp);
		rdiagimp[i] = ldiagimp[i];
	}
	ldiagimp[0] = { 0., 0. }; rdiagimp[N] = { 0., 0. };
	for (int i = 0; i < N; i++){
		middiagimp[i] = (Complex(2) + muimp) / muimp;
	}

	//initializing wavepacket for Crank Nicolson
	Complex expon;
	for (int i = 1; i < N+1; i++){
		expon = -(pow((x[i] - x0), 2.0) / (2 * pow(sig0, 2.0))) + (li*k0*x[i]);
		psi[i] = pow(pi*pow(sig0, 2.), -.25)*exp(expon);
		psic[i] = conj(psi[i]);
		rho[i] = real(psi[i] * psic[i]);
	}

	//initializing wavepacket for implicit and explicit schemes
	for (int i = 1; i < N + 1; i++){
		psiimp[i] = psi[i];
		psicimp[i] = psic[i];
		rhoimp[i] = rho[i];
		psiexp[i] = psi[i];
		psicexp[i] = psic[i];
		rhoexp[i] = rho[i];
	}


	//print1DVdoub(rho, N+2, "rho: ");
	
	ofstream fout("psi.txt");
	std::stringstream buffer;

	for (int ti = 0; ti < ts; ti++){
		for (int i = 0; i < N; i++){
			b[i] = psi[i] + psi[i + 1] * (mu - Complex(2) + .5*alpha*mu*dt*V[i]) + psi[i + 2];
		}
		psinew = tridagComp(ldiag, middiag, rdiag, b, psinew);
		for (int i = 1; i < N+1; i++){
			psi[i] = psinew[i-1];
			psic[i] = conj(psi[i]);
			rho[i] = real(psi[i] * psic[i]);
			buffer << x[i] << "\t" << rho[i] << "\n";
			area[ti] += dx*rho[i];
		}
		buffer << endl;
	}
	fout << buffer.str();
	fout.close();

	fout.open("psiimp.txt");
	std::stringstream bufferimp;
	for (int ti = 0; ti < ts; ti++){
		psinewimp = tridagComp(ldiagimp, middiagimp, rdiagimp, psiimp, psinewimp);
		for (int i = 1; i < N + 1; i++){
			psiimp[i] = psinewimp[i - 1];
			psicimp[i] = conj(psiimp[i]);
			rhoimp[i] = real(psiimp[i] * psicimp[i]);
			bufferimp << x[i] << "\t" << rhoimp[i] << "\n";
			areaimp[ti] += dx*rhoimp[i];
		}
		bufferimp << endl;
	}
	fout << bufferimp.str();
	fout.close();
	
	fout.open("psiexp.txt");
	std::stringstream bufferexp;
	for (int ti = 0; ti < ts; ti++){
		for (int i = 1; i < N + 1; i++){
			bufferexp << x[i] << "\t" << rhoexp[i] << "\n";
			psiexp[i] = psiexp[i] + (li * dt / pow(dx, 2)) * (psiexp[i + 1] - Complex(2) * psiexp[i] + psiexp[i - 1]);
			psicexp[i] = conj(psiexp[i]);
			rhoexp[i] = real(psiexp[i] * psicexp[i]);
//			double sigt = sig0*sqrt(1 + ti*ti*dt*dt / (pow(sig0, 4.0)));
//			rhoexp[i] = pow(pi*pow(sigt, 2.0), -.5)*exp(-1*pow(x[i] - x0 - k0*ti*dt, 2.0) / pow(sigt, 2.0));
			areaexp[ti] += dx*rhoexp[i];
		}
		bufferexp << endl;
	}
	fout << bufferexp.str();
	fout.close();



	fout.open("Vout.txt");
	std::stringstream bufferV;
	for (int i = 0; i < N; i++){
		bufferV << x[i + 1] << "\t" << V[i] << "\n";
	}
	fout << bufferV.str();
	fout.close();
	
	fout.open("Area.txt");
	std::stringstream bufferArea;
	for (int i = 0; i < ts; i++){
		bufferArea << t[i] << "\t" << area[i] << "\t" << areaimp[i] << "\t" << areaexp[i] << "\n";
	}
	fout << bufferArea.str();
	fout.close();

	keep_window_open();
}

