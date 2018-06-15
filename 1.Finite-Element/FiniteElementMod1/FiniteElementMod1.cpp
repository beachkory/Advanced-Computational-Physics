// FiniteElementMod1.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>
#include "nr3.h"
#include <fstream>
#include "sort.h"
#include "sparse.h"
#include "tridag.h"

using namespace std;

inline void keep_window_open()
{
	cin.clear();
	cout << "Please enter a character to exit\n";
	char ch;
	cin >> ch;
	return;
}
string printarray(VecDoub &arg, int length) {
	std::stringstream buffer;
	for (int n = 0; n<length; ++n)
		buffer << arg[n] << ',';
	cout << buffer.str() << '\n';
	return buffer.str();
}

int main()
{
	int n = 150;  // number of elements in arrays x and y
	Doub xmin = 0, xmax = 1;  //defining the domain
	Doub h = (xmax - xmin) / (n - 1); // grid spacing
	VecDoub y(n-2);  // defining the array y(x)
	VecDoub x(n);  // defining x-position array
	VecDoub b(n-2); // RHS of Ax=b (i.e. the source integrals)
	VecDoub middiag(n-2); // vector containing elements of the center diagonal in A
	VecDoub ldiag(n-2); // vector containing left bottom diagonal in A
	VecDoub rdiag(n-2); // vector containing right top diagonal in A
	x[0] = xmin;
	for (int i = 1; i < n; i++) {
		x[i] = x[i - 1] + h;
	}
	for (int i = 0; i < n-2; i++) {
		middiag[i] = 2 / h;
	}
	for (int i = 1; i < n-1; i++) {
		b[i-1]=-1*pow(h, 3.0) * (6 * pow(i, 2.0) + 1);
	}
	for (int i = 0; i < n - 3; i++) {
		ldiag[i+1] = -1 / h;
		rdiag[i] = -1 / h;
	}

	ofstream fout("finite_element.gp");
	fout << printarray(x, n) << endl;
	fout << printarray(tridag(ldiag, middiag, rdiag, b, y), n - 2) << endl;
	fout.close();
	keep_window_open();
}