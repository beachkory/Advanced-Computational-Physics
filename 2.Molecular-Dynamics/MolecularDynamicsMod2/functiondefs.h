inline void keep_window_open()
{
	cin.clear();
	cout << "Please enter a character to exit\n";
	char ch;
	cin >> ch;
	return;
}
void print1Dint(int arg[], int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
	//return buffer.str();
}
void print1Ddoub(double arg[], int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
	//return buffer.str();
}
void printarray(double * arg, int cols, int rows, string name) {
	std::stringstream buffer;
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			buffer << *(arg + j*rows + i) << ' ';
		}
		buffer << '\n';
	}
	cout << name << buffer.str();
	//return buffer.str();
}

//Initializes velocities with a Maxwell distribution
double gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if (!available) {
		do {
			v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq) / rsq);
		gset = v1 * fac;
		available = true;
		return v2*fac;
	}
	else {
		available = false;
		return gset;
	}
}

//returns the number of particles within a given radius threshold
const int threshNint(double arg[], double rthresh, double rmin, int N){ //outputs minimum value of array (not including 0)
	int threshN = 0;
	for (int i = 0; i < N; i++){
		if ((arg[i] <= rthresh) && (arg[i] != 0) && (arg[i] >= rmin)){
			threshN += 1;
		}
	}
	const int threshNconst = threshN;
	return threshNconst;
}
//returns an array of the distances of the neighbors within the threshold
vector<double> threshdist(double arg[], double rthresh, double rmin, int arraysize){
	const int Nthresh = threshNint(arg, rthresh, rmin, arraysize);
	vector<double> threshvector(Nthresh);
	int index = 0;
	for (int i = 0; i < arraysize; i++){
		if ((arg[i] <= rthresh) && (arg[i] != 0) && (arg[i] >= rmin)){
			threshvector[index] = arg[i];
			index += 1;
		}
	}
	return threshvector;
}

//returns an array that contains the indicies of the neighbors within the threshod
vector<int> threshindex(double arg[], double rthresh, double rmin, int arraysize){
	const int Nthresh = threshNint(arg, rthresh, rmin, arraysize);
	vector<int> threshvector(Nthresh);
	int index = 0;
	for (int i = 0; i < arraysize; i++){
		if ((arg[i] <= rthresh) && (arg[i] != 0) && (arg[i] >= rmin)){
			threshvector[index] = i;
			index += 1;
		}
	}
	return threshvector;
}

//outputs minimum value of array (not including 0)
int minindex(double arg[], int N){
	int index = 0;
	if (arg[index] == 0){
		index += 1;
	}
	double minimum = arg[index];
	for (int i = 0; i < N; i++){
		if ((arg[i] < minimum) && (arg[i] != 0)){
			minimum = arg[i];
			index = i;
		}
	}
	return index;
}

//outputs minimum value of array (not including 0)
double minnumber(double arg[], int N){
	int index = 0;
	if (arg[index] == 0){
		index += 1;
	}
	double minimum = arg[index];
	for (int i = 0; i < N; i++){
		if ((arg[i] < minimum) && (arg[i] != 0)){
			minimum = arg[i];
			index = i;
		}
	}
	return minimum;
}