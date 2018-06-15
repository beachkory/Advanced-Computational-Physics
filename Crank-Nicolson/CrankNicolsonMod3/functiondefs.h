inline void keep_window_open()
{
	cin.clear();
	cout << "Please enter a character to exit\n";
	char ch;
	cin >> ch;
	return;
}
void printarray(vector<vector<bool>> arg, int cols, int rows) {
	std::stringstream buffer;
	for (int i = 0; i < arg.size(); i++)
	{
		for (int j = 0; j < arg[i].size(); j++)
		{
			buffer << arg[i][j] << ' ';
		}
		buffer << '\n';
	}
	cout << buffer.str();
	//return buffer.str();
}
void print1Dint(vector<int> arg, int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
	//return buffer.str();
}
void print1DCdoub(VecComplex arg, int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
}
void print1DVdoub(VecDoub arg, int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
}
void print1Ddoub(double *arg, int length, string name) {
	std::stringstream buffer;
	for (int n = 0; n < length; ++n){
		buffer << arg[n] << ',';
	}
	cout << name << buffer.str() << '\n';
}
VecComplex_O tridagComp(VecComplex_I &a, VecComplex_I &b, VecComplex_I &c, VecComplex_I &r, VecComplex_O &u)
{
	Int j, n = a.size();
	Complex bet;
	VecComplex gam(n);
	if (b[0] == Complex(0)) throw("Error 1 in tridag");
	u[0] = r[0] / (bet = b[0]);
	for (j = 1; j<n; j++) {
		gam[j] = c[j - 1] / bet;
		bet = b[j] - a[j] * gam[j];
		if (bet == Complex(0)) throw("Error 2 in tridag");
		u[j] = (r[j] - a[j] * u[j - 1]) / bet;
	}
	for (j = (n - 2); j >= 0; j--)
		u[j] -= gam[j + 1] * u[j + 1];
	return u;
}