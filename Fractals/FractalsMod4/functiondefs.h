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

//creates an array containing the coordinates of all the diagonal neighbors of (x,y)
vector<int> neighbors(int xin, int yin, int N, vector<int> &arr){
	//  xlup = xin - 1, ylup = yin - 1,
	//	xrup = xin + 1, yrup = yin - 1,
	//	xldown = xin - 1, yldown = yin + 1,
	//	xrdown = xin + 1, yrdown = yin + 1; //neighboring diagonal blocks
	//  xup = xin, yup = yin - 1,
	//	xl = xin - 1, yl = yin,
	//	xdown = xin, ydown = yin + 1,
	//	xr = xin + 1, yr = yin;		// neighboring edge blocks
	arr = { xin - 1, yin - 1, xin + 1, yin - 1, xin - 1, yin + 1, xin + 1, yin + 1, xin, yin - 1, xin - 1, yin, xin, yin + 1, xin + 1, yin }; // arrays of xy coordinates of neighbors
	for (int i = 0; i < 16; i ++){
		arr[i] -= (N)*int(floor(float(arr[i]) / N));
	}
	return arr;
}

//returns true if diagonal neighbor of (x,y) is occupied by a particle
bool neighboristrue(int xin, int yin, int N, vector<vector<bool>> &mat, vector<int> &arr){
	vector<int> vec;
	vec = neighbors(xin, yin, N, arr);
	bool out = false;
	for (int i = 0; i < 16; i += 2){
		out = (out || mat[vec[i]][vec[i + 1]]);
	}
	return out;
}