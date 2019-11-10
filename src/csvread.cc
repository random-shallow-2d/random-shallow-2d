#include <complex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int main() {
  ifstream myfile;
  myfile.open("foo.csv", ios::in);
  string line;
  istringstream is;
  complex<double> c;
  if (myfile.is_open()) {
    for (int l = 1; l <= 8; l++) {
      getline(myfile, line, ';');
      istringstream is(line);
      is >> c;
      cout << c << endl;
    }
  }

  return 0;
}
