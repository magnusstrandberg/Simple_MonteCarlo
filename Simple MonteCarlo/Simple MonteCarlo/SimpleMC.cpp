#include <iostream>
#include <fstream>
#include <string>

using namespace std;
void fileReader(string);


int main() {
	cout << "Hellooou \n";
	string filelocation = "C:\\Users\\Magnus\\source\\repos\\Simple_MonteCarlo\\Simple MonteCarlo\\Simple MonteCarlo\\input.txt";
	fileReader (filelocation);

	std::cin.ignore(); //So we can see shit
	return 0;
}

void fileReader(string FileLocation) {
	string line;
	ifstream input(FileLocation);
	if (input.is_open())
	{
		while (getline(input, line))
		{
			
			if (!line.empty() && line.at(0) != '%')
			{
				cout << line << '\n';
			}
			
		}
		input.close();
	}

	else cout << "Unable to open file";

	return;
}