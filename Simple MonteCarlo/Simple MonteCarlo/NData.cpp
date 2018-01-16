#include "NData.h"
#include <iostream>
#include <fstream>
#include <string>


NData::NData()
{
}


NData::~NData()
{
}

void NData::readdatafile(std::string filename , int it)
{
	materialdata some_material;
	
	std::string line;
	std::ifstream input(filename);
	
	if (input.is_open())
	{
		getline(input, line);

		std::string::size_type split = line.find(' ');
		std::string::size_type old;

		some_material.Symbol = line.substr(0, split);

		old = split + 1;
		split = line.find(' ', old);

		some_material.Z = std::stoi(line.substr(old,split));

		old = split + 1;
		split = line.find(' ', old);

		some_material.A = std::stoi(line.substr(old, split));

		old = split + 1;
		split = line.find(' ', old);

		some_material.Aw = std::stod(line.substr(old, split));

		old = split + 1;
		split = line.find(' ', old);

		some_material.T = std::stod(line.substr(old, split));

		getline(input, line);

		some_material.Fissionable = std::stoi(line);

		if (some_material.Fissionable)
		{
			for (int i = 0; i < some_material.Fissionable; i++)
			{
				getline(input, line);
				std::string key, value;

				std::string::size_type split = line.find(' ');

				key = line.substr(0, split);
				value = line.substr(split + 1);
				some_material.Fission_nubar[std::stod(key)] = std::stod(value);
			}
		}

		while (getline(input, line))
		{
			MT a_mt;
			split = line.find(' ');
			a_mt.MT_number = std::stoi(line.substr(0, split));

			old = split + 1;
			split = line.find(' ', old);

			a_mt.Q_value = std::stod(line.substr(old, split));

			old = split + 1;
			split = line.find(' ', old);

			int amount = std::stoi(line.substr(old, split));

			for (int i = 0; i < amount; i++)
			{
				getline(input, line);
				std::string key, value;

				std::string::size_type split = line.find(' ');

				key = line.substr(0, split);
				value = line.substr(split + 1);
				a_mt.Ec_pair[std::stod(key)] = std::stod(value);
			}
			some_material.MT_data.push_back(a_mt);

		}

		materials[it] = some_material;

		input.close();
	}

	else std::cout << "Unable to open file" + it;

	return;
}

void NData::readlist(std::string filename)
{
	std::string line;
	std::ifstream input(filename);

	std::vector <std::string> data;

	while (getline(input, line))
	{
		data.push_back(line);
	}

	materials.resize(data.size());

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < data.size(); i++)
	{
		readdatafile(data[i], i);
	}

	return;
}


