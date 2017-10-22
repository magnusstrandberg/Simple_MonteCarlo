#include "Input.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

using namespace std;


Input::Input()
{
}


Input::~Input()
{
}


//legacy
void Input::defaultValues() {
	N = -999;
	M = -999;
	Needle = -999;
	Coin = -999;
	l = -999;
	L = -999;
}


void Input::fileReader(const std::string filelocation)
{
	Input::defaultValues();
	string line;
	ifstream input(filelocation);
	if (input.is_open())
	{
		while (getline(input, line))
		{

			if (!line.empty() && line.at(0) != '%')
			{
				//cout << line << "\n";
				dataParser(line);
			}

		}
		input.close();
	}

	else cout << "Unable to open file";

	return;
}

void Input::dataParser(const std::string line)
{
	bool parameter = true;
	string para, value;
	
	string::size_type split = line.find('=');

	para = line.substr(0, split);
	para.erase(std::remove(para.begin(), para.end(), ' '), para.end());

	value = line.substr(split + 1);
	value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
	


	if (!para.empty() && !value.empty()) { createDataSet(para, value); }
	return;
}

void Input::createDataSet(const std::string para, std::string value)
{
	
	if (para == "Surf")
	{
		surfCreator(value);
		return;
	}
	if (para == "Cell")
	{

		return;
	}
	if (para == "Subspace")
	{

		return;
	}
	//Simple numerical values can be handled directly, removing ;
	value.erase(std::remove(value.begin(), value.end(), ';'), value.end());
	if (para == "N") { 
		N = stoi(value);
		return;
	}
	
	if (para == "M") { 
		M = stoi(value);
		return;
	}
	if (para == "Needle") {
		Needle = stoi(value);
		return;
	}
	if (para == "Coin") {
		Coin = stoi(value);
		return;
	}
	if (para == "L") {
		L = stod(value);
		return;
	}
	if (para == "l") {
		l = stod(value);
		return;
	}

	return;

}

/*
A=0	B=1
C=2	D=3
E=4 F=5
G=6 H=7
I=8 J=9
*/

void Input::surfCreator(const std::string data)
{
	
	string type, values;

	string::size_type split = data.find(';');
	string::size_type old;

	type = data.substr(0, split);
	values = data.substr(split + 1);
	
	Surf_input some_surface;

	if (type == "para_plane")
	{
		surfParPlane(values, some_surface);
		return;
		
	}
	if (type == "plane")
	{
		surfPlane(values, some_surface);
		return;
	}
	if (type == "cylinder")
	{
		surfCylinder(values, some_surface);
		return;
	}

	if (type == "sphere")
	{
		surfSphere(values, some_surface);
		return;
	}
	if (type == "cubid")
	{
		surfCubid(values, some_surface);
		return;
	}

	if (type == "square_prism_inf")
	{
		surfSquarePrismInf(values, some_surface);
		return;
	}

	return;
	
	
}

void Input::surfSphere(const std::string values, Surf_input sphe)
{
	//type;ComplexID;Rank;radius;centre point;
	double centre_point[3], radius;
	sphe.complement.push_back(0);
	sphe.surf_ids.push_back(1);
	sphe.type = sphere;

	string::size_type split;
	string::size_type old;
	split = values.find(';');
	//comp ID
	sphe.complex_id = stoi(values.substr(0, split));

	old = split + 1;
	split = values.find(';', old);
	//Subspace rank
	sphe.subspace_rank = stoi(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//radius
	radius = stod(values.substr(old, split));

	string::size_type old_point, point;

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		centre_point[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	std::vector <double> components(10, 0);

	components[0] = 1;
	components[1] = 1;
	components[2] = 1;

	components[6] = (-2 * centre_point[0]);
	components[7] = (-2 * centre_point[1]);
	components[8] = (-2 * centre_point[2]);

	components[9] = (pow(centre_point[0], 2) + pow(centre_point[1], 2)
		+ pow(centre_point[2], 2) - pow(radius, 2));

	sphe.Surfparams.push_back(components);


	Complex_surf_input.push_back(sphe);
	return;
}

void Input::surfCylinder(const std::string values, Surf_input cyl)
{
	//type;ComplexID;Rank;height(0 is infinite);parallel to; centre point (at bottom if height); radius;
	double height, centre_point[3], radius;
	std::vector <bool> comp_cyl;
	int parallel_to;

	string::size_type split;
	string::size_type old;

	split = values.find(';');
	//comp ID
	cyl.complex_id = stoi(values.substr(0, split));

	old = split + 1;
	split = values.find(';', old);
	//Subspace rank
	cyl.subspace_rank = stoi(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//height
	height = stod(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//parallel to
	parallel_to = stoi(values.substr(old, split));

	string::size_type old_point, point;

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		centre_point[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	old = split + 1;
	split = values.find(';', old);

	radius = stod(values.substr(old, split));

	comp_cyl.push_back(0);
	std::vector <double> cylcomponents(10, 0);


	std::vector <double> botcomponents(10, 0);
	std::vector <double> topcomponents(10, 0);


	if (parallel_to == 1)
	{
		cylcomponents[1] = 1;
		cylcomponents[2] = 1;

		cylcomponents[7] = (-2 * centre_point[1]);
		cylcomponents[8] = (-2 * centre_point[2]);

		cylcomponents[9] = (pow(centre_point[1], 2) + pow(centre_point[2], 2) - pow(radius, 2));

		if (height != 0)
		{
			botcomponents[6] = 1;
			botcomponents[9] = (-1 * centre_point[0]);
			comp_cyl.push_back(1);
			topcomponents[6] = 1;
			topcomponents[9] = (-1 * (centre_point[0] + height));
			comp_cyl.push_back(0);
		}
	}
	else if (parallel_to == 2)
	{
		cylcomponents[0] = 1;
		cylcomponents[2] = 1;

		cylcomponents[6] = (-2 * centre_point[0]);
		cylcomponents[8] = (-2 * centre_point[2]);

		cylcomponents[9] = (pow(centre_point[0], 2) + pow(centre_point[2], 2) - pow(radius, 2));

		if (height != 0)
		{
			botcomponents[7] = 1;
			botcomponents[9] = (-1 * centre_point[1]);
			comp_cyl.push_back(1);
			topcomponents[7] = 1;
			topcomponents[9] = (-1 * (centre_point[1] + height));
			comp_cyl.push_back(0);
		}

	}
	else
	{
		cylcomponents[0] = 1;
		cylcomponents[1] = 1;

		cylcomponents[6] = (-2 * centre_point[0]);
		cylcomponents[7] = (-2 * centre_point[1]);

		cylcomponents[9] = (pow(centre_point[0], 2) + pow(centre_point[1], 2) - pow(radius, 2));

		if (height != 0)
		{
			botcomponents[8] = 1;
			botcomponents[9] = (-1 * centre_point[2]);
			comp_cyl.push_back(1);
			topcomponents[8] = 1;
			topcomponents[9] = (-1 * (centre_point[2] + height));
			comp_cyl.push_back(0);
		}
	}



	cyl.Surfparams.push_back(cylcomponents);
	cyl.surf_ids.push_back(1);
	cyl.type = cylinder_inf;
	cyl.complement = comp_cyl;

	if (height != 0) {
		cyl.Surfparams.push_back(botcomponents);
		cyl.Surfparams.push_back(topcomponents);
		cyl.surf_ids.push_back(2);
		cyl.surf_ids.push_back(3);
		cyl.type = cylinder;
	}

	Complex_surf_input.push_back(cyl);
	return;

}

void Input::surfParPlane(const std::string values, Surf_input para_plane)
{
	//type;ComplexID;Rank;xyz;cross point
	//standard info
	para_plane.type = plane;
	para_plane.complement.push_back(0);
	//Find complex ID
	string::size_type split;
	string::size_type old;
	
	split = values.find(';');

	para_plane.complex_id = stoi(values.substr(0, split));
	//Find Rank
	old = split + 1;
	split = values.find(';', old);

	para_plane.subspace_rank = stoi(values.substr(old, split));

	//Components
	old = split + 1;
	split = values.find(';', old);

	std::vector <double> components(10, 0);

	if (stoi(values.substr(old, split)) == 1)
	{
		components[6] = 1;
	}
	else if (stoi(values.substr(old, split)) == 2)
	{
		components[7] = 1;
	}
	else if (stoi(values.substr(old, split)) == 3)
	{
		components[8] = 1;
	}

	old = split + 1;
	split = values.find(';', old);

	components[9] = (-1 * stod(values.substr(old, split)));
	//adding compnents to parameter list
	para_plane.Surfparams.push_back(components);
	para_plane.surf_ids.push_back(1);

	//adding all to surf input vector

	Complex_surf_input.push_back(para_plane);

	return;
}

void Input::surfPlane(const std::string values, Surf_input threeplane)
{
	//type;ComplexID;Rank;point1;point2;point3;
	//standard info
	threeplane.type = plane;
	threeplane.complement.push_back(0);

	string::size_type split;
	string::size_type old;

	//Find complex ID
	split = values.find(';');

	threeplane.complex_id = stoi(values.substr(0, split));
	//Find Rank
	old = split + 1;
	split = values.find(';', old);

	threeplane.subspace_rank = stoi(values.substr(old, split));

	double point1[3];
	double point2[3];
	double point3[3];
	std::vector <double> components;

	string::size_type old_point, point;

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		point1[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		point2[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		point3[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	components = generalPlanParams(point1,point2,point3);
	threeplane.Surfparams.push_back(components);
	threeplane.surf_ids.push_back(1);

	Complex_surf_input.push_back(threeplane);

	return;
}

std::vector<double> Input::generalPlanParams(double point1[], double point2[], double point3[])
{
	std::vector <double> components(10, 0);

	double G, H, I, J;

	G = ((point2[1] * point3[2])
		- (point3[1] * point2[2])
		- (point1[1] * (point3[2] - point2[2]))
		+ (point1[2] * (point3[1] - point2[1]))
		);

	H = ((point2[2] * point3[0])
		- (point3[2] * point2[0])
		- (point1[2] * (point3[0] - point2[0]))
		+ (point1[0] * (point3[2] - point2[2]))
		);

	I = ((point2[0] * point3[1])
		- (point3[0] * point2[1])
		- (point1[0] * (point3[1] - point2[1]))
		+ (point1[1] * (point3[0] - point2[0]))
		);

	J = (((-1 * point1[0]) * ((point2[1] * point3[2]) - (point3[1] * point2[2])))
		+ (point1[1] * ((point2[0] * point3[2]) - (point3[0] * point2[2])))
		- (point1[2] * ((point2[0] * point3[1]) - (point3[0] * point2[1]))))
		;

	components[6] = G;
	components[7] = H;
	components[8] = I;
	components[9] = J;

	return components;
}

void Input::surfCubid(const std::string values, Surf_input cuboid )
{
	//type;ComplexID;Rank;xlength;ylength;zlength;centre point;
	double start_point[3], xlength, ylength, zlength;
	cuboid.type = cubid;

	string::size_type split;
	string::size_type old;

	split = values.find(';');
	//comp ID
	cuboid.complex_id = stoi(values.substr(0, split));

	old = split + 1;
	split = values.find(';', old);
	//Subspace rank
	cuboid.subspace_rank = stoi(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//xlength
	xlength = stod(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//ylength
	ylength = stod(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//zlength
	zlength = stod(values.substr(old, split));

	string::size_type old_point, point;

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		start_point[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	std::vector <double> xcomponents(10, 0);
	std::vector <double> x2components(10, 0);
	std::vector <double> ycomponents(10, 0);
	std::vector <double> y2components(10, 0);
	std::vector <double> zcomponents(10, 0);
	std::vector <double> z2components(10, 0);

	xcomponents[6] = 1;
	xcomponents[9] = (-1 * start_point[0]);
	cuboid.Surfparams.push_back(xcomponents);
	cuboid.complement.push_back(1);
	cuboid.surf_ids.push_back(1);
	x2components[6] = 1;
	x2components[9] = (-1 * (start_point[0] + xlength));
	cuboid.Surfparams.push_back(x2components);
	cuboid.complement.push_back(0);
	cuboid.surf_ids.push_back(2);

	ycomponents[7] = 1;
	ycomponents[9] = (-1 * start_point[0]);
	cuboid.Surfparams.push_back(ycomponents);
	cuboid.complement.push_back(1);
	cuboid.surf_ids.push_back(3);
	y2components[7] = 1;
	y2components[9] = (-1 * (start_point[0] + ylength));
	cuboid.Surfparams.push_back(y2components);
	cuboid.complement.push_back(0);
	cuboid.surf_ids.push_back(4);

	zcomponents[8] = 1;
	zcomponents[9] = (-1 * start_point[0]);
	cuboid.Surfparams.push_back(zcomponents);
	cuboid.complement.push_back(1);
	cuboid.surf_ids.push_back(5);
	z2components[8] = 1;
	z2components[9] = (-1 * (start_point[0] + zlength));
	cuboid.Surfparams.push_back(z2components);
	cuboid.complement.push_back(0);
	cuboid.surf_ids.push_back(6);

	Complex_surf_input.push_back(cuboid);

	return;
}

void Input::surfSquarePrismInf(const std::string values, Surf_input prism)
{
	//type;ComplexID;Rank;a;b;parallel to;centre point;
	double centre_point[3], a, b;
	int parallel;
	prism.type = cubid;

	string::size_type split;
	string::size_type old;

	split = values.find(';');
	//comp ID
	prism.complex_id = stoi(values.substr(0, split));

	old = split + 1;
	split = values.find(';', old);
	//Subspace rank
	prism.subspace_rank = stoi(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//a
	a = stod(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//b
	b = stod(values.substr(old, split));

	old = split + 1;
	split = values.find(';', old);
	//parallel
	parallel = stoi(values.substr(old, split));

	string::size_type old_point, point;

	old = split + 1;
	split = values.find(';', old);
	old_point = old;

	for (int i = 0; i < 3; i++)
	{
		point = values.find(',', old_point);
		centre_point[i] = stod(values.substr(old_point, point));
		old_point = point + 1;
	}

	std::vector <double> acomponents(10, 0);
	std::vector <double> a2components(10, 0);
	std::vector <double> bcomponents(10, 0);
	std::vector <double> b2components(10, 0);

	if (parallel == 1)
	{
		acomponents[7] = 1;
		acomponents[9] = (-1 * (centre_point[1] - (a / 2)));
		prism.Surfparams.push_back(acomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(1);
		a2components[7] = 1;
		a2components[9] = (-1 * (centre_point[1] + (a / 2)));
		prism.Surfparams.push_back(a2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(2);

		bcomponents[8] = 1;
		bcomponents[9] = (-1 * (centre_point[2] - (b / 2)));
		prism.Surfparams.push_back(bcomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(3);
		b2components[8] = 1;
		b2components[9] = (-1 * (centre_point[2] + (b / 2)));
		prism.Surfparams.push_back(b2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(4);
	}
	else if (parallel == 2)
	{
		acomponents[8] = 1;
		acomponents[9] = (-1 * (centre_point[2] - (a / 2)));
		prism.Surfparams.push_back(acomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(1);
		a2components[8] = 1;
		a2components[9] = (-1 * (centre_point[2] + (a / 2)));
		prism.Surfparams.push_back(a2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(2);

		bcomponents[6] = 1;
		bcomponents[9] = (-1 * (centre_point[0] - (b / 2)));
		prism.Surfparams.push_back(bcomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(3);
		b2components[6] = 1;
		b2components[9] = (-1 * (centre_point[0] + (b / 2)));
		prism.Surfparams.push_back(b2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(4);

	}
	else
	{
		acomponents[6] = 1;
		acomponents[9] = (-1 * (centre_point[0] - (a / 2)));
		prism.Surfparams.push_back(acomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(1);
		a2components[6] = 1;
		a2components[9] = (-1 * (centre_point[0] + (a / 2)));
		prism.Surfparams.push_back(a2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(2);

		bcomponents[7] = 1;
		bcomponents[9] = (-1 * (centre_point[1] - (b / 2)));
		prism.Surfparams.push_back(bcomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(3);
		b2components[7] = 1;
		b2components[9] = (-1 * (centre_point[2] + (b / 2)));
		prism.Surfparams.push_back(b2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(4);
	}

	Complex_surf_input.push_back(prism);
	return;
}


//legacy
void Input::printData()
{
	cout << N << "\n";
	cout << M << "\n";
}
//legacy
int Input::checkInputCompletness()
{
	if (N == -999) {return 1;}
	if (M == -999) {return 1;}
	if (Needle == -999) { return 1; }
	if (Coin == -999) { return 1; }
	if (L == -999) { return 1; }
	if (l == -999) { return 1; }
	return 0;
}
