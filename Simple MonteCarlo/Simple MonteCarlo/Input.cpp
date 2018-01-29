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
		cellCreator(value);
		return;
	}
	if (para == "Subspace")
	{
		subspaceCreator(value);
		return;
	}
	if (para == "Transform")
	{
		transformCreator(value);
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
	
	return;

}

/*
A=0	B=1
C=2	D=3
E=4 F=5
G=6 H=7
I=8 J=9
*/

void Input::subspaceCreator(const std::string data)
{
	//ID;boundrycellID;Ranges
	subspace_input some_subspace;
	
	string::size_type split = data.find(';');
	string::size_type old;

	//ID
	some_subspace.subspace_id = stoi(data.substr(0, split));


	old = split + 1;
	split = data.find(';', old);
	//Cell ID
	some_subspace.boundery_cell_id = stoi(data.substr(old, split));

	string::size_type old_point, point;

	old = split + 1;
	split = data.find(';', old);
	old_point = old;

	
	point = data.find(',', old_point);
	some_subspace.subspacerange.x.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;
	point = data.find(',', old_point);
	some_subspace.subspacerange.x.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;

	point = data.find(',', old_point);
	some_subspace.subspacerange.y.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;
	point = data.find(',', old_point);
	some_subspace.subspacerange.y.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;

	point = data.find(',', old_point);
	some_subspace.subspacerange.z.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;
	point = data.find(',', old_point);
	some_subspace.subspacerange.z.push_back(stod(data.substr(old_point, point)));
	old_point = point + 1;
	
	


	subspace_input_data.push_back(some_subspace);
	return;
}

void Input::cellCreator(const std::string data)
{
	string::size_type split = data.find(';');
	string::size_type old;

	//ID;Name;Rank;Nrcomp;complements(surf1,complement1,surf2,complement2..);material;(laticetype;subspace inside;[bottom]/[hexa_pitch,centre];)
	Cell_input some_cell;

	//Name
	some_cell.cell_name = data.substr(0, split);
	


	old = split + 1;
	split = data.find(';', old);
	//Cell ID
	some_cell.Cell_id = stoi(data.substr(old, split));
	

	old = split + 1;
	split = data.find(';', old);
	//Subspace Rank
	some_cell.subspace_rank = stoi(data.substr(old, split));

	old = split + 1;
	split = data.find(';', old);
	//Subspace Rank
	int componentsNr = stoi(data.substr(old, split));


	string::size_type old_point, point;


	old = split + 1;
	split = data.find(';', old);
	old_point = old;

	for (int i = 0; i < componentsNr; i++)
	{
		cell_comp a_comp;
		point = data.find(',', old_point);
		a_comp.comp_surface_id = stoi(data.substr(old_point, point));
		old_point = point + 1;
		point = data.find(',', old_point);
		a_comp.cell_complement = stoi(data.substr(old_point, point));
		old_point = point + 1;
		some_cell.cell_complements.push_back(a_comp);
	}

	old = split + 1;
	split = data.find(';', old);

	some_cell.material = stoi(data.substr(old, split));

	old = split + 1;
	split = data.find(';', old);

	//Contains Latice (0 is no, 1 is square, 2 hexa)
	int latice = stoi(data.substr(old, split));

	some_cell.latice_info.type = latice;

	if (latice > 0)
	{

		old = split + 1;
		split = data.find(';', old);

		//subspace inside
		some_cell.latice_info.subspace_inside = stoi(data.substr(old, split));
		if (latice ==1)
		{
			old = split + 1;
			split = data.find(';', old);
			old_point = old;
			for (size_t i = 0; i < 3; i++)
			{
				point = data.find(',', old_point);
				some_cell.latice_info.bottom.push_back(stod(data.substr(old_point, point)));
				old_point = point + 1;
			}

		}
		else if (latice == 2)
		{
			old = split + 1;
			split = data.find(';', old);
			//hexa pitch
			some_cell.latice_info.hexa_pitch = stod(data.substr(old, split));
			old = split + 1;
			split = data.find(';', old);
			old_point = old;
			for (size_t i = 0; i < 3; i++)
			{
				point = data.find(',', old_point);
				some_cell.latice_info.bottom.push_back(stod(data.substr(old_point, point)));
				old_point = point + 1;
			}

		}
		else
		{
			some_cell.latice_info.hexa_pitch = -1;
		}
	}

	Cell_input_data.push_back(some_cell);
	return;
}

void Input::surfCreator(const std::string data)
{
	
	string type, values;

	string::size_type split = data.find(';');

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

	if (type == "Hex")
	{
		surfHexPrism(values, some_surface);
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
	//type;ComplexID;Rank;xlength;ylength;zlength;bottom point;
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
		b2components[9] = (-1 * (centre_point[1] + (b / 2)));
		prism.Surfparams.push_back(b2components);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(4);
	}

	Complex_surf_input.push_back(prism);
	return;
}


void Input::surfHexPrism(const std::string values, Surf_input prism)
{
	//type;ComplexID;Rank;height(0 is infinite);parallel to; centre point (at bottom if height); radius;
	double height, centre_point[3], radius;
	std::vector <bool> comp_cyl;
	int parallel_to;

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

	/*
	A=0	B=1
	C=2	D=3
	E=4 F=5
	G=6 H=7
	I=8 J=9
	*/

	if (parallel_to == 3)
	{
		/*
		std::vector <double> parpcomponents(10, 0);
		std::vector <double> parmcomponents(10, 0);

		parpcomponents[6] = 1;
		parmcomponents[6] = 1;

		parpcomponents[9] = (-1*(centre_point[0] + radius));
		parmcomponents[9] = (-1*(centre_point[0] - radius));

		prism.Surfparams.push_back(parpcomponents);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(1);
		prism.Surfparams.push_back(parmcomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(2);

		double tmppoint1[3];
		double tmppoint2[3];
		double tmppoint3[3];

		//tilted planes created with help of threepoint plane function
		
		//12 a clock
		tmppoint1[0] = centre_point[0];
		tmppoint1[1] = centre_point[1]+radius;
		tmppoint1[2] = centre_point[2];
		//12 a clock but one step up
		tmppoint2[0] = centre_point[0];
		tmppoint2[1] = centre_point[1] + radius;
		tmppoint2[2] = centre_point[2] + 1.0;
		//2 aclock
		tmppoint3[0] = centre_point[0] + (radius*(sqrt(3)/2));
		tmppoint3[1] = centre_point[1] + (radius*0.5);
		tmppoint3[2] = centre_point[2];

		std::vector <double> tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(3);

		//6 aclock
		tmppoint1[0] = centre_point[0];
		tmppoint1[1] = centre_point[1] - radius;
		tmppoint1[2] = centre_point[2];
		//6 aclock but one step up
		tmppoint2[0] = centre_point[0];
		tmppoint2[1] = centre_point[1] - radius;
		tmppoint2[2] = centre_point[2] + 1.0;
		//4 aclock
		tmppoint3[0] = centre_point[0] + (radius*(sqrt(3) / 2));
		tmppoint3[1] = centre_point[1] - (radius*0.5);
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(4);

		//6 aclock
		tmppoint1[0] = centre_point[0];
		tmppoint1[1] = centre_point[1] - radius;
		tmppoint1[2] = centre_point[2];
		//6 aclock but one step up
		tmppoint2[0] = centre_point[0];
		tmppoint2[1] = centre_point[1] - radius;
		tmppoint2[2] = centre_point[2] + 1.0;
		//8 aclock
		tmppoint3[0] = centre_point[0] - (radius*(sqrt(3) / 2));
		tmppoint3[1] = centre_point[1] - (radius*0.5);
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(5);

		//12 aclock
		tmppoint1[0] = centre_point[0];
		tmppoint1[1] = centre_point[1] + radius;
		tmppoint1[2] = centre_point[2];
		//12 aclock but one step up
		tmppoint2[0] = centre_point[0];
		tmppoint2[1] = centre_point[1] + radius;
		tmppoint2[2] = centre_point[2] + 1.0;
		//10 aclock
		tmppoint3[0] = centre_point[0] - (radius*(sqrt(3) / 2));
		tmppoint3[1] = centre_point[1] + (radius*0.5);
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(6);
		*/
		
		std::vector <double> parpcomponents(10, 0);
		std::vector <double> parmcomponents(10, 0);

		parpcomponents[7] = 1;
		parmcomponents[7] = 1;

		parpcomponents[9] = (-1 * (centre_point[1] + radius));
		parmcomponents[9] = (-1 * (centre_point[1] - radius));

		prism.Surfparams.push_back(parpcomponents);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(1);
		prism.Surfparams.push_back(parmcomponents);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(2);

		double tmppoint1[3];
		double tmppoint2[3];
		double tmppoint3[3];

		//tilted planes created with help of threepoint plane function

		//3 a clock
		tmppoint1[0] = centre_point[0] + radius;
		tmppoint1[1] = centre_point[1];
		tmppoint1[2] = centre_point[2];
		//3 a clock but one step up
		tmppoint2[0] = centre_point[0] + radius;
		tmppoint2[1] = centre_point[1];
		tmppoint2[2] = centre_point[2] + 1.0;
		//1 aclock
		tmppoint3[0] = centre_point[0]  + (radius*0.5);
		tmppoint3[1] = centre_point[1] + (radius*(sqrt(3) / 2));
		tmppoint3[2] = centre_point[2];

		std::vector <double> tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(3);


		//5 aclock
		tmppoint3[0] = centre_point[0] + (radius*0.5);
		tmppoint3[1] = centre_point[1] - (radius*(sqrt(3) / 2));
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(4);

		//9 aclock
		tmppoint1[0] = centre_point[0] - radius;
		tmppoint1[1] = centre_point[1];
		tmppoint1[2] = centre_point[2];
		//9 aclock but one step up
		tmppoint2[0] = centre_point[0] - radius;
		tmppoint2[1] = centre_point[1];
		tmppoint2[2] = centre_point[2] + 1.0;
		//7 aclock
		tmppoint3[0] = centre_point[0] - (radius*0.5);
		tmppoint3[1] = centre_point[1] - (radius*(sqrt(3) / 2));
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(1);
		prism.surf_ids.push_back(5);

		//11 aclock
		tmppoint3[0] = centre_point[0] - (radius*0.5);
		tmppoint3[1] = centre_point[1] + (radius*(sqrt(3) / 2));
		tmppoint3[2] = centre_point[2];

		tmp_parameters = generalPlanParams(tmppoint1, tmppoint2, tmppoint3);

		prism.Surfparams.push_back(tmp_parameters);
		prism.complement.push_back(0);
		prism.surf_ids.push_back(6);


		if (height != 0)
		{
			std::vector <double> botcomp(10, 0);
			std::vector <double> topcomp(10, 0);

			botcomp[8] = 1;
			botcomp[9] = (-1 * centre_point[2]);
			prism.Surfparams.push_back(botcomp);
			prism.complement.push_back(1);
			topcomp[8] = 1;
			topcomp[9] = (-1 * (centre_point[2] + height));
			prism.Surfparams.push_back(topcomp);
			prism.complement.push_back(0);
			prism.type = prism_hex;
		}
		else
		{
			prism.type = prism_hex_inf;
		}
		
	}


	Complex_surf_input.push_back(prism);

	return;
}

void Input::transformCreator(const std::string data)
{
	string::size_type split = data.find(';');
	string::size_type old;
	//Targets ID
	int targetID;
	targetID = stoi(data.substr(0, split));


	old = split + 1;
	split = data.find(';', old);
	//Target type
	int type;
	type = stoi(data.substr(old, split));

	old = split + 1;
	split = data.find(';', old);

	if (type == 1)
	{
		size_t i = 0;
		for (size_t j = 0; j < Complex_surf_input.size(); j++)
		{
			if (Cell_input_data[j].Cell_id == targetID)
			{
				i = j;
				break;
			}
		}
		//is moved 1 or 0(no)
		if (stoi(data.substr(old, split)) == 1)
		{
			Cell_input_data[i].moved = true;
			old = split + 1;
			split = data.find(';', old);
			//x
			Cell_input_data[i].place.push_back(stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//y
			Cell_input_data[i].place.push_back(stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//z
			Cell_input_data[i].place.push_back(stod(data.substr(old, split)));

					
		}
		old = split + 1;
		split = data.find(';', old);
		//rotated yes no
		if (stoi(data.substr(old, split)) == 1)
		{
			Cell_input_data[i].rotated = true;
			old = split + 1;
			split = data.find(';', old);
			//x
			Cell_input_data[i].angles.push_back(-1 * stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//y
			Cell_input_data[i].angles.push_back(-1 * stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//z
			Cell_input_data[i].angles.push_back(-1 * stod(data.substr(old, split)));
				
			
		}
	}
	else if (type == 2)
	{
		size_t i = 0;
		for (size_t j = 0; j < Complex_surf_input.size(); j++)
		{
			if (Complex_surf_input[j].complex_id == targetID)
			{
				i = j;
				break;
			}
		}
		if (stoi(data.substr(old, split)) == 1)
		{
			Complex_surf_input[i].moved = true;
			old = split + 1;
			split = data.find(';', old);
			//x
			Complex_surf_input[i].place.push_back(stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//y
			Complex_surf_input[i].place.push_back(stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//z
			Complex_surf_input[i].place.push_back(stod(data.substr(old, split)));

		}

		old = split + 1;
		split = data.find(';', old);

		if (stoi(data.substr(old, split)) == 1)
		{
			Complex_surf_input[i].rotated = true;
			old = split + 1;
			split = data.find(';', old);
			//x
			Complex_surf_input[i].angles.push_back(-1 * stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//y
			Complex_surf_input[i].angles.push_back(-1 * stod(data.substr(old, split)));

			old = split + 1;
			split = data.find(';', old);
			//z
			Complex_surf_input[i].angles.push_back(-1 * stod(data.substr(old, split)));
		}

	}
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

	return 0;
}
