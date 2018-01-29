#include "NData.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <random>
#include "stat.h"


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

void NData::BuildComps(std::string filename)
{
	std::string line;
	std::ifstream input(filename);

	std::vector <std::string> data;

	getline(input, line);

	int Comps_nr = stoi(line);

	for (int i = 0; i < Comps_nr; i++)
	{
		getline(input, line);

		Composit some_comp;

		std::string::size_type split = line.find(' ');
		std::string::size_type old;

		some_comp.Name = line.substr(0, split);

		old = split + 1;
		split = line.find(' ', old);

		int mix_nr = std::stoi(line.substr(old, split));

		for (int j = 0; j < mix_nr; j++)
		{
			getline(input, line);
			Mix_rates some_mix;
			std::string::size_type split = line.find(' ');
			std::string::size_type old;

			some_mix.Material_index = GetmaterialIndex(line.substr(0, split));

			old = split + 1;
			split = line.find(' ', old);

			some_mix.rate = std::stod(line.substr(old, split));
			some_comp.components.push_back(some_mix);
		}

		mixes.push_back(some_comp);
	}
	


	return;
}

void NData::AddGeom(Universe Geo, int top_index) 
{
	Geometry = Geo;
	top_uni = top_index;
	return;
}

void NData::CalcNuclearDensities()
{
	
}

double NData::EtoCrossUnknown(std::string Symbol, int MTn, double E)
{
	for (int i = 0; i < materials.size(); i++)
	{
		if (materials[i].Symbol == Symbol)
		{
			for (int j = 0; j < materials[i].MT_data.size(); j++)
			{
				if (materials[i].MT_data[j].MT_number == MTn)
				{
					std::map<double, double>::iterator itb = materials[i].MT_data[j].Ec_pair.begin();
						std::map<double, double>::iterator ite = materials[i].MT_data[j].Ec_pair.end();
					ite--;
					if (itb->first < E && E < ite->first)
					{
						std::map<double, double>::iterator it, it1;
						it = materials[i].MT_data[j].Ec_pair.lower_bound(E);

						it1 = it;
						it--;

						if (it != materials[i].MT_data[j].Ec_pair.begin())
						{
							//(sigma_i1-sigma_i)/(E_i1-E_i) *(E - E_i) + sigma_i
							double crosssection;
							crosssection = (((it1->second - it->second) / (it1->first - it->first))*(E - it->first)) + it->second;
							return crosssection;
						}
					}
				}
			}
		}
	}

	return -1;
}

double NData::EtoCross(int material_index, int mt_index, double E)
{
	std::map<double,double>::iterator itb = materials[material_index].MT_data[mt_index].Ec_pair.begin();
	std::map<double, double>::iterator ite =  materials[material_index].MT_data[mt_index].Ec_pair.end();
	ite--;
	if (itb->first < E && E < ite->first)
	{
		std::map<double, double>::iterator it, it1;
		it = materials[material_index].MT_data[mt_index].Ec_pair.lower_bound(E);

		it1 = it;
		it--;
		//(sigma_i1-sigma_i)/(E_i1-E_i) *(E - E_i) + sigma_i
		double crosssection;
		crosssection = (((it1->second - it->second) / (it1->first - it->first))*(E - it->first)) + it->second;
		return crosssection;

	}
	else
	{
		return 0;
	}
}


double NData::TotCross(int material_index, double E)
{
	double totcrosssection=0;
	for (int i = 0; i < materials[material_index].MT_data.size(); i++)
	{
		totcrosssection += EtoCross(material_index, i, E);
	}

	return totcrosssection;
}

int NData::GetmaterialIndex(std::string Symbol) 
{
	for (int i = 0; i < materials.size(); i++)
	{
		if (materials[i].Symbol == Symbol)
		{
			return i;
		}
	}
}

int NData::GetMixIndex(std::string Symbol)
{
	for (int i = 0; i < mixes.size(); i++)
	{
		if (mixes[i].Name == Symbol)
		{
			return i;
		}
	}
}

std::vector <int> NData::GetMt_index(int material_index)
{
	int size = materials[material_index].MT_data.size();
	std::vector <int> MT_index_pairs(size);

	for (int i = 0; i < size; i++)
	{
		MT_index_pairs[i] = materials[material_index].MT_data[i].MT_number;
	}
	return MT_index_pairs;
}

void NData::PrintLogLog(int material_index, int mt_index, bool total) 
{
	double LogMin = log10(1.005e-11 );
	double LogMax = log10(20 - 1e-8);
	int divider = 500;

	double delta = (LogMax - LogMin) / (divider - 1);

	std::vector <double> energypoints(500,LogMin);

	for (int i = 1; i < divider; i++)
	{
		energypoints[i] += delta*i;
	}

	for (int i = 0; i < divider; i++)
	{
		energypoints[i] = pow(10, energypoints[i]);
	}

	if (total)
	{
		std::string logfile;

		logfile = "loglog" + materials[material_index].Symbol + "_total" + ".txt";

		std::ofstream log;
		log.open(logfile, std::ios::out | std::ios::trunc);
		log << materials[material_index].Symbol + "\n";
		log << "Total cross section. \n";
		log << "****************************************\n";

		for (size_t i = 0; i < energypoints.size(); i++)
		{
			log << energypoints[i]
				<< ";"
				<< TotCross(material_index, energypoints[i])
				<< "\n";
		}
	}
	else
	{
		std::string logfile;

		logfile = "loglog" + materials[material_index].Symbol + "_" + std::to_string(mt_index) + ".txt";

		std::ofstream log;
		log.open(logfile, std::ios::out | std::ios::trunc);
		log << "Mt number: " << mt_index << "\n";
		log << "****************************************\n";
		
		for (size_t i = 0; i < energypoints.size(); i++)
		{
			log << energypoints[i]
				<< ";"
				<< EtoCross(material_index, mt_index, energypoints[i])
				<< "\n";
		}
	}

	return;
}

void NData::PrintinelasticLogLog(int material_index)
{
	std::vector <int> inelastics;
	std::vector <int> Mt_print;
	for (int i = 51; i < 92; i++)
	{
		for (int j = 0; j < materials[material_index].MT_data.size(); j++)
		{
			if (materials[material_index].MT_data[j].MT_number == i)
			{
				inelastics.push_back(j);
				Mt_print.push_back(i);
			}
		}
	}

	double LogMin = log10(1.005e-11);
	double LogMax = log10(20 - 1e-8);
	int divider = 500;

	double delta = (LogMax - LogMin) / (divider - 1);

	std::vector <double> energypoints(500, LogMin);

	for (int i = 1; i < divider; i++)
	{
		energypoints[i] += delta*i;
	}

	for (int i = 0; i < divider; i++)
	{
		energypoints[i] = pow(10, energypoints[i]);
	}

	std::string logfile;

	logfile = "loglog" + materials[material_index].Symbol + "_inelastic.txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << "Inelastics: \n";
	log << "****************************************\n";

	for (size_t i = 0; i < energypoints.size(); i++)
	{
		log << energypoints[i]
			<< ";";
		for (int j = 0; j < inelastics.size(); j++)
		{
			log << EtoCross(material_index, inelastics[j], energypoints[i]) << ";";
		}

		log << "\n";

	}

	log.close();

	logfile = "loglog" + materials[material_index].Symbol + "_inelasticPrint.txt";

	log.open(logfile, std::ios::out | std::ios::trunc);

	for (size_t i = 0; i < Mt_print.size(); i++)
	{
		log << Mt_print[i];

		log << "\n";

	}

	log.close();



}

double NData::SumMTsofE(int material_index, std::vector <int> Mts, double E)
{
	double cross = 0;
	for (int i = 0; i < Mts.size(); i++)
	{
		cross += EtoCross(material_index, Mts[i], E);
	}

	return cross;
}

void NData::PrintLogLogTotalComp(std::string Symbol)
{
	double LogMin = log10(1.005e-11);
	double LogMax = log10(20 - 1e-8);
	int divider = 500;

	double delta = (LogMax - LogMin) / (divider - 1);

	std::vector <double> energypoints(500, LogMin);

	for (int i = 1; i < divider; i++)
	{
		energypoints[i] += delta*i;
	}

	for (int i = 0; i < divider; i++)
	{
		energypoints[i] = pow(10, energypoints[i]);
	}

	int index = GetMixIndex(Symbol);

	std::string logfile;

	logfile = "loglog_" + mixes[index].Name + "_total" + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << materials[index].Symbol + "\n";
	log << "Total cross section. \n";
	log << "****************************************\n";

	for (size_t i = 0; i < energypoints.size(); i++)
	{
		log << energypoints[i] << ";";
		double cross = 0;
		for (int j = 0; j < mixes[index].components.size(); j++)
		{
			
			cross += (mixes[index].components[j].rate * TotCross(mixes[index].components[j].Material_index, energypoints[i]));

		}
		log << cross;
		log	<< "\n";
	}

}

void NData::ElasticReaction(int Neutron_index, int material_index)
{
	double E = Bank[Neutron_index].E.back();
	double Vt = 0;

	double Dirt[] = { 0,0,0 };
	double Vvt[] = { 0,0,0 };

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	if (materials[material_index].Z == 1 && E <= 200e-6)
	{
		double T_ev = materials[material_index].T;
		T_ev = T_ev * 8.6167e-11;

		double Et = SampleMaxwell(T_ev);

		Vt = sqrt((2 * Et * 1.602e-19) / (1.660e-27 * materials[material_index].Aw));
		IsotropicDirection(Dirt);
		Vvt[0] = Vt*Dirt[0];
		Vvt[1] = Vt*Dirt[1];
		Vvt[2] = Vt*Dirt[2];
	}

	//velocity lab
	double Vn = sqrt((2 * E * 1.602e-19) / (1.660e-27));
	//Vector velocity lab
	double Vvn[3];
	Vvn[0] = Vn*Bank[Neutron_index].Dir[0];
	Vvn[1] = Vn*Bank[Neutron_index].Dir[1];
	Vvn[2] = Vn*Bank[Neutron_index].Dir[2];
	//Vector central mass
	double Vvcm[3];
	Vvcm[0] = (Vvn[0] + materials[material_index].Aw * Vvt[0]) / (materials[material_index].Aw + 1);
	Vvcm[1] = (Vvn[1] + materials[material_index].Aw * Vvt[1]) / (materials[material_index].Aw + 1);
	Vvcm[2] = (Vvn[2] + materials[material_index].Aw * Vvt[2]) / (materials[material_index].Aw + 1);
	//Vector neutron C frame
	double VvnC[3];
	VvnC[0] = (round(Vvn[0]) - round(Vvcm[0]));
	VvnC[1] = (round(Vvn[1]) - round(Vvcm[1]));
	VvnC[2] = (round(Vvn[2]) - round(Vvcm[2]));

	//||V|| Netron Central frame 
	double VnC = sqrt(pow(VvnC[0],2) + pow(VvnC[1],2) + pow(VvnC[2],2));

	double DirC[] = { 0,0,0 };
	
	IsotropicDirection(DirC);
	// Vector new direction Central Frame 
	VvnC[0] = VnC * DirC[0];
	VvnC[1] = VnC * DirC[1];
	VvnC[2] = VnC * DirC[2];
	//Vector new direction
	double Vvn_prime[3];
	Vvn_prime[0] = VvnC[0] + Vvcm[0];
	Vvn_prime[1] = VvnC[1] + Vvcm[1];
	Vvn_prime[2] = VvnC[2] + Vvcm[2];
	// ||V_prime|| Lab
	double Vn_prime = sqrt(pow(Vvn_prime[0],2) + pow(Vvn_prime[1],2) + pow(Vvn_prime[2],2));
	double DirN[3];
	DirN[0] = Vvn_prime[0] / Vn_prime;
	DirN[1] = Vvn_prime[1] / Vn_prime;
	DirN[2] = Vvn_prime[2] / Vn_prime;

	double E_newJ = 0.5 * 1.660e-27*pow(Vn_prime,2);
	double E_new = E_newJ / 1.602e-19;
	Bank[Neutron_index].E.push_back(E_new);
	Bank[Neutron_index].speed = Vn_prime;
	Bank[Neutron_index].Dir[0] = DirN[0];
	Bank[Neutron_index].Dir[1] = DirN[1];
	Bank[Neutron_index].Dir[2] = DirN[2];

	return;

}

void NData::InElasticReactions(int Neutron_index, int material_index, double Q)
{
	double E = Bank[Neutron_index].E.back();
	double Vt = 0;

	double Dirt[] = { 0,0,0 };
	double Vvt[] = { 0,0,0 };

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	if (materials[material_index].Z == 1 && E <= 200e-6)
	{
		double T_ev = materials[material_index].T;
		T_ev = T_ev * 8.6167e-11;

		double Et = SampleMaxwell(T_ev);

		Vt = sqrt((2 * Et * 1.602e-19) / (1.660e-27 * materials[material_index].Aw));
		IsotropicDirection(Dirt);
		Vvt[0] = Vt*Dirt[0];
		Vvt[1] = Vt*Dirt[1];
		Vvt[2] = Vt*Dirt[2];
	}

	//velocity lab
	double Vn = Bank[Neutron_index].speed;
	//Vector velocity lab
	double Vvn[3];
	Vvn[0] = Vn*Bank[Neutron_index].Dir[0];
	Vvn[1] = Vn*Bank[Neutron_index].Dir[1];
	Vvn[2] = Vn*Bank[Neutron_index].Dir[2];
	//Vector central mass
	double Vvcm[3];
	Vvcm[0] = (Vvn[0] + materials[material_index].Aw * Vvt[0]) / (materials[material_index].Aw + 1);
	Vvcm[1] = (Vvn[1] + materials[material_index].Aw * Vvt[1]) / (materials[material_index].Aw + 1);
	Vvcm[2] = (Vvn[2] + materials[material_index].Aw * Vvt[2]) / (materials[material_index].Aw + 1);
	//Vector neutron C frame
	double VvnC[3];
	VvnC[0] = (round(Vvn[0]) - round(Vvcm[0]));
	VvnC[1] = (round(Vvn[1]) - round(Vvcm[1]));
	VvnC[2] = (round(Vvn[2]) - round(Vvcm[2]));

	//||V||^2 Netron Central frame 
	double VnC = pow(VvnC[0], 2) + pow(VvnC[1], 2) + pow(VvnC[2], 2);
	// ||V'|| 
	VnC = sqrt(VnC + (2 * materials[material_index].Aw * Q * 1.602e-19) / ((materials[material_index].Aw + 1)*1.660e-27));

	double DirC[] = { 0,0,0 };

	IsotropicDirection(DirC);
	// Vector new direction Central Frame 
	VvnC[0] = VnC * DirC[0];
	VvnC[1] = VnC * DirC[1];
	VvnC[2] = VnC * DirC[2];
	//Vector new direction
	double Vvn_prime[3];
	Vvn_prime[0] = VvnC[0] + Vvcm[0];
	Vvn_prime[1] = VvnC[1] + Vvcm[1];
	Vvn_prime[2] = VvnC[2] + Vvcm[2];
	// ||V_prime|| Lab
	double Vn_prime = sqrt(pow(Vvn_prime[0], 2) + pow(Vvn_prime[1], 2) + pow(Vvn_prime[2], 2));
	double DirN[3];
	DirN[0] = Vvn_prime[0] / Vn_prime;
	DirN[1] = Vvn_prime[1] / Vn_prime;
	DirN[2] = Vvn_prime[2] / Vn_prime;

	double E_newJ = 0.5 * 1.660e-27*pow(Vn_prime, 2);
	double E_new = E_newJ / 1.602e-19;
	Bank[Neutron_index].E.push_back(E_new);
	Bank[Neutron_index].speed = Vn_prime;
	Bank[Neutron_index].Dir[0] = DirN[0];
	Bank[Neutron_index].Dir[1] = DirN[1];
	Bank[Neutron_index].Dir[2] = DirN[2];

	return;
}

void NData::FissionReaction(int Neutron_index, int material_index)
{
	Bank[Neutron_index].died = true;
	Bank[Neutron_index].cause_of_death = 2;

	double E = Bank[Neutron_index].E.back();

	int born;

	std::map<double, double>::iterator it, it1;
	it = materials[material_index].Fission_nubar.lower_bound(E);

	it1 = it;
	it--;
	//(sigma_i1-sigma_i)/(E_i1-E_i) *(E - E_i) + sigma_i
	double nubar;
	nubar = (((it1->second - it->second) / (it1->first - it->first))*(E - it->first)) + it->second;
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	if (dis(gen) < (nubar - floor(nubar)))
	{
		born = ceil(nubar);
	}
	else
	{
		born = floor(nubar);
	}

	Bank[Neutron_index].children = born;

	for (int i = 0; i < born; i++)
	{
		Neutron newN;
		newN.E.push_back(SampleMaxwell(1.2895));
		newN.speed = sqrt((2 * newN.E.back() * 1.602e-19) / (1.660e-27));
		newN.time = Bank[Neutron_index].time;
		newN.generation = Bank[Neutron_index].generation + 1;
		double dir[] = { 0,0,0 };
		IsotropicDirection(dir);
		newN.Pos[0] = Bank[Neutron_index].Pos[0];
		newN.Pos[1] = Bank[Neutron_index].Pos[1];
		newN.Pos[2] = Bank[Neutron_index].Pos[2];
		Bank[Neutron_index].offspring.push_back(newN);
	}

	return;

}

void NData::CaptureReaction(int Neutron_index)
{
	Bank[Neutron_index].died = true;
	Bank[Neutron_index].cause_of_death = 1;
	return;
}

double NData::SampleMaxwell(double T_ev)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	double shi1, shi2, shi3, shi4, R, Et;
	do
	{
		shi1 = dis(gen);
		shi2 = dis(gen);
		R = (pow(shi1,2) + pow(shi2,2));
	} while (R > 1);

	shi3 = dis(gen);
	shi4 = dis(gen);

	Et = -T_ev * (pow(shi1,2)*log(shi3) / R + log(shi4));

	return Et;
}

double NData::SampleLength(double crosstot)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	return -1 * log(dis(gen)) / crosstot;
}

void NData::FindReaction(int Material_index, int neutron_index)
{
	double E = Bank[neutron_index].E.back();

	//Calculate total crosssection of composition.
	double total_cross = TotCross(Material_index,E);
	std::vector <double> Cross_N;
	for (int i = 0; i < materials[Material_index].MT_data.size(); i++)
	{
		Cross_N.push_back(EtoCross(Material_index,i,E));
	}

	std::vector<std::pair <int,double> > Rates(1);
	Rates[0].second = 0.0;
	Rates[0].first = -1;
	for (int i = 0; i < Cross_N.size(); i++)
	{
		std::pair <int, double> tmp;
		double tmpCross = Cross_N[i] / total_cross;
		tmpCross = tmpCross + Rates[i].second;
		tmp.second = tmpCross;
		tmp.first = i;
		Rates.push_back(tmp);
	}

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	//dis(gen)

	double choosyboi = dis(gen);

	//Choose Mt
	int ReactionMt;
	int Index_inelastic;
	for (int i = 1; i < Rates.size(); i++)
	{
		if (choosyboi <= Rates[i].second)
		{
			ReactionMt = materials[Material_index].MT_data[Rates[i].first].MT_number;
			Index_inelastic = i;
			break;
		}
	}
	//Move to reaction based on Mt-nr
	if (ReactionMt == 2)
	{
		ElasticReaction(neutron_index, Material_index);
	}
	else if (ReactionMt == 18)
	{
		FissionReaction(neutron_index, Material_index);
	}
	else if (51 <= ReactionMt && ReactionMt <= 90)
	{
		double Q = materials[Material_index].MT_data[Rates[Index_inelastic].first].Q_value;
		InElasticReactions(neutron_index, Material_index, Q);
	}
	else if (ReactionMt == 91)
	{
		//No, idea if same 51-90
	}
	else if (102 <= ReactionMt && ReactionMt <= 107)
	{
		CaptureReaction(neutron_index);
	}

	return;

}

void NData::FindNucleus(int Mix_nr, int neutron_index)
{

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	//dis(gen)

	double E = Bank[neutron_index].E.back();
	if (E < 9.9E-12 )
	{
		Bank[neutron_index].died = true;
		return;
	}
	//Calculate total crosssection of composition.
	double total_cross=0;
	std::vector <double> Cross_N;
	for (int i = 0; i < mixes[Mix_nr].components.size(); i++)
	{
		double tmp = mixes[Mix_nr].components[i].rate * TotCross(mixes[Mix_nr].components[i].Material_index, E);
		total_cross += tmp;
		Cross_N.push_back(tmp);
	}
	

	double track_length = SampleLength(total_cross);
	//Here needs to be check for surface changes...
	Bank[neutron_index].path_length += track_length;
	double V = sqrt((2 * E * 1.602e-19) / (1.660e-27));
	Bank[neutron_index].time += track_length/ V;

	//Calculate Rates
	std::vector<double> Rates;
	Rates.push_back(0.0);
	for (int i = 0; i < Cross_N.size(); i++)
	{
		double tmp;
		tmp = Cross_N[i] / total_cross;
		Rates.push_back(tmp+Rates.back());
	}



	double choosyboi = dis(gen);

	for (int i = 1; i <= Rates.size(); i++)
	{
		if (choosyboi < Rates[i])
		{
			FindReaction(mixes[Mix_nr].components[i-1].Material_index, neutron_index);
			return;
		}
	}


}

void NData::FindNucleusGeometry(int neutron_index)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);
	//dis(gen)

	double E = Bank[neutron_index].E.back();
	if (E < 9.9E-12)
	{
		Bank[neutron_index].died = true;
		return;
	}
	//Calculate total crosssection of composition.
	double total_cross = 0;

	Geometry.CellInUniverse(Bank[neutron_index].Pos, top_uni, Bank[neutron_index].target);
	std::string material_tmp= Geometry.subspaces[Bank[neutron_index].target[0]].printCellName(Bank[neutron_index].target[0]);
	
	//Escaped neutrons are dead
	if (material_tmp == "Outside")
	{
		Bank[neutron_index].died = true;
		Bank[neutron_index].cause_of_death = 0;
		return;
	}

	int index = GetMixIndex(material_tmp);

	std::vector <double> Cross_N;
	for (int i = 0; i < mixes[index].components.size(); i++)
	{
		double tmp = mixes[index].components[i].rate * TotCross(mixes[index].components[i].Material_index, E);
		total_cross += tmp;
		Cross_N.push_back(tmp);
	}

	double track_length = SampleLength(total_cross);
	double track_length_to_border = Geometry.LineLengthUniverse(Bank[neutron_index].Pos, Bank[neutron_index].Dir, top_uni);
	double V = sqrt((2 * E * 1.602e-19) / (1.660e-27));
	if (track_length >= track_length_to_border )
	{
		double epsilon = 1e-8;
		Bank[neutron_index].path_length += (track_length_to_border + epsilon);
		Bank[neutron_index].time += track_length_to_border / V;
		MoveNeutron(neutron_index, track_length_to_border);
		return;
	}
	
	Bank[neutron_index].path_length += track_length;
	Bank[neutron_index].time += track_length / V;
	MoveNeutron(neutron_index, track_length);

	std::vector<double> Rates;
	Rates.push_back(0.0);
	for (int i = 0; i < Cross_N.size(); i++)
	{
		double tmp;
		tmp = Cross_N[i] / total_cross;
		Rates.push_back(tmp + Rates.back());
	}



	double choosyboi = dis(gen);

	for (int i = 1; i <= Rates.size(); i++)
	{
		if (choosyboi < Rates[i])
		{
			FindReaction(mixes[index].components[i - 1].Material_index, neutron_index);
			return;
		}
	}
	

}

void NData::PopulateBank(int Amount, double startingE)
{
	Neutron tmp;
	tmp.E.push_back(startingE);
	double dir[] = { 0,0,0 };
	int tmp_pos[] = { 0,0,0 };
	tmp.Pos[0] = tmp_pos[0];
	tmp.Pos[1] = tmp_pos[1];
	tmp.Pos[2] = tmp_pos[2];
	tmp.speed = sqrt((2 * startingE * 1.602e-19) / (1.660e-27));
	tmp.time = 0.0;
	tmp.generation = 0;

	for (int i = 0; i < Amount; i++)
	{
		IsotropicDirection(dir);

		tmp.Dir[0] = dir[0];
		tmp.Dir[1] = dir[1];
		tmp.Dir[2] = dir[2];
		Bank.push_back(tmp);
	}
	return;
}

void NData::IsotropicDirection(double * dir)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(0, 1);

	double omega = dis(gen) * 2 - 1;
	double phi = dis(gen) * 2 * 3.14159265359;
	
	dir[0] = sqrt(1 - pow(omega, 2))*cos(phi);
	dir[1] = sqrt(1 - pow(omega, 2))*sin(phi);
	dir[2] = omega;

	return;
}


void NData::MoveNeutron(int neutron_index, double length)
{
	Bank[neutron_index].Pos[0] = Bank[neutron_index].Pos[0] + (length * Bank[neutron_index].Dir[0]);
	Bank[neutron_index].Pos[1] = Bank[neutron_index].Pos[1] + (length * Bank[neutron_index].Dir[1]);
	Bank[neutron_index].Pos[2] = Bank[neutron_index].Pos[2] + (length * Bank[neutron_index].Dir[2]);
	
	return;
}

void NData::Teller()
{
	std::vector <Neutron> tmp_bank;

	for (int i = 0; i < Bank.size(); i++)
	{
		if (Bank[i].children != 0)
		{
				tmp_bank.insert(tmp_bank.end(), Bank[i].offspring.begin(), Bank[i].offspring.end());
		}
	}

	Graveyard.insert(Graveyard.end(), Bank.begin(), Bank.end());

	EmptyBank();

	Bank = tmp_bank;

	return;
}

void NData::EmptyBank()
{
	Bank.erase(Bank.begin(), Bank.end());
}

void NData::SlowdownPlot(int N, int index_comp) 
{

	PopulateBank(N, 14);
	int material_index = mixes[index_comp].components[0].Material_index;

	#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < Bank.size(); i++)
	{
		for (int j = 0; j < 200; j++)
		{
			ElasticReaction(i, material_index);
		}
	}

	std::vector <double> means;
	std::vector <double> StdDiv;

	for (int i = 0; i < 200; i++)
	{
		std::vector <double> Energies;
		for (int j = 0; j < Bank.size(); j++)
		{
			Energies.push_back(Bank[j].E[i]);
		}

		means.push_back(getMean(Energies));
		StdDiv.push_back(getStddiv(Energies,means.back()));
	}

	std::string logfile;

	logfile = mixes[index_comp].Name + "_slowdown" + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << materials[index_comp].Symbol + "\n";
	log << "mean;StdDiv \n";
	log << "****************************************\n";

	for (size_t i = 0; i < means.size(); i++)
	{
		log << means[i] << ";";
		log << StdDiv[i];
		log << "\n";
	}


	EmptyBank();
	return;
}

void NData::SlowdownPlotInelastic(int N, int index_comp)
{

	PopulateBank(N, 14);
	

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < Bank.size(); i++)
	{
		for (int j = 0; j < 200; j++)
		{
			if (!Bank[i].died)
			{
				FindNucleus(index_comp, i);
			}
		}
	}

	std::vector <double> means;
	std::vector <double> StdDiv;

	for (int i = 0; i < 50; i++)
	{
		std::vector <double> Energies;
		for (int j = 0; j < Bank.size(); j++)
		{
			if (i < Bank[j].E.size())
			{
				Energies.push_back(Bank[j].E[i]);
			}
		}

		means.push_back(getMean(Energies));
		StdDiv.push_back(getStddiv(Energies, means.back()));
	}

	std::string logfile;

	logfile = mixes[index_comp].Name + "_slowdown_inelastic" + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << materials[index_comp].Symbol + "\n";
	log << "mean;StdDiv \n";
	log << "****************************************\n";

	for (size_t i = 0; i < means.size(); i++)
	{
		log << means[i] << ";";
		log << StdDiv[i];
		log << "\n";
	}


	EmptyBank();
	return;
}

void NData::Multiplication(int N, int index_comp)
{


	PopulateBank(N, 1);

#pragma omp parallel for schedule(dynamic,1)
	for (int i = 0; i < N; i++)
	{
		int Guard = 0;
		while (!Bank[i].died && Guard < 1000)
		{
			FindNucleus(index_comp, i);
			Guard++;
		}
	}

	std::vector <double> children;

	for (int i = 0; i < N; i++)
	{
		children.push_back(Bank[i].children);
	}

	double mean = getMean(children);
	double StdDiv = getStddiv(children, mean);

	Teller();

	std::string logfile;

	logfile = mixes[index_comp].Name + "_test" + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << materials[index_comp].Symbol + "\n";
	log << "mean;StdDiv \n";
	log << mean << ";" << StdDiv;
	
	return;
}

void NData::TestEnrichment(int N, int index_comp)
//Uranium must be first in Composit and must be so that 235 is first
{
	std::vector <double> old;
	for (size_t i = 0; i < 2; i++)
	{
		old.push_back(mixes[index_comp].components[i].rate);
	}
	
	std::vector <double> means;
	std::vector <double> StdDiv;
	std::vector <double> Enrichment;

	while (mixes[index_comp].components[0].rate/mixes[index_comp].components[1].rate < 0.20)
	{
		PopulateBank(N, 1);

		#pragma omp parallel for schedule(dynamic,1)
		for (int i = 0; i < N; i++)
		{
			int Guard = 0;
			while (!Bank[i].died && Guard < 1000)
			{
				FindNucleus(index_comp, i);
				Guard++;
			}
		}

		std::vector <double> children;

		for (int i = 0; i < N; i++)
		{
			children.push_back(Bank[i].children);
		}

		means.push_back(getMean(children));
		StdDiv.push_back(getStddiv(children, means.back()));
		Enrichment.push_back(mixes[index_comp].components[0].rate);
		std::cout << Enrichment.back() << "\n";
		EmptyBank();

		mixes[index_comp].components[0].rate = mixes[index_comp].components[0].rate + 0.005;
		mixes[index_comp].components[1].rate = mixes[index_comp].components[1].rate - 0.005;

	}

	for (size_t i = 0; i < 2; i++)
	{
		mixes[index_comp].components[i].rate = old[i];
	}

	


	std::string logfile;

	logfile = mixes[index_comp].Name + "_Enrichmenttest" + ".txt";

	std::ofstream log;
	log.open(logfile, std::ios::out | std::ios::trunc);
	log << materials[index_comp].Symbol + "\n";
	log << "mean;StdDiv;Enrichment \n";
	for (int i = 0; i < means.size(); i++)
	{
		log << means[i] << ";" << StdDiv[i] << ";" << Enrichment[i] <<"\n";
	}
	return;
}

void NData::TestGeo(int N)
{
	PopulateBank(N, 1);
	
	while (Bank[0].died == false)
	{
		FindNucleusGeometry(0);
	}

	int jumps = Bank[0].E.size();
	return;
}

void NData::ExternalSource()
{
	int M = 10;
	int N = 100;

	
	for (int k = 0; k < M; k++)
	{
		PopulateBank(N, 1);
		int Guard_out = 0;
		do
		{

			#pragma omp parallel for schedule(dynamic,1)
			for (int i = 0; i < N; i++)
			{
				int Guard = 0;
				while (!Bank[i].died && Guard < 1000)
				{
					FindNucleusGeometry(i);
					Guard++;
				}
			}
			Teller();
			Guard_out++;
			if (Bank.size() == 0)
			{
				break;
			}
		} while (Guard_out < 5);
		EmptyBank();
	}

	return;

}

