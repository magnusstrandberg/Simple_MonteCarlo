#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <map>

struct MT
{
	int MT_number;
	double Q_value;
	std::map<double, double> Ec_pair;
};

struct materialdata
{
	int Z, A;
	double T, Aw;
	std::string Symbol;
	int Fissionable;
	std::map<double, double> Fission_nubar;
	std::vector <MT> MT_data;

};

struct Mix_rates
{
	int Material_index;
	double rate;
};

struct Composit
{
	std::string Name;
	std::vector <Mix_rates> components;
};

struct Neutron 
{
	std::vector <double> E;
	double Dir[3];
	bool died = false;
	int children = 0;
};



class NData
{
public:
	NData();
	~NData();
	void readdatafile(std::string filename, int it);
	void readlist(std::string filename);
	void BuildComps(std::string filename);
	double EtoCrossUnknown(std::string Symbol, int MTn, double E);
	double EtoCross(int material_index, int mt_index, double E);
	double TotCross(int material_index, double E);
	int GetmaterialIndex(std::string Symbol);
	int GetMixIndex(std::string Symbol);
	std::vector<int> GetMt_index(int material_index);

	void PrintLogLog(int material_index, int mt_index, bool total);
	void PrintinelasticLogLog(int material_index);
	double SumMTsofE(int material_index, std::vector<int> Mts, double E);
	void PrintLogLogTotalComp(std::string Symbol);

	void ElasticReaction(int Neutron_index, int material_index);
	void InElasticReactions(int Neutron_index, int material_index, double Q);
	void FissionReaction(int Neutron_index, int material_index);
	void CaptureReaction(int Neutron_index);

	double SampleMaxwell(double T_ev);

	void FindReaction(int Material_index, int neutron_index);
	void FindNucleus(int Mix_nr, int neutron_index);
	void PopulateBank(int Amount, double startingE);

	void IsotropicDirection(double * dir);

	void Teller(int Neutron_index, int Material_index);
	
	
	std::vector <materialdata> materials;
	std::vector <Composit> mixes;
	std::vector <Neutron> Bank;
	std::vector <Neutron> Graveyard;
};

