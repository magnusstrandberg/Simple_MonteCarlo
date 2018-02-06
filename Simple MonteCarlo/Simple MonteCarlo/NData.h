#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "universe.h"

struct MT
{
	int MT_number;
	double Q_value;
	std::map<double, double> Ec_pair;
};

struct materialdata
{
	int Z, A;
	double T, Aw, Nuclear_density;
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

struct Track_length
{
	int target[2]; //for cell
	double path_length;
	double macrocross;
	int material_index;
	double neutron_dens;
	//0 escape, 1 surface cross, Mt-nr 
	int reaction_type;

};

struct Neutron 
{
	std::vector <double> E;
	double Dir[3];
	double Pos[3];
	std::vector <std::vector <double> > Position_history;
	bool died = false;
	// 0 escape,1 captured,2fission
	int cause_of_death;
	int target[2];
	int children = 0;
	bool parent_was_thermal = false;
	//reaction in second,0 escape,1 captured,2fission,3elastic,4inellastic
	//What hit What Happen
	std::vector <std::pair <int, int> > whwh;
	std::vector <Track_length> TLE;
	std::vector <Neutron> offspring;
	double time,speed, time_born;
	double path_length = 0;
	int generation;
};



class NData
{
public:
	NData();
	~NData();
	void readdatafile(std::string filename, int it);
	void readlist(std::string filename);
	void BuildComps(std::string filename);
	void AddGeom(Universe Geo, int top_index);
	void CalcNuclearDensities();
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

	double SampleLength(double crosstot);

	void FindReaction(int Material_index, int neutron_index);
	void FindNucleus(int Mix_nr, int neutron_index);
	void FindNucleusGeometry(int neutron_index);
	void PopulateBank(int Amount, double startingE);

	void IsotropicDirection(double * dir);

	void MoveNeutron(int neutron_index, double length);

	void Teller();

	void SlowdownPlot(int M, int index_comp);

	void SlowdownPlotInelastic(int N, int index_comp);

	void Multiplication(int N, int M, int index_comp);

	void TestEnrichment(int N, int M,int index_comp);

	void FissionTime(int N, int index_comp);

	void fourfactor(int N, int M, int index_comp);

	void TestGeo(int N);

	void Tally();

	void Analog(std::vector<std::vector<std::vector<Track_length>>> data, int N);

	void TrackLengthEstimator(std::vector<std::vector<std::vector<Track_length>>> data, int N);

	void ExternalSource(int N, int M);


	void EmptyBank();

	void EmptyGraveyard();
	
	
	std::vector <materialdata> materials;
	std::vector <Composit> mixes;
	std::vector <Neutron> Bank;
	std::vector <Neutron> Graveyard;
	std::vector <std::vector <std::vector <Track_length> > > TLE_M;
	Universe Geometry;
	int top_uni;
};

