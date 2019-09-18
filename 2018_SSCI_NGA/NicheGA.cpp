/**
*  Copyright (c) 2018 SCUT
*  All rights reserved.
*
*  @author      Ting Huang (cshting@mail.scut.edu.cn)
*  @date		2018-03-31
*
*  @brief       The niching genetic algorithm (NMA) is proposed to solve 25 MSTSPs.
*
*  @version
*    - V4
*
*  @command
*		[begin Index] [end index] [runs] [NP] [crossover rate] [mutation rate] [dir] [selection strategy] [crossover strategey] [mutation strategy]
*		 0 24 50 150 0.9 0.01 .\Results roulette PMX RSM 6
*
*	@reference
*		T. Huang, Y. Gong and J. Zhang, “Seeking Multiple Solutions of Combinatorial optimization Problems: A Proof of Principle Study,” 2018 IEEE Symposium Series on Computational Intelligence (SSCI), Bangalore, India, 2018, pp. 1212-1218.
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <limits>
#include <sstream>  
#include <string>  
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cassert>

using namespace std;
vector < vector<int> > dist;
vector<vector<int> > disimilar_dist_matrix;
vector<vector<int> > share_dist_matrix;
vector<vector<int> > share_dist_matrix2;
vector<vector<int> > share_matrix;
char file_name_char[128];
void convert_file_name(const string &);
string input_dir;
string output_dir;
int MIN_LEN;
long long FES = 0;
class GeneticAlgorithm;
class Problem;
class Chromosome;


int rndint(int Lbound, int Ubound)
{
	return (rand() % (Ubound - Lbound)) + Lbound;
}

double rnddouble(double Lbound, double Ubound)
{
	return rand() / (double)(RAND_MAX / (Ubound - Lbound)) + Lbound;
}
class cardinate
{
public:
	int x;
	int y;
};
struct index_len
{
public:
	int index;
	int len;
};
bool sort_index_asc_second(const index_len& a, const index_len & b)
{
	return a.len < b.len;
}
class City
{
public:
	vector<cardinate> city_vec;
	int city_num;
public:
	int distance(cardinate, cardinate);
	City(const string &problem_name);
	City() :city_num(0) {};
};
City::City(const string &problem_name) {
	string file_string = input_dir + problem_name + ".tsp";
	ifstream ifile_tsp;
	convert_file_name(file_string);
	ifile_tsp.open(file_name_char, ios::in);
	if (!ifile_tsp)
		exit(-1);
	city_vec.clear();
	cardinate cardinate_tmp;
	while (!ifile_tsp.eof())
	{
		ifile_tsp >> cardinate_tmp.x >> cardinate_tmp.y;
		city_vec.push_back(cardinate_tmp);
	}
	city_vec.pop_back();
	city_num = city_vec.size();
	ifile_tsp.close();

	dist.clear();
	vector<int> dist_vec(city_num, 0);
	for (int i = 0; i < city_num; i++)
		dist.push_back(dist_vec);
	for (int i = 0; i < city_num; i++)
	{
		for (int j = 0; j < city_num; j++)
			dist[i][j] = dist[j][i] = distance(city_vec[i], city_vec[j]);
	}
}
int City::distance(cardinate a, cardinate b)
{
	return  (int)(sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)) + 0.5);
}

class Chromosome
{
public:
	vector<int> path;
	int length;
	double fitness;
public:
	Chromosome() :length(0) {};
	Chromosome(vector<int> p);
	Chromosome(const Chromosome & b);

	void evaluate();
	int distance(int a, int b) { return  (int)dist[a][b]; };
	friend class Problem;
};
bool len_descend(const Chromosome& a, const Chromosome& b)
{
	return a.length < b.length;
}
bool operator<(const Chromosome& c1, const Chromosome& c2)
{
	return c1.fitness < c2.fitness;
}
bool operator==(const Chromosome& c1, const Chromosome& c2)
{
	return c1.length == c2.length;
}
void Chromosome::evaluate()
{
	length = 0;
	for (int i = 1; i < path.size(); i++)
		length += distance(path[i], path[i - 1]);
	length += distance(path[path.size() - 1], path[0]);
	fitness = length;
	FES++;

}
Chromosome::Chromosome(vector<int> p)
{
	path = p;
	length = 0;
	fitness = 0.0;
}
Chromosome::Chromosome(const Chromosome & b)
{
	path = b.path;
	length = b.length;
	fitness = b.fitness;
}

class GeneticAlgorithm
{
private:
	vector<Chromosome> chroms;
	City city_list;
	int NP;
	int current_NP;
	int city_num;
	int CF = 5;

	long long maxFES;
	double mutation_ratio;
	double crossover_ratio;
	int max_iterations;

public:

	int iter;
	friend class Problem;
	GeneticAlgorithm(int np, int cn, double cr, double mr, int MAX_ITER, City c, long long maxfes);
	GeneticAlgorithm(const GeneticAlgorithm & b);

	void initiate();
	void speciation(const GeneticAlgorithm &, const int species_size);
	void crossover();
	void crossover(const string &);
	void mutation(const string &);
	void mutation();
	void pass_next_generation_speciation_V2(GeneticAlgorithm &, const int &);
	void length_evaluate();

	void remove_redudant_report2_file_V3(GeneticAlgorithm &, double threshold, const int &run, const int &track_fes, const string &problem_name, const int species_size);

	bool sort_cmp_len(const vector<Chromosome> &a, const vector<Chromosome> &b) {
		return a[0].length < b[0].length;
	}
	bool terminate() { return FES >= maxFES; };
	int share_distance(const vector<int>& a, const vector<int>& b);
	int share_distance(const Chromosome& a, const Chromosome& b);

	void report();
	void report_2_file(const int &run, const string &problem_name);
};

class index_value
{
public:
	int index;
	int value;
};
bool sort_index_asc_value(const index_value &a, const index_value &b)
{
	return a.value < b.value;
}
void GeneticAlgorithm::remove_redudant_report2_file_V3(GeneticAlgorithm & Parent, double threshold, const int & run, const int &track_fes, const string & problem_name, int species_size)
{
	vector<Chromosome> chroms_no_redudent;
	vector<Chromosome> chroms_tmp = chroms;
	int min_len_tmp = chroms_tmp[0].length;
	int min_index = 0;
	for (int i = 0; i < chroms_tmp.size(); i++)
	{
		if (chroms_tmp[i].length < min_len_tmp)
		{
			min_len_tmp = chroms_tmp[i].length;
			min_index = i;
		}
	}
	chroms_no_redudent.push_back(chroms[min_index]);
	int threshold_min_len = min_len_tmp * (1 + threshold);


	vector<index_value> index_tmp_vec(species_size);
	index_value index_tmpp;
	int species_size_tmp = species_size;
	for (int i = 0; i < chroms_tmp.size(); i += species_size_tmp)
	{
		if (i + species_size_tmp >= chroms_tmp.size())
			species_size_tmp = chroms_tmp.size() - i;
		for (int j = i; j < i + species_size_tmp; j++)
		{
			index_tmp_vec[j - i].index = j;
			index_tmp_vec[j - i].value = chroms_tmp[j].length;
		}

		sort(index_tmp_vec.begin(), index_tmp_vec.end(), sort_index_asc_value);



		for (int j = 0; j < species_size_tmp; j++)
		{
			bool exist_flag = false;
			int max_share_dist = -999999;
			for (int k = 0; k < chroms_no_redudent.size(); k++)
			{
				int share_dist_tmp = share_distance(chroms_tmp[index_tmp_vec[j].index], chroms_no_redudent[k]);
				if (max_share_dist < share_dist_tmp)
					max_share_dist = share_dist_tmp;
				if (city_num == share_dist_tmp)
				{
					exist_flag = true;
					break;
				}
			}
			if (exist_flag == true)
				continue;
			if (chroms_tmp[index_tmp_vec[j].index].length == min_len_tmp)
			{
				chroms_no_redudent.push_back(chroms_tmp[index_tmp_vec[j].index]);
			}
			else if (chroms_tmp[index_tmp_vec[j].index].length < threshold_min_len && max_share_dist < city_num - 2 * log(city_num))
			{
				chroms_no_redudent.push_back(chroms_tmp[index_tmp_vec[j].index]);
			}
			else {
				break;
			}
		}
	}
	string file_string = output_dir + problem_name + "_RUN_" + to_string(run) + ".alg_solution";
	ofstream ofile_record;
	convert_file_name(file_string);
	ofile_record.open(file_name_char, ios::out);
	for (int i = 0; i < chroms_no_redudent.size(); i++)
	{
		ofile_record << chroms_no_redudent[i].length << "\t";
		for (int j = 0; j < city_num; j++)
			ofile_record << chroms_no_redudent[i].path[j] << "\t";
		ofile_record << endl;
	}
	ofile_record.close();
}
int GeneticAlgorithm::share_distance(const Chromosome & a, const Chromosome & b)
{
	return share_distance(a.path, b.path);
}

int GeneticAlgorithm::share_distance(const vector<int>& path1, const vector<int>& path2)
{
	int city_num = path1.size();
	vector<int> share_dist_vec(city_num, 0);
	int share_distance_count = 0;
	int p11, p12, p21, p22;
	for (int i = 0; i < city_num; i++)
	{
		p11 = i; p12 = (i + 1) % city_num;
		for (int j = 0; j < city_num; j++)
		{
			p21 = j; p22 = (j + 1) % city_num;

			if (path1[p11] == path2[p21] && path1[p12] == path2[p22] || path1[p11] == path2[p22] && path1[p12] == path2[p21])
			{
				share_distance_count++;
				break;
			}
		}

	}
	assert(share_distance_count <= city_num);
	return share_distance_count;
}
bool fitness_desc(const Chromosome& a, const Chromosome& b)
{

	return a.fitness > b.fitness; //asc

}
int find_place_index(vector<int> path, int num)
{
	for (int i = 0; i < path.size(); i++)
	{
		if (path[i] == num)
			return i;
	}
	return -1;
}
void GeneticAlgorithm::crossover(const string & crossover_strategy)
{
	int  p1, p2, p;
	for (int i = 0; i < NP; i += 2)
	{

		if (rnddouble(0.0, 1.0) >= crossover_ratio)
			continue;
		if (crossover_strategy == "CX")
		{//Cycle Crossover operator
			int p = rndint(0, city_num - 1);
			vector<bool> flag_parent1(city_num, false);
			vector<bool> flag_parent2(city_num, false);
			vector<int> path1(city_num, -1);
			vector<int> path2(city_num, -1);
			int inherite_index, find_place, find_index, count_tmp_found;
			//offspring 1
			count_tmp_found = 0;
			inherite_index = p;
			while (count_tmp_found < city_num)
			{
				for (int k = 0; k < city_num; k++)
				{
					if (flag_parent1[inherite_index] == true)
						inherite_index = (inherite_index + 1) % city_num;

					while (flag_parent1[inherite_index] != true)
					{

						path1[inherite_index] = chroms[i].path[inherite_index];
						flag_parent1[inherite_index] = true;
						count_tmp_found++;
						find_place = chroms[i + 1].path[inherite_index];
						inherite_index = find_place_index(chroms[i].path, find_place);
					}
				}
			}
			//offspring2
			count_tmp_found = 0;
			inherite_index = p;
			while (count_tmp_found < city_num)
			{
				for (int k = 0; k < city_num; k++)
				{
					if (flag_parent2[inherite_index] == true)
						inherite_index = (inherite_index + 1) % city_num;
				}
				while (flag_parent2[inherite_index] != true)
				{

					path2[inherite_index] = chroms[i + 1].path[inherite_index];
					flag_parent2[inherite_index] = true;
					count_tmp_found++;
					find_place = chroms[i].path[inherite_index];
					inherite_index = find_place_index(chroms[i + 1].path, find_place);
				}
			}
			chroms[i].path = path1;
			chroms[i + 1].path = path2;
		}
		if (crossover_strategy == "PMX")
		{//partially mapped crossover 
			vector<bool> exist_flag_p1(city_num, false);
			vector<bool> exist_flag_p2(city_num, false);
			vector<int> path1(city_num, -1);
			vector<int> path2(city_num, -1);
			vector<int> exist_city;

			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			if (p1 > p2)
				p2 += city_num;
			//offspring 1
			exist_city.clear();
			for (int kk = p1; kk < p2; kk++)
			{
				int k = kk % city_num;
				path1[k] = chroms[i].path[k];
				exist_city.push_back(path1[k]);
				exist_flag_p1[k] = true;
			}
			int inherite_index, inherite_city, find_place_p1, find_index, place_index, count_tmp_found;
			int wait_place_index_p2, wait_place_index_p1;
			for (int kk = p1; kk < p2; kk++)
			{
				int k = kk % city_num;
				find_index = find_place_index(exist_city, chroms[i + 1].path[k]);
				if (find_index != -1)
					continue;
				inherite_index = k;
				inherite_city = chroms[i + 1].path[inherite_index];
				place_index = inherite_index;
				do {
					find_place_p1 = chroms[i].path[place_index];
					place_index = find_place_index(chroms[i + 1].path, find_place_p1);
				} while (exist_flag_p1[place_index] == true);

				path1[place_index] = inherite_city;
				exist_flag_p1[place_index] = true;
			}
			for (int k = 0; k < city_num; k++)
			{
				if (exist_flag_p1[k] == false)
				{
					exist_flag_p1[k] = true;
					path1[k] = chroms[i + 1].path[k];
				}
			}
			// offspring2
			exist_city.clear();
			for (int kk = p1; kk < p2; kk++)
			{
				int k = kk % city_num;
				path2[k] = chroms[i + 1].path[k];
				exist_city.push_back(path2[k]);
				exist_flag_p2[k] = true;
			}
			for (int kk = p1; kk < p2; kk++)
			{
				int k = kk % city_num;
				find_index = find_place_index(exist_city, chroms[i].path[k]);
				if (find_index != -1)
					continue;
				inherite_index = k;
				inherite_city = chroms[i].path[inherite_index];
				place_index = inherite_index;
				do {
					find_place_p1 = chroms[i + 1].path[place_index];
					place_index = find_place_index(chroms[i].path, find_place_p1);
				} while (exist_flag_p2[place_index] == true);

				path2[place_index] = inherite_city;
				exist_flag_p2[place_index] = true;
			}
			for (int k = 0; k < city_num; k++)
			{
				if (exist_flag_p2[k] == false)
				{
					exist_flag_p2[k] = true;
					path2[k] = chroms[i].path[k];
				}
			}
			chroms[i].path = path1;
			chroms[i + 1].path = path2;
		}
	}
}
void GeneticAlgorithm::mutation(const string & mutation_strategy)
{
	int p1, p2;

	for (int i = 0; i < NP; i++)
	{
		if (rnddouble(0.0, 1.0) >= mutation_ratio)
			continue;
		if (mutation_strategy == "RSM")
		{//reverse the sequence
			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			if (p1 > p2)
			{
				int tmp = p1;
				p1 = p2;
				p2 = tmp;
			}
			reverse(chroms[i].path.begin() + p1, chroms[i].path.begin() + p2);
		}
		if (mutation_strategy == "PSM")
		{// point PSM
			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			int city_tmp = chroms[i].path[p1];
			chroms[i].path[p1] = chroms[i].path[p2];
			chroms[i].path[p2] = city_tmp;
		}
	}
}

bool sort_index_len_asc(const vector<int> &a, const vector<int> &b)
{
	return a[1] < b[1];
}

void GeneticAlgorithm::pass_next_generation_speciation_V2(GeneticAlgorithm &Offspring, const int &species_size)
{
	int share_distance_tmp = -1;
	int species_size_tmp = species_size;
	for (int s = 0; s < NP; s += species_size_tmp)
	{
		if (species_size_tmp + s >= NP)
			species_size_tmp = NP - s;

		int current_seed = s;
		int count = 0;
		for (int i = s; i < s + species_size_tmp; i++)
		{
			count = 0;
			int max_share_distance_index = -1;
			int max_share_distance = -1;
			for (int j = s; j < s + species_size_tmp; j++)
			{
				share_distance_tmp = share_distance(Offspring.chroms[i], chroms[j]);
				if (share_distance_tmp > max_share_distance)
				{
					max_share_distance_index = j;
					max_share_distance = share_distance_tmp;
				}
			}
			if (Offspring.chroms[i].length < chroms[max_share_distance_index].length)
			{
				chroms[max_share_distance_index] = Offspring.chroms[i];
			}
		}
	}
}

void GeneticAlgorithm::length_evaluate()
{
	for (int i = 0; i < chroms.size(); i++)
		chroms[i].evaluate();
	iter++;
}

void GeneticAlgorithm::speciation(const GeneticAlgorithm & Parent, const int species_size)
{

	vector<int> path_tmp(city_num);
	for (int j = 0; j < city_num; j++)
	{
		path_tmp[j] = j;
	}
	vector<index_len> chroms_index(NP);
	vector<int> shuffle_index(species_size - 1);
	for (int i = 0; i < species_size - 1; i++)
	{
		shuffle_index[i] = i;
	}
	vector<int> tmp(city_num, 0);
	for (int i = 0; i < NP; i++)
	{
		chroms_index[i].index = i;
		chroms_index[i].len = Parent.chroms[i].length;
	}
	sort(chroms_index.begin(), chroms_index.end(), sort_index_asc_second);
	vector<bool> flag_selected(NP, false);
	vector<index_len> chroms_index_share_distance(NP);
	int tmp_species_size = species_size;
	for (int i = 0; i < NP; i += tmp_species_size)
	{
		if (tmp_species_size + i >= NP)
			tmp_species_size = NP - i;

		int seed_index = 0;
		while (flag_selected[chroms_index[seed_index].index] == true)
		{
			seed_index++;
		}
		flag_selected[chroms_index[seed_index].index] = true;
		chroms[i] = Parent.chroms[chroms_index[seed_index].index];
		for (int j = 0; j < NP; j++)
		{
			int current_index = chroms_index[j].index;
			chroms_index_share_distance[j].index = current_index;
			if (flag_selected[current_index] != true)
				chroms_index_share_distance[j].len = city_num - share_distance(chroms[i], Parent.chroms[current_index]);
			else
				chroms_index_share_distance[j].len = 999999;
		}
		sort(chroms_index_share_distance.begin(), chroms_index_share_distance.end(), sort_index_asc_second);


		random_shuffle(shuffle_index.begin(), shuffle_index.end());
		int count_share = 0;
		for (int j = 0; j < tmp_species_size - 1; j++)
		{
			chroms[i + j + 1] = Parent.chroms[chroms_index_share_distance[shuffle_index[j]].index];
			flag_selected[chroms_index_share_distance[shuffle_index[j]].index] = true;
			count_share += chroms_index_share_distance[shuffle_index[j]].len;
		}
		if (count_share == 0)
		{
			random_shuffle(path_tmp.begin(), path_tmp.end());
			for (int j = 0; j < tmp_species_size - 1; j++)
			{
				chroms[i + j + 1] = path_tmp;
				chroms[i + j + 1].evaluate();
			}
		}
	}
}

GeneticAlgorithm::GeneticAlgorithm(int np, int cn, double cr, double mr, int MAX_ITER, City c, long long maxfes)
{
	NP = np;
	current_NP = 0;
	city_num = cn;
	mutation_ratio = mr;
	crossover_ratio = cr;
	max_iterations = MAX_ITER;
	iter = 0;
	city_list = c;
	FES = 0;
	maxFES = maxfes;
}

GeneticAlgorithm::GeneticAlgorithm(const GeneticAlgorithm & b)
{
	NP = b.NP;
	current_NP = b.current_NP;
	chroms = b.chroms;
	city_num = b.city_num;
	mutation_ratio = b.mutation_ratio;
	crossover_ratio = b.crossover_ratio;
	max_iterations = b.max_iterations;
	iter = b.iter;
	city_list = b.city_list;
	CF = b.CF;
	maxFES = b.maxFES;
}
void GeneticAlgorithm::initiate()
{
	int ratio;
	vector<int> path_tmp(city_num);
	for (int i = 0; i < city_num; i++)
	{
		path_tmp[i] = i;
	}
	int next_city;
	for (int i = 0; i < NP; i++)
	{
		random_shuffle(path_tmp.begin(), path_tmp.end());
		chroms.push_back(Chromosome(path_tmp));
	}
	vector<int> share_distance_vec(city_num, 0);
}

void GeneticAlgorithm::crossover()
{
	int p1, p2;
	int count_p, selected_city;
	Chromosome parent[2];

	for (int i = 0; i < NP; i += 2)
	{
		if (rnddouble(0, 1) <= crossover_ratio)
		{
			vector<bool> exist_flag_p1(city_num, false);
			vector<bool> exist_flag_p2(city_num, false);
			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			if (p1 > p2)
			{
				int tmp = p1;
				p1 = p2;
				p2 = tmp;
			}
			parent[0] = chroms[i];
			parent[1] = chroms[i + 1];

			for (int j = p1; j < p2; j++)
			{
				chroms[i].path[j] = parent[1].path[j];  exist_flag_p1[chroms[i].path[j]] = true;
				chroms[i + 1].path[j] = parent[0].path[j];  exist_flag_p2[chroms[i + 1].path[j]] = true;
			}
			count_p = p2;
			for (int j = 0; j < city_num - (p2 - p1); j++)
			{
				do {
					selected_city = parent[0].path[count_p];
					count_p = (count_p + 1) % city_num;
				} while (exist_flag_p1[selected_city]);
				chroms[i].path[(p2 + j) % city_num] = selected_city;
				exist_flag_p1[selected_city] = true;
			}
			count_p = p2;
			for (int j = 0; j < city_num - (p2 - p1); j++)
			{
				do {
					selected_city = parent[1].path[count_p];
					count_p = (count_p + 1) % city_num;
				} while (exist_flag_p2[selected_city]);
				chroms[i + 1].path[(p2 + j) % city_num] = selected_city;
				exist_flag_p1[selected_city] = true;
			}
		}
	}
}

void GeneticAlgorithm::mutation()
{
	int p1, p2;
	for (int i = 0; i < NP; i++)
	{
		if (rnddouble(0, 1) <= mutation_ratio)
		{
			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			if (p1 > p2)
			{
				int tmp = p1;
				p1 = p2;
				p2 = tmp;
			}
			reverse(chroms[i].path.begin() + p1, chroms[i].path.begin() + p2);
		}
	}
}
int share_distance(const vector<int> & path, const Chromosome & b)
{
	return share_distance(path, b.path);
}

class Problem
{
public:
	City cities;
	vector<vector<int> > solutions;
	string file_string;
	string problem_name;
	int min_len;
public:
	Problem(City &c, vector<vector<int> > &s, const int &l, string & n);
	Problem(City &c, const string &problem_name);
	friend class GeneticAlgorithm;
};
Problem::Problem(City &c, vector<vector<int> > &s, const int &l, string & n) {
	cities = c;
	solutions = s;
	MIN_LEN = l;
	min_len = l;
	problem_name = n;
};
Problem::Problem(City &city, const string &problem_name_n)
{
	cities = city;
	problem_name = problem_name_n;
}

void convert_file_name(const string & file_name)
{
	ostringstream ostr;
	ostr << file_name;
	string str = ostr.str();
	int m = 0;
	for (m = 0; m <= str.length(); m++)
		file_name_char[m] = str[m];
	file_name_char[m] = '\0';
}
int main(int argc, char* argv[])
{
	int Func_index[] = { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39, 40, 41 };
	int MaxFES[] = { 60000,60000,60000,60000,60000,60000,60000,60000,60000,60000,60000,60000, // 1- 12
		1200000,1200000,1200000,1200000,1200000,1200000,1200000 , //13-19
		1200000, 1200000, 1200000, 1200000, 1200000, 1200000, 0, 0, 0, 0,0,	 // 20-30
		0, 0, 0, 0, 0, 0, 0	,0, 0, 0, 0				// 31-38
	};

	string Func_name[] = {
		"simple1_9", "simple2_10", "simple3_10", "simple4_11", "simple5_12", "simple6_12",			//1-6
		"geometry1_10", "geometry2_12", "geometry3_10", "geometry4_10", "geometry5_10", "geometry6_15",	//7-12
		"composite1_28", "composite2_34", "composite3_22", "composite4_33", "composite5_35",  "composite7_42", "composite8_60" ,	 //13-19
		"composite6_39", "composite8_45", "composite9_48", "composite10_55", "composite11_59","composite13_66", "", "", "", "",	"",//20 - 25
		"att48","eil51","berlin52","st70","eil76", "pr76", "eil101","kroA100", "lin105", "kroA150", "lin318"	//31-41
	};
	input_dir = "./benchmark_MSTSP/";

	int FUNC_BEGIN = atoi(argv[1]);
	int FUNC_END = atoi(argv[2]);
	int RUNS = atoi(argv[3]);
	int MAX_ITER = -1;
	int NP = atoi(argv[4]);
	double  crossover_ratio = atof(argv[5]);
	double mutation_ratio = atof(argv[6]);
	output_dir = "./" + string(argv[7]) + "/";
	string selection_operator = string(argv[8]);
	string crossover_operator = string(argv[9]);
	string mutation_operator = string(argv[10]);
	int species_size = atoi(argv[11]);
	int track_gap;
	string file_string = output_dir + "info.param";
	ofstream ofile_record;
	convert_file_name(file_string);
	ofile_record.open(file_name_char, ios::out);
	ofile_record << "NP: " << NP << endl
		<< "maxIteration: " << MAX_ITER << endl
		<< "RUNS: " << RUNS << endl
		<< "crossove Rate: " << crossover_ratio << endl
		<< "mutation Rate: " << mutation_ratio << endl
		<< "species size: " << species_size << endl;
	ofile_record.close();

	srand(time(NULL));
	string problem_name;
	double refine_niche_radius;
	for (int func = FUNC_BEGIN; func <= FUNC_END; func++)
	{
		problem_name = Func_name[func];
		if (problem_name == "")
			continue;
		City cities(problem_name);
		Problem TSP_pro(cities, problem_name);

		if (func >= 0 && func < 29)
		{
			refine_niche_radius = 0.01;
			MAX_ITER = MaxFES[func] / NP;
		}
		else	if (func >= 30)
		{
			MaxFES[func] = cities.city_num * 50000;
			MAX_ITER = MaxFES[func] / NP;
			refine_niche_radius = 0.05;
		}

		string fs_file_time;
		fstream fs_time;
		clock_t  clockBegin, clockEnd, clockEnd_gap;
		string file_string_time = output_dir + problem_name + ".alg_time";
		ofstream ofile_record_time;
		convert_file_name(file_string_time);
		ofile_record_time.open(file_name_char, ios::out);
		if (ofile_record_time.is_open() == false)
			return -1;

		for (int run = 0; run < RUNS; run++)
		{
			cout << problem_name << ": Begin " << run << endl;
			int track_count = 1;
			clockBegin = clock();
			int time_gap = 0;

			GeneticAlgorithm Parent(NP, cities.city_num, crossover_ratio, mutation_ratio, MAX_ITER, cities, MaxFES[func]);
			Parent.initiate();
			Parent.length_evaluate();
			GeneticAlgorithm Offspring = Parent;
			while (!Offspring.terminate())
			{
				Offspring.speciation(Parent, species_size);
				Offspring.crossover(crossover_operator);
				Offspring.mutation(mutation_operator);
				Offspring.length_evaluate();
				Parent.pass_next_generation_speciation_V2(Offspring, species_size);
			}
			Offspring.remove_redudant_report2_file_V3(Parent, 0.01, run, MaxFES[func], problem_name, species_size);
			clockEnd = clock();
			ofile_record_time << (double)(clockEnd - clockBegin - time_gap) / CLOCKS_PER_SEC << endl;
			cout << problem_name << ": Finished " << run << endl;
		} // for (int run = 0; run < RUNS; run++)
		ofile_record_time.close();
	}
#ifdef _MSC_VER
	system("pause");
#endif
	return 0;
}