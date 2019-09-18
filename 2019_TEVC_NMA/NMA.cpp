/**
*  Copyright (c) 2018 SCUT
*  All rights reserved.
*
*  @author      Ting Huang (cshting@mail.scut.edu.cn)
*  @date		2019-01-29
*
*  @brief       The niching memetic algorithm (NMA) is proposed to solve MSTSPs.
*
*  @version
*	- V6: minor modifications
*
*  @command
*		[begin Index] [end index] [runs] [NP] [crossover rate] [mutation rate] [dir] [selection strategy] [crossover strategey] [mutation strategy]
*		 0 24 50 150 0.9 0.01 .\Results roulette critical_PMX NicheEM 2-opt
*
*	@reference
*		T. Huang, Y. Gong, S. Kwong, H. Wang and J. Zhang, “A Niching Memetic Algorithm for Multi-Solution Traveling Salesman Problem,” IEEE Transactions on Evolutionary Computation. DOI: 10.1109/TEVC.2019.2936440.
*/

typedef int distType;

#include <iostream>
#include <vector>
#include <ctime>
#include <cmath>
#include <limits>
#include <sstream>  
#include <cstring>  
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <cassert>
#include <numeric>
#include <list>
#include <typeinfo>

#undef _MSC_VER

using namespace std;
vector < vector<distType> > dist;
vector < vector<distType> > neighbor_matrix;
vector<vector<distType> > disimilar_dist_matrix;
vector<vector<int> > share_dist_matrix;
vector<vector<int> > share_dist_matrix2;
vector<vector<int> > share_matrix;
char file_name_char[128];
void convert_file_name(const string &);
string input_dir;
string output_dir;
distType MIN_LEN;
long long FES = 0;
int SPECIES_SIZE = 0;
class GeneticAlgorithm;
class Problem;
class FIFO_list;
class Chromosome;
class index_value_V2;
const int MAXN = 80;
int STAGNT_NUM = 50;
double CURR_SP_EDGE_MATRIX[MAXN][MAXN];			// record the critical edge of current species
double INCRE_GEN_EDGE_MATRIX[MAXN][MAXN];		//record the edge set of sucessive generations
bool FLAG_VISITED[MAXN];
distType BESTLEN = -1;
distType PRE_BESTLEN = -1;

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
class index_len
{
public:
	int index;
	distType len;
};
bool sort_index_asc_second(const index_len& a, const index_len & b)
{
	return a.len < b.len;
}
bool sort_index_des_second(const index_len& a, const index_len & b)
{
	return a.len > b.len;
}


class index_value_V2
{
public:
	int index;
	double value;
};
bool sort_index_value_V2_asc(const index_value_V2 &a, const index_value_V2 &b)
{
	return a.value < b.value;
}
bool sort_index_value_V2_des(const index_value_V2 &a, const index_value_V2 &b)
{
	return a.value > b.value;
}

class City
{
public:
	vector<cardinate> city_vec;
	int city_num;
public:
	distType distance(cardinate, cardinate);
	City(const string &problem_name);
	City() :city_num(0) {};
};
City::City(const string &problem_name) {
	string file_string = input_dir + problem_name + ".tsp";
	ifstream ifile_tsp;
	convert_file_name(file_string);
	ifile_tsp.open(file_name_char, ios::in);
	if (!ifile_tsp)
	{
		cout << file_name_char << " not exist" << endl;
		system("pause");
		exit(-1);
	}
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

	for (int i = 0; i < dist.size(); i++)
	{
		dist[i].clear();
		neighbor_matrix[i].clear();
	}
	dist.clear();
	neighbor_matrix.clear();

	vector<distType> dist_vec(city_num, 0);
	for (int i = 0; i < city_num; i++)
		dist.push_back(dist_vec);
	for (int i = 0; i < city_num; i++)
	{
		for (int j = 0; j < city_num; j++)
			dist[i][j] = dist[j][i] = distance(city_vec[i], city_vec[j]);
	}
	vector<index_len> distlen_vec;
	index_len distlen_tmp;
	for (int i = 0; i < city_num; i++)
	{
		distlen_vec.clear();
		for (int j = 0; j < city_num; j++)
		{
			distlen_tmp.index = j;
			distlen_tmp.len = dist[i][j];
			distlen_vec.push_back(distlen_tmp);
		}
		sort(distlen_vec.begin(), distlen_vec.end(), sort_index_asc_second);
		for (int j = 0; j < city_num; j++)
		{
			dist_vec[j] = distlen_vec[j].index;
		}
		neighbor_matrix.push_back(dist_vec);
	}
}


distType City::distance(cardinate a, cardinate b)
{
	return  (distType)(sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y)) + 0.5);
}
class FIFO_item
{
public:
	int index1;
	int index2;
	int time;
	int successive_generation;
};
class FIFO_list
{
public:
	vector <FIFO_item> ciritcal_edge_vec;

	bool update_existed_critical_edge(const FIFO_item &FIFO_item);
	bool update_existed_critical_index(const int &index);
	void append_item(const FIFO_item &FIFO_item);
	void report();
	void final_report();
	void update_successive_generation(const int &FES);
};
void FIFO_list::update_successive_generation(const int & FES)
{
	vector <FIFO_item>::iterator erase_iter;
	for (vector <FIFO_item>::iterator iter = ciritcal_edge_vec.begin(); iter != ciritcal_edge_vec.end(); )
	{
		if ((*iter).time == FES)
		{
			(*iter).successive_generation = 0;
		}
		else
		{
			(*iter).successive_generation++;
		}

		if ((*iter).successive_generation >= STAGNT_NUM)
		{

			INCRE_GEN_EDGE_MATRIX[(*iter).index1][(*iter).index2] = 0;
			iter = ciritcal_edge_vec.erase(iter);
		}
		else
		{
			iter++;
		}
	}

}
void FIFO_list::append_item(const FIFO_item & FIFO_item)
{
	ciritcal_edge_vec.push_back(FIFO_item);
}
void FIFO_list::report()
{
	cout << "FIFO_list " << ciritcal_edge_vec.size() << endl;
	for (vector <FIFO_item>::iterator iter = ciritcal_edge_vec.begin(); iter != ciritcal_edge_vec.end(); iter++)
	{
		cout << (*iter).index1 << "-" << (*iter).index2 << ":" << (*iter).time << ":" << (*iter).successive_generation << endl;
	}
	cout << endl;
}
void FIFO_list::final_report()
{

	for (vector <FIFO_item>::iterator iter = ciritcal_edge_vec.begin(); iter != ciritcal_edge_vec.end(); iter++)
	{
		cout << (*iter).index1 << "-" << (*iter).index2 << ":" << (*iter).successive_generation << "/" << INCRE_GEN_EDGE_MATRIX[(*iter).index1][(*iter).index2] << endl;
	}

	if (!ciritcal_edge_vec.empty())
		cout << endl << "****************** critical edge *****************" << endl;
}
bool FIFO_list::update_existed_critical_edge(const FIFO_item & in_FIFO_item)
{
	for (vector <FIFO_item>::iterator iter = ciritcal_edge_vec.begin(); iter != ciritcal_edge_vec.end(); iter++)
	{
		if (in_FIFO_item.index1 == (*iter).index1 && in_FIFO_item.index2 == (*iter).index2)
		{
			(*iter).time = in_FIFO_item.time;
			return true;
		}

	}
	return false;
}


bool FIFO_list::update_existed_critical_index(const int &index)
{
	int oldest_time = -1;
	int oldest_index = -1;
	int existed_count = 0;
	vector <FIFO_item>::iterator erase_iter;
	int i = 0;
	for (vector <FIFO_item>::iterator iter = ciritcal_edge_vec.begin(); iter != ciritcal_edge_vec.end(); iter++)
	{
		if ((*iter).index1 == index || (*iter).index2 == index)
		{
			existed_count++;
			if (oldest_time == -1 || (*iter).time < oldest_time)
			{
				oldest_time = (*iter).time;
				oldest_index = i;
				erase_iter = iter;
			}
		}
		i++;
	}
	if (existed_count >= 2)
	{
		assert(INCRE_GEN_EDGE_MATRIX[(*erase_iter).index1][(*erase_iter).index2] > 0);
		INCRE_GEN_EDGE_MATRIX[(*erase_iter).index1][(*erase_iter).index2] = 0;
		ciritcal_edge_vec.erase(erase_iter);
	}
	else
	{
		return false;
	}
}
class Chromosome
{
public:
	vector<int> path;
	distType length;
	double fitness;
	distType pre_length;
	int successive_generation;
public:
	Chromosome() :length(0) {};
	Chromosome(vector<int> p);
	Chromosome(const Chromosome & b);
	void evaluate();
	distType distance(int a, int b) { return  (distType)dist[a][b]; };
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

	if (pre_length == 0)
		pre_length = length;

}
Chromosome::Chromosome(vector<int> p)
{
	path = p;
	length = 0;
	fitness = 0.0;
	pre_length = 0.0;
	successive_generation = 0;
}
Chromosome::Chromosome(const Chromosome & b)
{
	path = b.path;
	length = b.length;
	fitness = b.fitness;
	successive_generation = b.successive_generation;
	pre_length = b.pre_length;
}
class Species_info
{
public:
	int begin_index;
	int end_index;
	int species_size;
};
class GeneticAlgorithm
{
private:
	vector<Chromosome> chroms;
	vector<Species_info> sp_info;

	distType min_len, max_len;
	int max_critical_threshold_gen;

	vector<vector<distType> > neighbor_maxti;
	City city_list;
	int NP;
	int current_NP;
	int city_num;
	double init_ratio;
	int specie_size;

	long long maxFES;
	double mutation_ratio;
	double crossover_ratio;
	int max_iterations;
	string crossover_opertor;
	string origin_crossover_operator;
	string mutation_operator;
	string local_search_operator;
public:

	int iter;
	friend class Problem;
	GeneticAlgorithm(const GeneticAlgorithm & b);
	GeneticAlgorithm(int np, int cn, double cr, double mr, int MAX_ITER, City c, long long maxfes,
		int sps_s, double init_r, const string &cr_op, const string &mt_op, const string &lc_op);

	void initiate_V2(Chromosome &);
	void initiate_V2();

	void speciation_V3(const GeneticAlgorithm &, const int, const int);
	void speciation_V3_TSPLIB(const GeneticAlgorithm &);


	void selection(GeneticAlgorithm &, const string &);
	int localSearch(Chromosome &);
	void localSearch();
	void localSearch_V3_Critical();
	int localSearch_V3_Critical_operation(Chromosome &);

	void crossover();
	void mutation(Chromosome &chrom);
	void mutation();

	void pass_next_generation_speciation_V3(GeneticAlgorithm &, const int &);

	void report_population(int num, const string &);
	void length_evaluate();
	void remove_redudant_report2_file_V6(double threshold, const int &run, const int &track_fes, const string &problem_name, const int species_size);
	bool sort_cmp_len(const vector<Chromosome> &a, const vector<Chromosome> &b) {
		return a[0].length < b[0].length;
	}
	bool terminate() { return FES >= maxFES; };
	int share_distance(const vector<int>& a, const vector<int>& b);
	int share_distance(const Chromosome& a, const Chromosome& b);
	void collect_run_distance(const int &begin, const int &species_size);
	void report();

	void record_edge_matrix();
	void update_edge_matrix();
	void report_critical(const int &FES);

	void empty_edge_matrix();
	bool exist_critical_edge(const int &p1, const int &p2, const Chromosome &chrom);
	bool exist_critical_edge(const int &p1, const int &p2);
	void add_value_matrix(double matrix[][MAXN], const int & index1, const int & index2, const double & value);

	FIFO_list critical_list;
	void update_sucessive();
	void report_sucessive(const int num);

};

void GeneticAlgorithm::report_sucessive(const int num)
{
	cout << "/////////////////////////" << endl;
	for (int i = 0; i < num; i++)
	{
		cout << chroms[i].pre_length << "/" << chroms[i].length << ":" << chroms[i].successive_generation << "\t";
	}
	cout << endl;
}

void GeneticAlgorithm::update_sucessive()
{
	for (int i = 0; i < NP; i++)
	{
		if (chroms[i].length == chroms[i].pre_length)
		{
			chroms[i].successive_generation++;
		}
		else
		{
			chroms[i].successive_generation = 0;
		}
		chroms[i].pre_length = chroms[i].length;
	}
}
void GeneticAlgorithm::add_value_matrix(double matrix[][MAXN], const int & index1, const int & index2, const double & value)
{
	if (index1 < index2)

		matrix[index1][index2] += value;
	else
		matrix[index2][index1] += value;
}
void GeneticAlgorithm::record_edge_matrix()
{
	for (int i = 0; i < city_num; i++)
	{
		memset(CURR_SP_EDGE_MATRIX[i], 0, sizeof(double)*MAXN);
	}

	distType bestLen = -1;
	vector <int> best_index_list;

	for (int sp = 0; sp < sp_info.size(); sp++)
	{
		int	i = sp_info[sp].begin_index;
		if (chroms[sp_info[sp].begin_index].successive_generation >= STAGNT_NUM)
		{
			bestLen = chroms[i].length;
			bool flag_exist = false;
			for (int k = 0; k < best_index_list.size(); k++)
			{
				double tmp_sharing_distance = share_distance(chroms[i], chroms[best_index_list[k]]);
				if (tmp_sharing_distance == city_num)
				{
					flag_exist = true;
					break;
				}
			}
			if (!flag_exist)
				best_index_list.push_back(i);
		}
	}
	if (best_index_list.size() <= 1)
		return;
	for (vector<int>::iterator iter = best_index_list.begin(); iter < best_index_list.end(); iter++)
	{
		int i = *iter;
		for (int j = 0; j < city_num - 1; j++)
		{
			add_value_matrix(CURR_SP_EDGE_MATRIX, chroms[i].path[j], chroms[i].path[j + 1], 1);
		}
		add_value_matrix(CURR_SP_EDGE_MATRIX, chroms[i].path[city_num - 1], chroms[i].path[0], 1);
	}
	for (int i = 0; i < city_num; i++)
	{
		for (int j = i + 1; j < city_num; j++)
		{
			CURR_SP_EDGE_MATRIX[i][j] /= (double)best_index_list.size();
		}
	}
}
bool GeneticAlgorithm::exist_critical_edge(const int & p1, const int & p2, const Chromosome & chrom)
{
	int city1, city2;
	if (chrom.path[p1] < chrom.path[p2])
	{
		city1 = chrom.path[p1];
		city2 = chrom.path[p2];
	}
	else
	{
		city2 = chrom.path[p1];
		city1 = chrom.path[p2];
	}
	if (INCRE_GEN_EDGE_MATRIX[city1][city2] >= 1)
		return true;
	else
		return false;
}

bool GeneticAlgorithm::exist_critical_edge(const int & city1, const int & city2)
{
	if (city1 < city2 && INCRE_GEN_EDGE_MATRIX[city1][city2] >= 1 ||
		city1 > city2 && INCRE_GEN_EDGE_MATRIX[city2][city1] >= 1)
		return true;
	else
		return false;
}

void GeneticAlgorithm::empty_edge_matrix()
{
	for (int i = 0; i < city_num; i++)
	{
		memset(CURR_SP_EDGE_MATRIX[i], 0, sizeof(double)*MAXN);
		memset(INCRE_GEN_EDGE_MATRIX[i], 0, sizeof(double) * MAXN);
	}
}
void GeneticAlgorithm::report_critical(const int &XFES)
{
	bool flag_exist = false;
	for (int i = 0; i < city_num; i++)
	{
		for (int j = i + 1; j < city_num; j++)
		{
			if (exist_critical_edge(i, j))
			{
				cout << i << "-" << j << ":" << INCRE_GEN_EDGE_MATRIX[i][j] << endl;
				flag_exist = true;
			}
		}
	}
	if (flag_exist)
		cout << "------------------" << FES << "/" << iter << "------------------" << endl;
	critical_list.final_report();
}

void GeneticAlgorithm::update_edge_matrix()
{
	FIFO_item tmp_FIFO_item;
	bool existed_item;
	for (int i = 0; i < city_num; i++)
	{
		for (int j = i + 1; j < city_num; j++)
		{
			if (CURR_SP_EDGE_MATRIX[i][j] == 1)
			{
				tmp_FIFO_item.index1 = i;
				tmp_FIFO_item.index2 = j;
				tmp_FIFO_item.time = FES;
				tmp_FIFO_item.successive_generation = 0;
				if (INCRE_GEN_EDGE_MATRIX[i][j] != 0)
				{
					existed_item = critical_list.update_existed_critical_edge(tmp_FIFO_item);
				}
				else
				{
					critical_list.update_existed_critical_index(i);
					critical_list.update_existed_critical_index(j);
					critical_list.append_item(tmp_FIFO_item);
				}
				INCRE_GEN_EDGE_MATRIX[i][j]++;
			}
		}
	}
	critical_list.update_successive_generation(FES);
}
void GeneticAlgorithm::collect_run_distance(const int & begin, const int & tmp_species_size)
{
	int species_size = tmp_species_size;
	if (begin + species_size > NP)
		species_size = NP - begin;
	double accu_share_dist, mean_share_dist, var_share_dist;
	accu_share_dist = mean_share_dist = var_share_dist = 0.0;
	vector<double> share_dist;
	int bestIndex = -1;
	for (int i = begin; i < begin + species_size; i++)
	{
		if (bestIndex == -1 || chroms[bestIndex].length > chroms[i].length)
			bestIndex = i;
	}

	for (int i = begin; i < begin + species_size; i++)
	{
		share_dist.push_back(city_num - share_distance(chroms[bestIndex], chroms[i]));
	}
	accu_share_dist = accumulate(share_dist.begin(), share_dist.end(), 0.0);
	mean_share_dist = accu_share_dist / (double)share_dist.size();

	accu_share_dist = 0.0;
	std::for_each(std::begin(share_dist), std::end(share_dist), [&](const double d) {
		accu_share_dist += (d - mean_share_dist)*(d - mean_share_dist);	});

	var_share_dist = sqrt(accu_share_dist / (share_dist.size() - 1));
	cout << "[best]: " << chroms[bestIndex].length << "\t" << "[mean]:" << mean_share_dist << "\t [var]:" << var_share_dist << endl;
}


int GeneticAlgorithm::localSearch(Chromosome &chrom)
{
	if (local_search_operator == "2-opt")
	{
		int count = 0;

		distType pre_dist = chrom.length;
		distType best_improve = 0;
		distType prove;
		int best_i, best_j;
		distType eliminate_len, add_len;
		for (int i = 0; i < city_num - 2; i++)
		{
			for (int j = i + 1; j < city_num - 1; j++)
			{
				prove = dist[chrom.path[(i - 1 + city_num) % city_num]][chrom.path[i]] + dist[chrom.path[j]][chrom.path[(j + 1) % city_num]]
					- dist[chrom.path[(i - 1 + city_num) % city_num]][chrom.path[j]] - dist[chrom.path[i]][chrom.path[(j + 1) % city_num]]; // improveÔ½´óÔ½ºÃ
				count += 4;
				if (prove >= best_improve)
				{
					best_improve = prove;
					best_i = i;
					best_j = j;
				}
			}
		}
		if (best_improve != 0)
		{
			reverse(chrom.path.begin() + best_i, chrom.path.begin() + best_j + 1);
			chrom.length -= best_improve;
		}
		return count;
	}
}
int GeneticAlgorithm::localSearch_V3_Critical_operation(Chromosome &chrom)
{
	int count = 0;

	distType pre_dist = chrom.length;
	distType best_improve = 0;
	distType prove;
	int best_i, best_j;
	distType eliminate_len, add_len;
	for (int i = 0; i < city_num - 2; i++)
	{
		int pre_i = (i - 1 + city_num) % city_num;
		if (exist_critical_edge(pre_i, i, chrom))
			continue;
		for (int j = i + 1; j < city_num - 1; j++)
		{
			int next_j = (j + 1) % city_num;
			if (exist_critical_edge(j, next_j, chrom))
				continue;
			prove = dist[chrom.path[(i - 1 + city_num) % city_num]][chrom.path[i]] + dist[chrom.path[j]][chrom.path[(j + 1) % city_num]]
				- dist[chrom.path[(i - 1 + city_num) % city_num]][chrom.path[j]] - dist[chrom.path[i]][chrom.path[(j + 1) % city_num]]; // improveÔ½´óÔ½ºÃ
			count += 4;
			if (prove >= best_improve)
			{
				best_improve = prove;
				best_i = i;
				best_j = j;
			}

		}
	}
	if (best_improve != 0)
	{
		reverse(chrom.path.begin() + best_i, chrom.path.begin() + best_j + 1);
		chrom.length -= best_improve;
	}
	return count;

}
void GeneticAlgorithm::localSearch_V3_Critical()
{
	vector<index_len> index_len_vec;
	index_len index_len_tmp;
	int count = 0;
	for (int i = 0; i < sp_info.size(); i++)
	{
		index_len_vec.clear();
		for (int j = sp_info[i].begin_index; j < sp_info[i].end_index; j++)
		{
			index_len_tmp.index = j;
			index_len_tmp.len = chroms[j].length;
			index_len_vec.push_back(index_len_tmp);
		}
		sort(index_len_vec.begin(), index_len_vec.end(), sort_index_asc_second);
		if (sp_info[i].species_size == 4)
		{
			count += localSearch_V3_Critical_operation(chroms[index_len_vec[0].index]);
		}
		else
		{
			for (int j = 0; j < index_len_vec.size() / 2; j++)
			{
				if (j == 0 || index_len_vec[j].len != index_len_vec[j - 1].len)
					count += localSearch_V3_Critical_operation(chroms[index_len_vec[j].index]);
			}
		}

	}
	FES += count / city_num;
}

void GeneticAlgorithm::localSearch()
{
	int count = 0;
	for (int i = 0; i < NP; i++)
	{
		count += localSearch(chroms[i]);
	}
	FES += (count / city_num);
}
class index_value
{
public:
	int index;
	distType value;
};
bool sort_index_asc_value(const index_value &a, const index_value &b)
{
	return a.value < b.value;
}
bool sort_index_des_value(const index_value &a, const index_value &b)
{
	return a.value > b.value;
}



void GeneticAlgorithm::remove_redudant_report2_file_V6(double threshold, const int & run, const int &track_fes, const string & problem_name, int species_size)
{
	vector<Chromosome> chroms_no_redudent;
	vector<Chromosome> chroms_tmp = chroms;
	vector<index_value> index_tmp_vec;
	vector<index_value> index_tmp_vec_all_species_best;
	index_value index_tmpp;
	int species_size_tmp;
	int current_index = 0;

	for (int sp = 0; sp < sp_info.size(); sp++)
	{
		index_tmp_vec.clear();
		species_size_tmp = sp_info[sp].species_size;

		for (int j = sp_info[sp].begin_index; j <= sp_info[sp].end_index; j++)
		{
			index_tmpp.index = j;
			index_tmpp.value = chroms_tmp[j].length;
			index_tmp_vec.push_back(index_tmpp);
		}
		sort(index_tmp_vec.begin(), index_tmp_vec.end(), sort_index_asc_value);

		vector<Chromosome> chrom_vec;
		for (int j = 0; j < index_tmp_vec.size(); j++)
		{
			chrom_vec.push_back(chroms_tmp[index_tmp_vec[j].index]);
		}


		for (int j = sp_info[sp].begin_index, i = 0; j <= sp_info[sp].end_index; i++, j++)
		{
			chroms_tmp[j] = chrom_vec[i];
		}

		index_tmpp.index = sp;
		index_tmpp.value = chroms_tmp[sp_info[sp].begin_index].length;
		index_tmp_vec_all_species_best.push_back(index_tmpp);
	}
	index_tmp_vec.clear();
	sort(index_tmp_vec_all_species_best.begin(), index_tmp_vec_all_species_best.end(), sort_index_asc_value);


	distType min_len_tmp = index_tmp_vec_all_species_best[0].value;
	distType threshold_min_len = min_len_tmp * (1 + threshold);

	for (int sp = 0; sp < index_tmp_vec_all_species_best.size(); sp++)
	{
		species_size_tmp = sp_info[index_tmp_vec_all_species_best[sp].index].species_size;
		int begin_index = sp_info[index_tmp_vec_all_species_best[sp].index].begin_index;
		for (int j = begin_index; j < begin_index + species_size_tmp; j++)
		{
			bool exist_flag = false;
			distType max_share_dist = -999999;
			for (int k = 0; k < chroms_no_redudent.size(); k++)
			{
				int share_dist_tmp = share_distance(chroms_tmp[j], chroms_no_redudent[k]);
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

			if (chroms_tmp[j].length == min_len_tmp)
			{
				chroms_no_redudent.push_back(chroms_tmp[j]);
			}
			else if (j == begin_index && chroms_tmp[j].length < threshold_min_len)
			{

				if (max_share_dist < (floor)(0.8*city_num))
				{
					chroms_no_redudent.push_back(chroms_tmp[j]);
					break;
				}
			}
			else
				break;
		}


	}
	string file_string;
	if (track_fes != 0)
	{
		file_string = output_dir + problem_name + "_RUN_" + to_string(run) + "_FES_" + to_string(track_fes) + ".alg_solution";
	}
	else
	{
		file_string = output_dir + problem_name + "_RUN_" + to_string(run) + ".alg_solution";
	}

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

void GeneticAlgorithm::crossover()
{

	int  p1, p2, p;
	for (int i = 0; i < NP; i += 2)
	{

		if (rnddouble(0.0, 1.0) >= crossover_ratio)
			continue;
		if (crossover_opertor == "critical_PMX")	//Partially-mapped crossover 1985  
		{
			vector<bool> exist_flag_p1(city_num, false);
			vector<bool> exist_flag_p2(city_num, false);
			vector<int> p1_mapping_p2(city_num, -1);
			vector<int> p2_mapping_p1(city_num, -1);

			vector<int> path1(city_num, -1);
			vector<int> path2(city_num, -1);

			p1 = rndint(0, city_num - 1);
			do
			{
				p2 = rndint(0, city_num - 1);
			} while (p1 == p2);
			if (p1 > p2)
				p2 += city_num;
			//offspring1
			vector<int> *path;
			vector<bool> *flag;
			Chromosome *chrom1, *chrom2;
			vector<int> *mapping;
			for (int m = 1; m <= 2; m++)
			{
				if (m == 1)
				{
					path = &path1;
					flag = &exist_flag_p1;
					chrom1 = &chroms[i];
					chrom2 = &chroms[i + 1];
					mapping = &p2_mapping_p1;
				}
				else
				{
					path = &path2;
					flag = &exist_flag_p2;
					chrom1 = &chroms[i + 1];
					chrom2 = &chroms[i];
					mapping = &p1_mapping_p2;
				}
				for (int kk = p1; kk < p2; kk++)
				{
					int k = kk % city_num;
					(*path)[k] = chrom2->path[k];
					(*flag)[(*path)[k]] = true;
					(*mapping)[(*path)[k]] = chrom1->path[k];
				}
				for (int k = 0; k < city_num; k++)
				{
					int kk1 = k;
					int kk2 = (k + 1) % city_num;
					if ((*flag)[chrom2->path[kk1]] == false || (*flag)[chrom2->path[kk2]] == false)
					{
						if (exist_critical_edge(chrom2->path[kk1], chrom2->path[kk2]))

						{
							if ((*flag)[chrom2->path[kk1]] == false)
							{
								(*path)[kk1] = chrom2->path[kk1];
								(*flag)[(*path)[kk1]] = true;
								(*mapping)[(*path)[kk1]] = chrom1->path[kk1];
							}

							if ((*flag)[chrom2->path[kk2]] == false)
							{

								(*path)[kk2] = chrom2->path[kk2];
								(*flag)[(*path)[kk2]] = true;
								(*mapping)[(*path)[kk2]] = chrom1->path[kk2];
							}
						}
					}
				}
				for (int k = 0; k < city_num; k++)
				{
					if ((*path)[k] != -1)
						continue;
					if ((*flag)[chrom1->path[k]] == false)
					{
						(*path)[k] = chrom1->path[k];
					}
					else
					{
						int inherite_city = chrom1->path[k];
						while ((*flag)[inherite_city] == true)
						{
							inherite_city = (*mapping)[inherite_city];
						}
						(*path)[k] = inherite_city;
					}
					(*flag)[(*path)[k]] = true;

				}
			}
			chroms[i].path = path1;
			chroms[i + 1].path = path2;
		}
	}
}
void GeneticAlgorithm::mutation()
{
	for (int i = 0; i < NP; i++)
	{
		if (rnddouble(0.0, 1.0) <= mutation_ratio)
			mutation(chroms[i]);
	}
}



void GeneticAlgorithm::mutation(Chromosome &chrom)
{
	int p1, p2, p3;
	if (mutation_operator == "RSM")
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
		reverse(chrom.path.begin() + p1, chrom.path.begin() + p2);
	}
	else 	if (mutation_operator == "PSM")
	{
		p1 = rndint(0, city_num - 1);
		do
		{
			p2 = rndint(0, city_num - 1);
		} while (p1 == p2);
		int city_tmp = chrom.path[p1];
		chrom.path[p1] = chrom.path[p2];
		chrom.path[p2] = city_tmp;
	}

	if (mutation_operator == "EM")		//Exchange mutation   
	{
		p1 = rndint(0, city_num - 1);
		do
		{
			p2 = rndint(0, city_num - 1);
		} while (p1 == p2 || p2 == (p1 - 1 + city_num) % city_num || p2 == (p1 + 1) % city_num);
		int city_tmp = chrom.path[p1];
		chrom.path[p1] = chrom.path[p2];
		chrom.path[p2] = city_tmp;
	}
	else 	if (mutation_operator == "NicheEM")		//Exchange mutation   
	{

		vector<bool> critical_flag(city_num, false);
		for (int j = 0; j < city_num; j++)
		{
			int nj = (j + 1) % city_num;
			if (exist_critical_edge(j, nj))
			{
				critical_flag[j] = true;
				critical_flag[nj] = true;
			}
		}
		vector<int> choose_city_vec;
		for (int j = 0; j < city_num; j++)
		{
			if (critical_flag[j] == false)
				choose_city_vec.push_back(j);
		}

		if (choose_city_vec.size() > 3)
		{
			do
			{
				p1 = rndint(0, choose_city_vec.size() - 1);
				p1 = choose_city_vec[p1];
				int count = 3;
				do
				{
					count--;
					p2 = rndint(0, choose_city_vec.size() - 1);
					p2 = choose_city_vec[p2];
				} while ((p1 == p2 || p2 == (p1 - 1 + city_num) % city_num || p2 == (p1 + 1) % city_num) && count >= 0);
			} while (p1 == p2 || p2 == (p1 - 1 + city_num) % city_num || p2 == (p1 + 1) % city_num);
		}
		else if (choose_city_vec.size() >= 2)
		{
			p1 = choose_city_vec[0];
			p2 = choose_city_vec.back();
		}

		if (choose_city_vec.size() >= 2)
		{

			int city_tmp = chrom.path[p1];
			chrom.path[p1] = chrom.path[p2];
			chrom.path[p2] = city_tmp;
		}

	}
	else 	if (mutation_operator == "DM")		// Displacement mutation  
	{
		p1 = rndint(0, city_num - 1);
		do
		{
			p2 = rndint(0, city_num - 1);
		} while (p1 == p2);
		vector<int> cut_vec;
		vector<int> left_vec;
		if (p1 > p2)
		{
			int tmp = p1;
			p1 = p2;
			p2 = tmp;
		}
		for (int k = 0; k < city_num; k++)
		{
			if (k >= p1 && k <= p2)
				cut_vec.push_back(chrom.path[k]);
			else
				left_vec.push_back(chrom.path[k]);
		}
		if (left_vec.size() == 0)
		{
			p3 = 0;
		}
		else if (left_vec.size() == 1 && p1 == 0)
		{
			p3 = 1;
		}
		else if (left_vec.size() == 1 && p1 == 1)
		{
			p3 = 0;
		}
		else
		{
			do
			{
				p3 = rndint(0, left_vec.size());
			} while (p3 == p1);
		}

		left_vec.insert(left_vec.begin() + p3, cut_vec.begin(), cut_vec.end());
		chrom.path = left_vec;
	}
	chrom.evaluate();
}


void GeneticAlgorithm::report_population(int num, const string & s)
{
	vector<index_len > index_pop(NP);
	index_len index_len_tmp;
	for (int i = 0; i < NP; i++)
	{
		index_len_tmp.index = i;
		index_len_tmp.len = chroms[i].length;
		index_pop[i] = index_len_tmp;
	}
	if (s == "sort")
		sort(index_pop.begin(), index_pop.end(), sort_index_asc_second);
	for (int i = 0; i < num; i++)
	{
		if (s != "sort" && i%SPECIES_SIZE == 0)
		{
			collect_run_distance(i, SPECIES_SIZE);
		}
		cout << chroms[index_pop[i].index].length << " : ";
		for (int j = 0; j < city_num; j++)
		{
			cout << chroms[index_pop[i].index].path[j] << " ";
		}cout << endl;
	}
}

void GeneticAlgorithm::pass_next_generation_speciation_V3(GeneticAlgorithm &Offspring, const int &species_size)
{
	sp_info = Offspring.sp_info;
	index_len tmp_index_len;
	int share_distance_tmp = -1;
	for (int s = 0; s < sp_info.size(); s++)
	{
		for (int i = sp_info[s].begin_index; i <= sp_info[s].end_index; i++)
		{
			int max_share_distance_index = -1;
			int max_share_distance = -1;
			vector<index_len> tmp_comp_vec;
			for (int j = sp_info[s].begin_index; j <= sp_info[s].end_index; j++)
			{
				share_distance_tmp = share_distance(Offspring.chroms[i], chroms[j]);
				if (share_distance_tmp > max_share_distance)
				{
					max_share_distance = share_distance_tmp;

					tmp_index_len.index = j;
					tmp_index_len.len = chroms[j].length;
					tmp_comp_vec.clear();
					tmp_comp_vec.push_back(tmp_index_len);
				}
				else
					if (share_distance_tmp == max_share_distance)
					{
						tmp_index_len.index = j;
						tmp_index_len.len = chroms[j].length;
						tmp_comp_vec.push_back(tmp_index_len);
					}
			}
			sort(tmp_comp_vec.begin(), tmp_comp_vec.end(), sort_index_des_second);
			if (Offspring.chroms[i].length < chroms[tmp_comp_vec[0].index].length)
			{
				chroms[tmp_comp_vec[0].index] = Offspring.chroms[i];
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
void GeneticAlgorithm::selection(GeneticAlgorithm &Parent, const string & s)
{
	int p1, p2;
	if (s == "tournamant")
	{
		for (int i = 0; i < NP; i++)
		{
			p1 = rndint(0, NP - 1);
			do {
				p2 = rndint(0, NP - 1);
			} while (p1 == p2);

			if (chroms[p1].length < chroms[p2].length) {
				chroms[i] = Parent.chroms[p1];
			}
			else {
				chroms[i] = Parent.chroms[p2];
			}
		}
	}
	if (s == "roulette")
	{
		distType overall_len = 0;
		distType max_len = -1;
		for (int i = 0; i < NP; i++)
		{
			overall_len += chroms[i].length;
			max_len = max_len > chroms[i].length ? max_len : chroms[i].length;
		}
		vector<double> probility_vec(NP, 0.0);
		double accum_probibilty = 0.0;
		for (int i = 0; i < NP; i++)
		{
			probility_vec[i] = (max_len - (double)chroms[i].length) / (double)(NP*max_len - overall_len);
			accum_probibilty += probility_vec[i];
		}
		for (int i = 0; i < NP; i++)
		{
			probility_vec[i] /= accum_probibilty;
		}
		for (int i = 1; i < NP; i++)
			probility_vec[i] += probility_vec[i - 1];

		for (int i = 0; i < NP; i++)
		{
			double p = rnddouble(0.0, 1.0);
			int j;
			for (j = NP - 1; j > 0; j--)
			{
				if (p > probility_vec[j])
				{

					break;
				}
			}
			if (j + 1 >= NP)
				j = NP - 1;
			chroms[i] = Parent.chroms[j + 1];
		}
	}
}

void GeneticAlgorithm::speciation_V3_TSPLIB(const GeneticAlgorithm & Parent)
{
	vector<index_len> chroms_index(NP);
	sp_info.clear();
	Species_info species_info_tmp;
	for (int i = 0; i < NP; i++)
	{
		chroms_index[i].index = i;
		chroms_index[i].len = Parent.chroms[i].length;
	}
	sort(chroms_index.begin(), chroms_index.end(), sort_index_asc_second);
	vector<bool> flag_selected(NP, false);

	int min_len = chroms_index[0].len;
	int mid_len = chroms_index[NP / 2].len;

	index_len idx_len;
	int tmp_species_size = 0;
	for (int i = 0; i < NP; i += tmp_species_size)
	{


		int seed_index = 0;
		while (flag_selected[chroms_index[seed_index].index] == true)
		{
			seed_index++;
		}
		flag_selected[chroms_index[seed_index].index] = true;
		chroms[i] = Parent.chroms[chroms_index[seed_index].index];
		if (30 + i >= NP)
			tmp_species_size = NP - i;
		else
		{
			double tmp_ratio = (chroms[i].length - min_len) / (double)(mid_len - min_len);
			tmp_ratio = tmp_ratio < 1.0 ? tmp_ratio : 1.0;
			tmp_species_size = floor(tmp_ratio*10.0) + 20;
			if (tmp_species_size % 2 == 1)
				tmp_species_size++;
		}

		species_info_tmp.begin_index = i;
		species_info_tmp.end_index = i + tmp_species_size - 1;
		species_info_tmp.species_size = tmp_species_size;
		sp_info.push_back(species_info_tmp);
		vector<index_len> chroms_index_share_distance;

		for (int j = 0; j < NP; j++)
		{
			int current_index = chroms_index[j].index;
			if (flag_selected[current_index] != true)
			{
				idx_len.index = current_index;
				idx_len.len = city_num - share_distance(chroms[i], Parent.chroms[current_index]);
				chroms_index_share_distance.push_back(idx_len);
			}
		}
		sort(chroms_index_share_distance.begin(), chroms_index_share_distance.end(), sort_index_asc_second);

		int count_share = 0;
		for (int j = 0; j < tmp_species_size - 1; j++)
		{
			chroms[i + j + 1] = Parent.chroms[chroms_index_share_distance[j].index];
			flag_selected[chroms_index_share_distance[j].index] = true;
			count_share += chroms_index_share_distance[j].len;
		}
		if (count_share == 0)
		{
			for (int j = 0; j < tmp_species_size - 1; j++)
			{
				mutation(chroms[i + j + 1]);
			}
		}
		else
		{
			random_shuffle(chroms.begin() + i + 1, chroms.begin() + i + tmp_species_size - 1);
		}
	}
}
void GeneticAlgorithm::speciation_V3(const GeneticAlgorithm & Parent, const int min_m, const int max_m)
{
	int total_count = 0;
	vector<index_len> chroms_index(NP);
	sp_info.clear();
	Species_info species_info_tmp;
	for (int i = 0; i < NP; i++)
	{
		chroms_index[i].index = i;
		chroms_index[i].len = Parent.chroms[i].length;
	}
	sort(chroms_index.begin(), chroms_index.end(), sort_index_asc_second);
	vector<bool> flag_selected(NP, false);

	int min_len = chroms_index[0].len;

	int beta_position = chroms.size() / min_m;
	int mid_len;
	if ((chroms.size() / min_m - beta_position) > 0.1)
		mid_len = chroms_index[beta_position - 1].len;
	else
		mid_len = (chroms_index[beta_position - 1].len + chroms_index[beta_position].len) / 2.0;

	index_len idx_len;
	int tmp_species_size = 0;
	for (int i = 0; i < NP; i += tmp_species_size)
	{
		int seed_index = 0;
		while (flag_selected[chroms_index[seed_index].index] == true)
		{
			seed_index++;
		}
		flag_selected[chroms_index[seed_index].index] = true;
		chroms[i] = Parent.chroms[chroms_index[seed_index].index];
		if (max_m + i >= NP)
			tmp_species_size = NP - i;
		else
		{
			double tmp_ratio = (chroms[i].length - min_len) / (double)(mid_len - min_len);
			tmp_ratio = tmp_ratio < 1.0 ? tmp_ratio : 1.0;
			tmp_species_size = floor(tmp_ratio*(max_m - min_m)) + min_m;
			if (tmp_species_size % 2 == 1)
				tmp_species_size++;
		}

		species_info_tmp.begin_index = i;
		species_info_tmp.end_index = i + tmp_species_size - 1;
		species_info_tmp.species_size = tmp_species_size;
		sp_info.push_back(species_info_tmp);

		vector<index_len> chroms_index_share_distance;

		for (int j = 0; j < NP; j++)
		{
			int current_index = chroms_index[j].index;
			if (flag_selected[current_index] != true)
			{
				idx_len.index = current_index;
				idx_len.len = city_num - share_distance(chroms[i], Parent.chroms[current_index]);
				chroms_index_share_distance.push_back(idx_len);
			}
		}
		sort(chroms_index_share_distance.begin(), chroms_index_share_distance.end(), sort_index_asc_second);
		int count_share = 0;
		for (int j = 0; j < tmp_species_size - 1; j++)
		{
			chroms[i + j + 1] = Parent.chroms[chroms_index_share_distance[j].index];
			flag_selected[chroms_index_share_distance[j].index] = true;
			count_share += chroms_index_share_distance[j].len;
		}
		if (count_share == 0)
		{
			for (int j = 0; j < tmp_species_size - 1; j++)
			{

				Chromosome  &chrom = chroms[i + j + 1];
				int p1, p2;
				p1 = rndint(0, city_num - 1);
				do
				{
					p2 = rndint(0, city_num - 1);
				} while (p1 == p2 || p2 == (p1 - 1 + city_num) % city_num || p2 == (p1 + 1) % city_num);

				distType prove = dist[chrom.path[(p1 - 1 + city_num) % city_num]][chrom.path[p1]] + dist[chrom.path[p1]][chrom.path[(p1 + 1) % city_num]]
					+ dist[chrom.path[(p2 - 1 + city_num) % city_num]][chrom.path[p2]] + dist[chrom.path[p2]][chrom.path[(p2 + 1) % city_num]];
				int city_tmp = chroms[i + j + 1].path[p1];
				chrom.path[p1] = chrom.path[p2];
				chrom.path[p2] = city_tmp;
				prove -= dist[chrom.path[(p1 - 1 + city_num) % city_num]][chrom.path[p1]] + dist[chrom.path[p1]][chrom.path[(p1 + 1) % city_num]]
					+ dist[chrom.path[(p2 - 1 + city_num) % city_num]][chrom.path[p2]] + dist[chrom.path[p2]][chrom.path[(p2 + 1) % city_num]];
				chrom.length -= prove;
				total_count += 8;
			}
		}
		else
		{
			random_shuffle(chroms.begin() + i + 1, chroms.begin() + i + tmp_species_size - 1);
		}


	}
	FES += (total_count / city_num);
}
GeneticAlgorithm::GeneticAlgorithm(int np, int cn, double cr, double mr, int MAX_ITER, City c, long long maxfes, int sps_s, double init_r, const string &cr_op, const string &mt_op, const string &lc_op)
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
	specie_size = sps_s;
	init_ratio = init_r;
	crossover_opertor = cr_op;
	origin_crossover_operator = crossover_opertor;
	mutation_operator = mt_op;
	local_search_operator = lc_op;
	max_critical_threshold_gen = 1;

	STAGNT_NUM = city_num;
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
	maxFES = b.maxFES;
	specie_size = b.specie_size;
	init_ratio = b.init_ratio;
	crossover_opertor = b.crossover_opertor;
	origin_crossover_operator = b.origin_crossover_operator;
	mutation_operator = b.mutation_operator;
	local_search_operator = b.local_search_operator;
	max_critical_threshold_gen = b.max_critical_threshold_gen;
}


void GeneticAlgorithm::initiate_V2(Chromosome & chrom)
{
	vector<bool> visited(city_num, false);
	int selected_city;
	int count = 0;
	do {

		selected_city = rndint(0, city_num - 1);
		if (visited[selected_city] == false)
		{
			chrom.path[count] = selected_city;
			visited[selected_city] = true;
			count++;
		}
	} while (count < init_ratio *city_num);

	for (int i = count; i < city_num; i++)
	{
		for (int j = 1; j < city_num; j++)
		{
			if (visited[neighbor_matrix[chrom.path[i - 1]][j]] == false)
			{
				selected_city = neighbor_matrix[chrom.path[i - 1]][j];
				break;
			}
		}
		visited[selected_city] = true;
		chrom.path[i] = selected_city;
	}
}
void GeneticAlgorithm::initiate_V2()
{
	vector<int> path_tmp(city_num, -1);

	for (int i = 0; i < NP; i++)
	{
		Chromosome chrom(path_tmp);
		initiate_V2(chrom);
		chroms.push_back(chrom);
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
	distType min_len;
public:
	Problem(City &c, vector<vector<int> > &s, const distType &l, string & n);
	Problem(City &c, const string &problem_name);
	friend class GeneticAlgorithm;
};
Problem::Problem(City &c, vector<vector<int> > &s, const distType &l, string & n) {
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
	int NP = atoi(argv[4]);
	double  crossover_ratio = atof(argv[5]);
	double mutation_ratio = atof(argv[6]);
	output_dir = "./" + string(argv[7]) + "/";
	string selection_operator = string(argv[8]);
	string crossover_operator = string(argv[9]);
	string mutation_operator = string(argv[10]);
	string ls_operator = string(argv[11]);
	int species_size = -1;
	int MAX_ITER = -1;

	double init_ratio = 0.5;
	double refine_niche_radius = 0.0;
	int track_gap;

	string file_string = output_dir + "info.param";
	ofstream ofile_record;
	convert_file_name(file_string);
	ofile_record.open(file_name_char, ios::out);
	ofile_record << "NP: " << NP << endl
		<< "RUNS: " << RUNS << endl
		<< "crossove Rate: " << crossover_ratio << endl
		<< "crossove Operator: " << crossover_operator << endl
		<< "mutation Rate: " << mutation_ratio << endl
		<< "mutation Operator: " << mutation_operator << endl
		<< "local Search: " << ls_operator << endl;
	ofile_record.close();

	if (crossover_operator != "critical_PMX")
		return -21;
	if (mutation_operator != "NicheEM")
		return -22;

	srand(time(NULL));
	string problem_name;

	for (int func = FUNC_BEGIN; func <= FUNC_END; func++)
	{
		problem_name = Func_name[func];
		if (problem_name == "")
			continue;
		City cities(problem_name);
		Problem TSP_pro(cities, problem_name);
		string fs_file_time;
		fstream fs_time;
		clock_t  clockBegin, clockEnd, clockEnd_gap;
		string file_string_time = output_dir + problem_name + ".alg_time";
		ofstream ofile_record_time;
		convert_file_name(file_string_time);
		ofile_record_time.open(file_name_char, ios::out);
		if (ofile_record_time.is_open() == false)
		{
			cout << file_name_char << " not found" << endl;
			return -1;
		}

		if (func >= 0 && func < 29)
		{
			MAX_ITER = MaxFES[func] / NP;
			refine_niche_radius = 0.01;
			init_ratio = 0.5;
		}
		else	if (func >= 30)
		{
			MaxFES[func] = cities.city_num * 50000;
			MAX_ITER = MaxFES[func] / NP;
			refine_niche_radius = 0.1;
			init_ratio = 0.5;
		}
		int NP_output = 16;
		for (int run = 0; run < RUNS; run++)
		{
			cout << problem_name << ": Begin " << run << endl;

			int track_count = 1;
			clockBegin = clock();
			int time_gap = 0;

			GeneticAlgorithm Parent(NP, cities.city_num, crossover_ratio, mutation_ratio, MAX_ITER, cities, MaxFES[func], species_size, init_ratio,
				crossover_operator, mutation_operator, ls_operator);

			Parent.initiate_V2();
			Parent.length_evaluate();
			GeneticAlgorithm Offspring = Parent;

			Offspring.empty_edge_matrix();
			while (!Offspring.terminate())
			{
				Offspring.speciation_V3(Parent, 4, 12);
				Offspring.record_edge_matrix();
				Offspring.update_edge_matrix();
				Offspring.crossover();
				Offspring.mutation();
				Offspring.length_evaluate();
				Offspring.localSearch_V3_Critical();
				Parent.pass_next_generation_speciation_V3(Offspring, species_size);
				Parent.update_sucessive();
			}
			clockEnd = clock();
			ofile_record_time << (double)(clockEnd - clockBegin) / CLOCKS_PER_SEC;
			ofile_record_time << endl;

			Parent.remove_redudant_report2_file_V6(refine_niche_radius, run, 0, problem_name, species_size);

			cout << problem_name << ": Finished " << run << endl;
		} // for (int run = 0; run < RUNS; run++)
		ofile_record_time.close();
	}
#ifdef _MSC_VER
	system("pause");
#endif
	return 0;
}
