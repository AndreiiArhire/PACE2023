/* EXACT */
// TO DO te uiti in componente conexe si daca upperbound < max din alta componenta lasi asa 
#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,fma")
#pragma GCC optimization("unroll-loops")
#include <iostream>
#include <fstream>
#include <string>
#include <deque>
#include <algorithm>
#include <sstream>
#include <unordered_set>
#include <random>       
#include <chrono>
#include <map>       
#include <unordered_map>       
#include <set>       
#include <bitset>

using namespace std;
int testNo, n, m;
const int MOD1 = 100000007;
const int PRIME1 = 10000019;
const int MOD2 = 100000037;
const int PRIME2 = 10000079;
const int EXIT_CODE = 1;
std::chrono::time_point<std::chrono::high_resolution_clock> begin_;


void ckeck_force_stop() {
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin_);
	double sec = elapsed.count() * 1e-9;
	if (sec > 30 * 60 - 10) {
		exit(EXIT_CODE);
	}
}



set<int> all_nodes_set;
vector<vector<int>> components;
vector<unordered_set<int>>edges;
vector<unordered_set<int>>red;
vector<unordered_set<int>>red_bound;
vector<unordered_set<int>>black;
vector<unordered_set<int>>black_bound;
vector<int> paired;
vector<int> paired2;
vector<int> poss;
string s;
vector<pair<int, int>>  best_solution, best_solution_all_graph;
vector<int> nodes, ordered_nodes;
unordered_set<int> neighbors, temp, nodes2;
vector < pair < int, pair < int, int > > > choices;
vector<int> bestDegree;
int best_red_degree = (int)1e9, best_red_degree_all_graph;
vector<int>temp_vector;
vector<  pair <int, int > > temp1;
vector< pair <int, int> > hashes;
vector<vector<char>>adjacency;
int STEP, LOWER_BOUND = -1, SW = 0;
struct partial_solution_hash {
	vector < pair <int, int> > node_pairs;
	int max_red_degree = 0;
	int number_of_red_edges = 0;
};

vector<partial_solution_hash> current_hash_solutions[2];

unordered_set<int> all_nodes, current_nodes;

vector<int> t1, t2, t3;
unordered_set<int> nbr;
vector<pair <int, int >> t4;
vector< vector<int> > temp_edges;
vector<int> temp_visited;
int lower_bound_tl;
struct pair_hash {
	template <class T1, class T2>
	std::size_t operator () (const std::pair<T1, T2>& p) const {
		auto a = std::hash<T1>{}(p.first);
		auto b = std::hash<T2>{}(p.second);
		// Cantor function for pairing hash
		return  1ll * (1ll * (a + b) * (a + b + 1) / 2 + 1ll * b) % 1000000009;
	}
};


unordered_map < pair <int, int>, int, pair_hash > hash_to_index;
//map < pair <int, int>, int > hash_to_index;

double getElapsed() {
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin_);
	double sec = elapsed.count() * 1e-9;
	return sec;
}



void init_data() {
	edges.resize(n + 5, unordered_set <int>());
	temp_edges.resize(n + 5, vector <int>());
	red.resize(n + 5, unordered_set <int>());
	red_bound.resize(n + 5, unordered_set <int>());
	black.resize(n + 5, unordered_set <int>());
	black_bound.resize(n + 5, unordered_set <int>());
	paired.resize(n + 5, 0);
	hashes.resize(3 * n + 5, make_pair(PRIME1, PRIME2));
	poss.resize(n + 5, 0);
	paired2.resize(n + 5, 0);
	temp_visited.resize(n + 5, 0);
	adjacency.resize(n + 5, vector<char>(n + 5, 0));
	for (int i = 1; i <= 3 * n; ++i) {
		hashes[i].first = 1ll * hashes[i - 1].first * PRIME1 % MOD1;
		hashes[i].second = 1ll * hashes[i - 1].second * PRIME2 % MOD2;
	}
}

void clear_data() {
	edges.clear();
	red.clear();
	black.clear();
	paired.clear();
	paired2.clear();
	nodes.clear();
	ordered_nodes.clear();
	neighbors.clear();
	choices.clear();
}

void readData() {
	string path_input = R"(tests/exact/exact_)";
	for (int l = 1; l <= 3 - to_string(testNo).size(); ++l) {
		path_input.push_back('0');
	}
	path_input += to_string(testNo);
	path_input += ".gr";
	//ifstream in(path_input);
	cin.tie();
	cout.tie();
	string problem_id;
	cin >> problem_id >> problem_id >> n >> m;
	init_data();
	for (int j = 1; j <= m; ++j) {
		int x = 0, y = 0;
		cin >> x >> y;
		edges[x].insert(y);
		edges[y].insert(x);
	}
	//in.close();
	if (n == 0) {
		exit(0);
	}
}


void split_into_components() {
	for (int i = 1; i <= n; ++i) {
		paired[i] = 0;
	}
	for (int i = 1; i <= n; ++i) {
		if (!paired[i]) {
			vector<int> curr_component;
			deque <int> dq;
			dq.push_back(i);
			while (!dq.empty()) {
				int nod = dq.front();
				paired[nod] = 1;
				curr_component.emplace_back(nod);
				dq.pop_front();
				for (auto it : edges[nod]) {
					if (!paired[it]) {
						paired[it] = 1;
						dq.push_back(it);
					}
				}
			}
			components.emplace_back(curr_component);
		}
	}
}

pair <int, vector<pair<int, int>>> get_initial_upper_bound(int cc) {
	vector< pair <int, int >> solution;
	int max_degree = 0, round = 0;
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	solution.clear();
	for (auto i : components[cc]) {
		black[i] = edges[i];
		red[i] = unordered_set<int>();
		nodes.push_back(i);
		nodes2.insert(i);
		paired[i] = 0;
		paired2[i] = 0;
	}
	while (nodes.size() > 1) {
		ckeck_force_stop();
		pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > best_option = { { {1e9, 1e9}, {1e9, 1e9} } , {1e9, 1e9} };
		for (int i = 0; i < nodes.size(); ++i) {
			for (int j = i + 1; j < nodes.size(); ++j) {
				int x = nodes[i];
				int y = nodes[j];
				int red_between = 0;
				if (red[x].count(y)) {
					red_between = 1;
				}
				if (((int)min(red[x].size(), red[y].size())) - red_between + abs(((int)(black[x].size() + red[x].size())) - ((int)(black[y].size() + red[y].size()))) > best_option.first.first.first) {
					continue;
				}

				if (black[x].size() + red[x].size() > black[y].size() + red[y].size()) {
					swap(x, y);
				}
				int common_both_black = 0, common_mix = 0, common_both_red = 0, edges_between = 0, max_red_degree = 0, new_red_change = 0, new_red_not_change = 0;
				int cnt_max = black[x].size() + red[x].size();
				if (black[x].count(y) || red[x].count(y)) {
					edges_between = 2;
				}
				for (auto it : black[x]) {
					if (black[y].count(it)) {
						++common_both_black;
					}
					if (red[y].count(it)) {
						++common_mix;
					}
					if (it != y && !red[y].count(it) && !black[y].count(it))
					{
						max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : black[y]) {
					if (it != x && !red[x].count(it) && !black[x].count(it)) {
						max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : red[x]) {
					if (black[y].count(it)) {
						++common_mix;
					}
					if (red[y].count(it)) {
						++common_both_red;
					}
				}
				max_red_degree = max(max_red_degree, (int)(red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between));
				new_red_not_change = (int)(red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between) - common_both_red;
				pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > curr_option = { { { max_red_degree, new_red_change }, {new_red_not_change, 0} }, { x, y } };
				best_option = min(best_option, curr_option);
			}
		}
		int x = best_option.second.first;
		int y = best_option.second.second;

		if (red[x].size() < red[y].size()) {
			swap(x, y);
		}
		red[x].erase(y);
		red[y].erase(x);
		black[y].erase(x);
		black[x].erase(y);
		for (auto it : red[y]) {
			red[it].erase(y);
			red[x].insert(it);
			red[it].insert(x);
			max_degree = max(max_degree, (int)red[it].size());
		}
		for (auto it : black[x]) {
			if (black[y].count(it)) {
				temp.insert(it);
			}
			else
			{
				black[it].erase(x);
				red[x].insert(it);
				red[it].insert(x);
				max_degree = max(max_degree, (int)red[it].size());
			}
		}
		for (auto it : black[y]) {
			if (!temp.count(it)) {
				red[x].insert(it);
				red[it].insert(x);
				max_degree = max(max_degree, (int)red[it].size());
			}
			black[it].erase(y);
		}
		black[x] = temp;
		red[y] = unordered_set<int>();
		black[y] = unordered_set<int>();
		max_degree = max(max_degree, best_option.first.first.first);
		solution.emplace_back(x, y);
		nodes2.erase(y);
		temp = unordered_set<int>();
		nodes.clear();
		for (auto it : nodes2) {
			nodes.emplace_back(it);
		}
		if (nodes.size() - 1 < max_degree) {
			while (nodes.size() > 1) {
				solution.emplace_back(nodes[0], nodes.back());
				nodes.pop_back();
			}
		}
	}
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	return make_pair(max_degree, solution);
}


pair < int, vector< pair <int, int > > > solve_bound(int cc, int hash_id, int step) {
	ckeck_force_stop();
	vector<pair <int, int >>solution;
	int max_degree = 0, round = 0;

	nodes = vector<int>();
	nodes2 = unordered_set<int>();

	for (auto i : components[cc]) {
		black_bound[i] = black[i];
		red_bound[i] = red[i];
		nodes2.insert(i);
	}

	for (auto it : current_hash_solutions[step][hash_id].node_pairs) {
		nodes2.erase(it.second);
	}
	for (auto it : nodes2) {
		nodes.emplace_back(it);
	}

	while (nodes.size() > 1) {
		ckeck_force_stop();
		pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > best_option = { { {1e9, 1e9}, {1e9, 1e9} } , {1e9, 1e9} };
		for (int i = 0; i < nodes.size(); ++i) {
			for (int j = i + 1; j < nodes.size(); ++j) {
				int x = nodes[i];
				int y = nodes[j];
				int red_between = 0;
				if (red_bound[x].count(y)) {
					red_between = 1;
				}
				if (((int)min(red_bound[x].size(), red_bound[y].size())) - red_between + abs(((int)(black_bound[x].size() + red_bound[x].size())) - ((int)(black_bound[y].size() + red_bound[y].size()))) > best_option.first.first.first) {
					continue;
				}
				if (black_bound[x].size() + red_bound[x].size() > black_bound[y].size() + red_bound[y].size()) {
					swap(x, y);
				}
				int common_both_black = 0, common_mix = 0, common_both_red = 0, edges_between = 0, max_red_degree = 0, new_red_change = 0, new_red_not_change = 0;
				int cnt_max = black_bound[x].size() + red_bound[x].size();
				if (black_bound[x].count(y) || red_bound[x].count(y)) {
					edges_between = 2;
				}
				for (auto it : black_bound[x]) {
					if (black_bound[y].count(it)) {
						++common_both_black;
					}
					if (red_bound[y].count(it)) {
						++common_mix;
					}
					if (it != y && !red_bound[y].count(it) && !black_bound[y].count(it))
					{
						max_red_degree = max(max_red_degree, (int)(red_bound[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : black_bound[y]) {
					if (it != x && !red_bound[x].count(it) && !black_bound[x].count(it)) {
						max_red_degree = max(max_red_degree, (int)(red_bound[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : red_bound[x]) {
					if (black_bound[y].count(it)) {
						++common_mix;
					}
					if (red_bound[y].count(it)) {
						++common_both_red;
					}
				}
				max_red_degree = max(max_red_degree, (int)(red_bound[x].size() + black_bound[x].size() + red_bound[y].size() + black_bound[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between));
				new_red_not_change = (int)(red_bound[x].size() + black_bound[x].size() + red_bound[y].size() + black_bound[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between) - common_both_red;
				pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > curr_option = { { { max_red_degree, new_red_change }, {new_red_not_change, 0} }, { x, y } };
				best_option = min(best_option, curr_option);
			}
		}
		int x = best_option.second.first;
		int y = best_option.second.second;

		if (red_bound[x].size() < red_bound[y].size()) {
			swap(x, y);
		}
		red_bound[x].erase(y);
		red_bound[y].erase(x);
		black_bound[y].erase(x);
		black_bound[x].erase(y);
		for (auto it : red_bound[y]) {
			red_bound[it].erase(y);
			red_bound[x].insert(it);
			red_bound[it].insert(x);
			max_degree = max(max_degree, (int)red_bound[it].size());
		}
		for (auto it : black_bound[x]) {
			if (black_bound[y].count(it)) {
				temp.insert(it);
			}
			else
			{
				black_bound[it].erase(x);
				red_bound[x].insert(it);
				red_bound[it].insert(x);
				max_degree = max(max_degree, (int)red_bound[it].size());
			}
		}
		for (auto it : black_bound[y]) {
			if (!temp.count(it)) {
				red_bound[x].insert(it);
				red_bound[it].insert(x);
				max_degree = max(max_degree, (int)red_bound[it].size());
			}
			black_bound[it].erase(y);
		}
		black_bound[x] = temp;
		red_bound[y] = unordered_set<int>();
		black_bound[y] = unordered_set<int>();
		max_degree = max(max_degree, best_option.first.first.first);
		solution.emplace_back(x, y);
		nodes2.erase(y);
		temp = unordered_set<int>();
		nodes.clear();
		for (auto it : nodes2) {
			nodes.emplace_back(it);
		}
		if (nodes.size() - 1 < max_degree) {
			while (nodes.size() > 1) {
				solution.emplace_back(nodes[0], nodes.back());
				nodes.pop_back();
			}
		}
	}

	for (auto i : components[cc]) {
		black_bound[i] = unordered_set<int>();
		red_bound[i] = unordered_set<int>();
	}
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	ckeck_force_stop();
	return make_pair(max_degree, solution);
}

inline int merge_nodes(const pair<int, int>& node_pair) {
	temp.clear();
	int x = node_pair.first, y = node_pair.second, max_red = 0;
	red[x].erase(y);
	red[y].erase(x);
	black[y].erase(x);
	black[x].erase(y);
	for (auto it : red[y]) {
		red[it].erase(y);
		red[x].insert(it);
		red[it].insert(x);
	}
	for (auto it : black[x]) {
		if (black[y].count(it)) {
			temp.insert(it);
		}
		else
		{
			black[it].erase(x);
			red[x].insert(it);
			red[it].insert(x);
			max_red = max(max_red, (int)red[it].size());
		}
	}
	for (auto it : black[y]) {
		if (!temp.count(it)) {
			red[x].insert(it);
			red[it].insert(x);
			max_red = max(max_red, (int)red[it].size());
		}
		black[it].erase(y);
	}
	max_red = max(max_red, (int)red[x].size());
	black[x] = temp;
	red[y] = unordered_set<int>();
	black[y] = unordered_set<int>();
	temp = unordered_set <int>();
	return max_red;
}


pair <int, int> get_solution_hash(const int& hash_id, const int& step, const pair <int, int>& node_pair) {
	pair <int, int> ret = { 0, 0 };
	temp_edges[node_pair.first].emplace_back(node_pair.second);
	temp_edges[node_pair.second].emplace_back(node_pair.first);
	int new_node1 = -1, new_node2 = -1;
	if (!all_nodes_set.count(node_pair.first)) {
		new_node1 = node_pair.first;
		all_nodes_set.insert(node_pair.first);
	}
	if (!all_nodes_set.count(node_pair.second)) {
		new_node2 = node_pair.second;
		all_nodes_set.insert(node_pair.second);
	}

	vector<int> sequence;
	deque <int> dq;
	for (auto it : all_nodes_set) {
		if (!temp_visited[it]) {
			vector<int> component;
			dq.emplace_back(it);
			temp_visited[it] = 1;
			while (!dq.empty()) {
				int nod = dq.front();
				dq.pop_front();
				component.emplace_back(nod);
				for (auto it : temp_edges[nod]) {
					if (!temp_visited[it]) {
						dq.push_back(it);
						temp_visited[it] = 1;
					}
				}
			}
			sort(component.begin(), component.end());
			for (auto it : component) {
				sequence.emplace_back(it);
			}
			sequence.emplace_back(0);
		}
	}

	for (int i = 0; i < sequence.size(); ++i) {
		temp_visited[sequence[i]] = 0;
		ret.first = 1ll * (1ll * ret.first + 1ll * sequence[i] * hashes[i].first) % MOD1;
		ret.second = 1ll * (1ll * ret.second + 1ll * sequence[i] * hashes[i].second) % MOD2;
	}
	all_nodes_set.erase(new_node1);
	all_nodes_set.erase(new_node2);
	temp_edges[node_pair.first].pop_back();
	temp_edges[node_pair.second].pop_back();
	return ret;
}

int max_clusters;
unordered_set <int> clusters;
unordered_set<int> nodes_left;


int build_graph_from_hash(const int& cc, const int& hash_id, const int& step) {
	int to_be_continued = 1;
	clusters = unordered_set <int>();
	for (auto i : components[cc]) {
		black[i] = edges[i];
		red[i] = unordered_set<int>();
		nodes_left.insert(i);
	}
	for (auto node_pair : current_hash_solutions[step][hash_id].node_pairs) {
		merge_nodes(node_pair);
		clusters.insert(node_pair.first);
		clusters.erase(node_pair.second);
		nodes_left.erase(node_pair.second);
	}
	if (clusters.size() > max_clusters) {
		clusters = unordered_set <int>();
		to_be_continued = 0;
		return to_be_continued;
	}
	if (SW == 1) {
		if (nodes_left.size() == 0) {
			if (current_hash_solutions[step][hash_id].max_red_degree < best_red_degree) {
				best_red_degree = current_hash_solutions[step][hash_id].max_red_degree;
				best_solution.clear();
				for (auto it : current_hash_solutions[step][hash_id].node_pairs) {
					best_solution.emplace_back(it);
				}
			}
		}
		return 1;
	}
	ckeck_force_stop();
	nodes_left = unordered_set <int>();
	pair < int, vector<pair <int, int>> > ret = make_pair(1e9, vector<pair<int, int>>());
	if ((LOWER_BOUND == -1 && current_hash_solutions[step][hash_id].max_red_degree >= best_red_degree - 1 /* && STEP >= 2*/) || (LOWER_BOUND != -1 && current_hash_solutions[step][hash_id].max_red_degree >= LOWER_BOUND)) {
		ret = solve_bound(cc, hash_id, step);
		int upper_bound = ret.first;
		int max_red_degree = max(ret.first, current_hash_solutions[step][hash_id].max_red_degree);
		if (max_red_degree < best_red_degree || (best_solution.empty() && max_red_degree == best_red_degree)) {
			best_red_degree = max_red_degree;
			best_solution.clear();
			for (auto it : current_hash_solutions[step][hash_id].node_pairs) {
				best_solution.emplace_back(it);
			}
			for (auto it : ret.second) {
				best_solution.emplace_back(it);
			}
		}
		if (upper_bound <= current_hash_solutions[step][hash_id].max_red_degree) {
			to_be_continued = 0;
		}
	}
	return to_be_continued;
}
int dist_cel_mult_2, muchie_rosie_creata;
pair <int, int> get_red_degree_by_simulation(int x, int y) {
	int  max_red_degree = 0, val = 0, muchii_rosii_noi = 0;
	dist_cel_mult_2 = 0;
	muchie_rosie_creata = 0;
	// x fixat -> 1
	// y fixat -> 2
	int common_both_black = 0, common_mix = 0, common_both_red = 0, edges_between = 0;
	int cnt_max = black[x].size() + red[x].size();
	if (red[x].count(y)) {
		muchii_rosii_noi--;
	}

	if (black[x].count(y) || red[x].count(y)) {
		dist_cel_mult_2 = 1;
		edges_between = 2;
	}
	for (auto it : black[x]) {
		if (black[y].count(it)) {
			dist_cel_mult_2 = 1;
			++common_both_black;
		}
		if (red[y].count(it)) {
			dist_cel_mult_2 = 1;
			++common_mix;
		}
		if (it != y && !red[y].count(it) && !black[y].count(it))
		{
			muchii_rosii_noi++;
			muchie_rosie_creata = 1;
			max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
		}
	}
	for (auto it : black[y]) {
		if (it != x && !red[x].count(it) && !black[x].count(it))
		{
			muchii_rosii_noi++;
			muchie_rosie_creata = 1;
			max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
		}
	}
	for (auto it : red[x]) {
		if (black[y].count(it)) {
			dist_cel_mult_2 = 1;
			++common_mix;
			val |= 2;
		}
		if (red[y].count(it)) {
			muchii_rosii_noi--;
			dist_cel_mult_2 = 1;
			++common_both_red;
		}
		if (it != y && !red[y].count(it))
		{
			val |= 2;
		}
	}
	for (auto it : red[y]) {
		if (it != x && !red[x].count(it)) {
			val |= 1;
		}
	}
	if (val == 3) {
		muchie_rosie_creata = 1;
	}
	max_red_degree = max(max_red_degree, (int)(red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between));
	return make_pair(max_red_degree, muchii_rosii_noi);
}


vector< pair <int, int> > reduce_component(int cc) {
	vector <pair <int, int>> solution;
	int change = 1;
	while (change) {
		change = 0;
		sort(components[cc].begin(), components[cc].end());
		for (auto node : components[cc]) {
			paired[node] = 0;
		}

		vector<pair <int, int>> to_be_united;
		set<int> not_erased;
		for (auto node : components[cc]) {
			if (paired[node]) {
				continue;
			}
			for (auto it : components[cc]) {
				not_erased.insert(it);
			}
			for (auto it : edges[node]) {
				if (it > node && !paired[it]) {
					set<int> nbrs1, nbrs2;
					for (auto it2 : edges[node]) {
						if (it2 != it) {
							nbrs1.insert(it2);
						}
					}
					for (auto it2 : edges[it]) {
						if (it2 != node) {
							nbrs2.insert(it2);
						}
					}
					if (nbrs1 == nbrs2) {
						paired[it] = 1;
						to_be_united.emplace_back(node, it);
					}
				}
			}
		}
		for (auto it : to_be_united) {
			solution.emplace_back(it);
			for (auto it2 : edges[it.second]) {
				edges[it2].erase(it.second);
			}
			edges[it.second] = unordered_set<int>();
			not_erased.erase(it.second);
			change = 1;
		}
		components[cc].clear();
		for (auto it : not_erased) {
			components[cc].emplace_back(it);
		}
		sort(components[cc].begin(), components[cc].end());
		vector< pair < vector<int>, int > > edges_nodes;
		for (auto node : components[cc]) {
			vector<int> curr;
			for (auto it : edges[node]) {
				curr.emplace_back(it);
			}
			sort(curr.begin(), curr.end());
			edges_nodes.emplace_back(curr, node);
		}
		sort(edges_nodes.begin(), edges_nodes.end());
		for (int i = 0; i + 1 < edges_nodes.size(); ++i) {
			int j = i + 1;
			for (; j < edges_nodes.size() && edges_nodes[j].first == edges_nodes[i].first; ++j) {
				auto it = make_pair(edges_nodes[i].second, edges_nodes[j].second);
				solution.emplace_back(it);
				for (auto it2 : edges[it.second]) {
					edges[it2].erase(it.second);
				}
				edges[it.second] = unordered_set<int>();
				not_erased.erase(it.second);
				change = 1;
			}
			i = j - 1;
		}
		components[cc].clear();
		for (auto it : not_erased) {
			components[cc].emplace_back(it);
		}
	}
	return solution;
}


vector<pair <int, int>> edge_pairs;
vector < int> deg_vec;
vector < int> red_edge;
unordered_set<int> new_red_edges;

void branch_and_bound(int cc) {
	partial_solution_hash emprty_solution;
	current_hash_solutions[0].emplace_back(emprty_solution);
	STEP = 0;
	ckeck_force_stop();
	for (char step = 0; !current_hash_solutions[step].empty() && best_red_degree > LOWER_BOUND; step ^= 1, ++STEP) {
		ckeck_force_stop();
		if (SW == 1) {
			sort(current_hash_solutions[step].begin(), current_hash_solutions[step].end(), [](partial_solution_hash i, partial_solution_hash j) { return i.max_red_degree < j.max_red_degree || (i.max_red_degree == j.max_red_degree && i.number_of_red_edges < j.number_of_red_edges) || (i.max_red_degree == j.max_red_degree && i.number_of_red_edges == j.number_of_red_edges && i.node_pairs < j.node_pairs); });
			while (current_hash_solutions[step].size() > 1000) {
				current_hash_solutions[step].pop_back();
			}
		}

		for (int hash_id = 0; hash_id < current_hash_solutions[step].size() && best_red_degree > LOWER_BOUND; ++hash_id) {
			auto solution_hash = current_hash_solutions[step][hash_id];

			ckeck_force_stop();
			if (current_hash_solutions[step][hash_id].max_red_degree >= best_red_degree) {
				continue;
			}

			if (!build_graph_from_hash(cc, hash_id, step)) {
				if (getElapsed() > lower_bound_tl) {
					break;
				}
				continue;
			}
			int number_of_red_edges = 0;
			for (auto i : components[cc]) {
				number_of_red_edges += red_bound[i].size();
			}
			number_of_red_edges /= 2;
			for (auto node_pair : current_hash_solutions[step][hash_id].node_pairs) {
				all_nodes_set.insert(node_pair.first);
				all_nodes_set.insert(node_pair.second);
				temp_edges[node_pair.first].emplace_back(node_pair.second);
				temp_edges[node_pair.second].emplace_back(node_pair.first);
			}

			current_nodes = unordered_set<int>();
			for (auto it : components[cc]) {
				current_nodes.insert(it);
			}

			for (auto node : solution_hash.node_pairs) {
				current_nodes.erase(node.second);
			}
			for (auto it : current_nodes) {
				t1.emplace_back(it);				// nu sunt
			}
			sort(t1.begin(), t1.end());

			for (int i = 0; i < t1.size(); ++i) {
				for (int j = i + 1; j < t1.size(); ++j) {
					auto r = get_red_degree_by_simulation(t1[i], t1[j]);
					int deg = r.first;
					int new_red_edges_count = r.second;
					if (!dist_cel_mult_2) {
						continue;
					}
					if (!muchie_rosie_creata) {
						edge_pairs.clear();
						deg_vec.clear();
						red_edge.clear();
						edge_pairs.emplace_back(t1[i], t1[j]);
						deg_vec.emplace_back(deg);
						red_edge.emplace_back(number_of_red_edges + new_red_edges_count);
						i = t1.size();
						break;
					}
					edge_pairs.emplace_back(t1[i], t1[j]);
					deg_vec.emplace_back(deg);
					red_edge.emplace_back(number_of_red_edges + new_red_edges_count);

				}
			}

			for (int l = 0; l < edge_pairs.size(); ++l) {
				auto pairr = edge_pairs[l];
				pair <int, int > curr_hash = get_solution_hash(hash_id, step, pairr);
				if (!hash_to_index.count(curr_hash)) {
					partial_solution_hash new_hash;
					new_hash.max_red_degree = max(current_hash_solutions[step][hash_id].max_red_degree, deg_vec[l]);
					new_hash.number_of_red_edges = red_edge[l];
					new_hash.node_pairs = solution_hash.node_pairs;
					new_hash.node_pairs.emplace_back(pairr);
					hash_to_index[curr_hash] = current_hash_solutions[1 ^ step].size();
					current_hash_solutions[1 ^ step].emplace_back(new_hash);
				}
				else
				{
					if (max(deg_vec[l], current_hash_solutions[step][hash_id].max_red_degree) < current_hash_solutions[1 ^ step][hash_to_index[curr_hash]].max_red_degree) {
						current_hash_solutions[1 ^ step][hash_to_index[curr_hash]].max_red_degree = max(deg_vec[l], current_hash_solutions[step][hash_id].max_red_degree);
						current_hash_solutions[1 ^ step][hash_to_index[curr_hash]].number_of_red_edges = red_edge[l];
					}
				}
			}


			edge_pairs.clear();
			deg_vec.clear();
			red_edge.clear();
			t1.clear();
			for (auto it : all_nodes_set) {
				temp_edges[it].clear();
			}
			all_nodes_set.clear();
			if (getElapsed() > lower_bound_tl) {
				hash_to_index.clear();
				current_hash_solutions[step].clear();
				current_hash_solutions[1 ^ step].clear();
				break;
			}
		}


		if (getElapsed() > lower_bound_tl) {

			hash_to_index.clear();
			current_hash_solutions[step].clear();
			current_hash_solutions[1 ^ step].clear();
			break;
		}

		hash_to_index.clear();
		current_hash_solutions[step].clear();
	}
	current_hash_solutions[0].clear();
	current_hash_solutions[1].clear();
}


int get_lowerbound(int cc, int stop) {
	/* upperbound is required */
	ckeck_force_stop();
	int lower_bound = 0;
	vector<int> nodes_to_be_erased;
	vector<unordered_set<int>>edges2 = edges;
	for (auto it : best_solution) {
		nodes_to_be_erased.emplace_back(it.second);
	}
	if (!best_solution.empty()) {
		nodes_to_be_erased.emplace_back(best_solution.back().first);
	}
	for (int i = 0; i < nodes_to_be_erased.size(); ++i) {
		poss[nodes_to_be_erased[i]] = i;
	}
	auto init_nod = nodes_to_be_erased;
	for (int r = min(15, (int)nodes_to_be_erased.size()); r <= nodes_to_be_erased.size() && lower_bound != stop; ++r) {
		ckeck_force_stop();
		for (int tip = 1; tip <= 6 && lower_bound != stop; ++tip) {
			ckeck_force_stop();
			auto temp_vector = init_nod;
			if (tip == 2) {
				temp_vector.clear();
				deque <int> dq;
				dq.emplace_back(nodes_to_be_erased.back());
				set <int> nods;
				nods.insert(nodes_to_be_erased.back());
				temp_vector.emplace_back(nodes_to_be_erased.back());
				while (!dq.empty()) {
					int nod = dq.front();
					dq.pop_front();
					for (auto it : edges[nod]) {
						if (!nods.count(it)) {
							nods.insert(it);
							dq.emplace_back(it);
							temp_vector.emplace_back(it);
						}
					}
				}
				reverse(temp_vector.begin(), temp_vector.end());
			}
			if (tip == 3) {
				temp_vector.clear();
				deque <int> dq;
				dq.emplace_back(nodes_to_be_erased.back());
				set <int> nods;
				nods.insert(nodes_to_be_erased.back());
				temp_vector.emplace_back(nodes_to_be_erased.back());
				while (!dq.empty()) {
					int nod = dq.front();
					dq.pop_front();
					for (auto it : edges[nod]) {
						if (!nods.count(it)) {
							nods.insert(it);
							dq.emplace_back(it);
							temp_vector.emplace_back(it);
						}
					}
				}
			}
			if (tip == 4) {
				sort(temp_vector.begin(), temp_vector.end(), [](int i, int j) { return edges[i].size() < edges[j].size(); });
			}
			if (tip == 5) {
				temp_vector.clear();
				set <pair <int, int> > prior;
				prior.insert(make_pair(poss[nodes_to_be_erased.back()], nodes_to_be_erased.back()));
				set <int> nods;
				nods.insert(nodes_to_be_erased.back());
				temp_vector.emplace_back(nodes_to_be_erased.back());
				while (!prior.empty()) {
					int nod = (*prior.begin()).second;
					prior.erase(prior.begin());
					for (auto it : edges[nod]) {
						if (!nods.count(it)) {
							nods.insert(it);
							prior.insert(make_pair(poss[it], it));
							temp_vector.emplace_back(it);
						}
					}
				}
			}
			if (tip == 6) {
				temp_vector.clear();
				set <pair <int, int> > prior;
				prior.insert(make_pair(-poss[nodes_to_be_erased.back()], nodes_to_be_erased.back()));
				set <int> nods;
				nods.insert(nodes_to_be_erased.back());
				temp_vector.emplace_back(nodes_to_be_erased.back());
				while (!prior.empty()) {
					int nod = (*prior.begin()).second;
					prior.erase(prior.begin());
					for (auto it : edges[nod]) {
						if (!nods.count(it)) {
							nods.insert(it);
							prior.insert(make_pair(-poss[it], it));
							temp_vector.emplace_back(it);
						}
					}
				}
			}
			nodes_to_be_erased = temp_vector;
			int to_be_elimate = (int)nodes_to_be_erased.size() - r;
			edges = edges2;
			components[cc].clear();
			for (int j = to_be_elimate; j < nodes_to_be_erased.size(); ++j) {
				components[cc].emplace_back(nodes_to_be_erased[j]);
			}
			for (int j = 0; j < to_be_elimate; ++j) {
				for (auto it : edges[nodes_to_be_erased[j]]) {
					edges[it].erase(nodes_to_be_erased[j]);
				}
				edges[nodes_to_be_erased[j]].clear();
			}

			best_red_degree = 1e9;
			best_solution = vector<pair <int, int>>();
			if (getElapsed() > lower_bound_tl) {
				break;
			}
			auto rr = get_initial_upper_bound(cc);
			int local_upper_bound = rr.first;

			if (getElapsed() > lower_bound_tl) {
				break;
			}
			best_red_degree = local_upper_bound;
			best_solution = rr.second;
			max_clusters = (int)components[cc].size();
			branch_and_bound(cc);
			if (getElapsed() > lower_bound_tl) {
				break;
			}
			lower_bound = max(lower_bound, best_red_degree);

			best_red_degree = 1e9;
			best_solution = vector<pair <int, int>>();
			components[cc] = nodes_to_be_erased;
			edges = edges2;
			best_red_degree = 1e9;
			best_solution.clear();
		}
	}
	ckeck_force_stop();
	return lower_bound;
}



int eval_state(int cc, vector<pair <int, int>>sol, int lim) {
	for (auto i : components[cc]) {
		black[i] = edges[i];
		red[i] = unordered_set<int>();
	}
	for (int i = 0; i < sol.size(); ++i) {
		int r = merge_nodes(sol[i]);
		if (r > lim) {
			return i;
		}
	}
	return sol.size();
}




vector< pair <int, int> > go_heuristic(int cc, vector< pair <int, int > > sol, int max_pos) {
	vector< pair <int, int >> solution;

	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	for (auto i : components[cc]) {
		black[i] = edges[i];
		red[i] = unordered_set<int>();
		nodes2.insert(i);
		paired[i] = 0;
		paired2[i] = 0;
	}
	for (int i = 0; i < max_pos - 1 - rand() % 3; ++i) {
		merge_nodes(sol[i]);
		nodes2.erase(sol[i].second);
		solution.emplace_back(sol[i]);
	}
	for (auto it : nodes2) {
		nodes.emplace_back(it);
	}

	while (nodes.size() > 1) {
		pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > best_option = { { {1e9, 1e9}, {1e9, 1e9} } , {1e9, 1e9} };
		for (int i = 0; i < nodes.size(); ++i) {
			for (int j = i + 1; j < nodes.size(); ++j) {
				int x = nodes[i];
				int y = nodes[j];
				int red_between = 0;
				if (red[x].count(y)) {
					red_between = 1;
				}
				if (((int)min(red[x].size(), red[y].size())) - red_between + abs(((int)(black[x].size() + red[x].size())) - ((int)(black[y].size() + red[y].size()))) > best_option.first.first.first) {
					continue;
				}

				if (black[x].size() + red[x].size() > black[y].size() + red[y].size()) {
					swap(x, y);
				}
				int common_both_black = 0, common_mix = 0, common_both_red = 0, edges_between = 0, max_red_degree = 0, new_red_change = 0, new_red_not_change = 0;
				int cnt_max = black[x].size() + red[x].size();
				if (black[x].count(y) || red[x].count(y)) {
					edges_between = 2;
				}
				for (auto it : black[x]) {
					if (black[y].count(it)) {
						++common_both_black;
					}
					if (red[y].count(it)) {
						++common_mix;
					}
					if (it != y && !red[y].count(it) && !black[y].count(it))
					{
						max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : black[y]) {
					if (it != x && !red[x].count(it) && !black[x].count(it)) {
						max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
						new_red_change++;
					}
				}
				for (auto it : red[x]) {
					if (black[y].count(it)) {
						++common_mix;
					}
					if (red[y].count(it)) {
						++common_both_red;
					}
				}
				max_red_degree = max(max_red_degree, (int)(red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between));
				new_red_not_change = (int)(red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix - edges_between) - common_both_red;
				pair < pair < pair <int, int>, pair <int, int> >, pair <int, int> > curr_option = { { { max_red_degree, new_red_change }, {new_red_not_change, 0} }, { x, y } };
				best_option = min(best_option, curr_option);
			}
		}
		int x = best_option.second.first;
		int y = best_option.second.second;

		if (red[x].size() < red[y].size()) {
			swap(x, y);
		}
		red[x].erase(y);
		red[y].erase(x);
		black[y].erase(x);
		black[x].erase(y);
		for (auto it : red[y]) {
			red[it].erase(y);
			red[x].insert(it);
			red[it].insert(x);
		}
		for (auto it : black[x]) {
			if (black[y].count(it)) {
				temp.insert(it);
			}
			else
			{
				black[it].erase(x);
				red[x].insert(it);
				red[it].insert(x);
			}
		}
		for (auto it : black[y]) {
			if (!temp.count(it)) {
				red[x].insert(it);
				red[it].insert(x);
			}
			black[it].erase(y);
		}
		black[x] = temp;
		red[y] = unordered_set<int>();
		black[y] = unordered_set<int>();
		solution.emplace_back(x, y);
		nodes2.erase(y);
		temp = unordered_set<int>();
		nodes.clear();
		for (auto it : nodes2) {
			nodes.emplace_back(it);
		}
	}
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	return  solution;

}




pair <int, vector<pair <int, int>>> find_heuristic_upper_bound(int cc, vector<pair <int, int>> sol, int lim) {
	nodes_left = unordered_set<int>();
	int max_pos = eval_state(cc, sol, lim);
	if (max_pos == sol.size()) {
		return  make_pair(lim, sol);
	}
	int l;
	for (l = 1; l <= 500000 && max_pos != sol.size() && getElapsed() < lower_bound_tl; ++l) {
		vector<pair <int, int>> sol2 = sol;
		if (l % 3 == 1) {
			for (int p = 1; p <= 1; ++p) {
				for (auto i : components[cc]) {
					nodes_left.insert(i);
				}
				int pos_to_change = rand() % (max_pos + 1);
				for (int i = 0; i <= pos_to_change; ++i) {
					nodes_left.erase(sol2[i].second);
				}
				vector<int> lef;
				for (auto i : nodes_left) {
					lef.emplace_back(i);
				}
				int random_el = lef[rand() % lef.size()];
				nodes_left.clear();
				sol2[pos_to_change].first = random_el;
			}
		}
		if (l % 3 == 2)
		{
			for (int p = 1; p <= 1; ++p) {
				int pos_to_change = rand() % (max_pos + 1);
				nodes_left.insert(sol2[pos_to_change].second);
				for (int i = pos_to_change + 1; i < sol2.size(); ++i) {
					nodes_left.insert(sol2[i].first);
					nodes_left.insert(sol2[i].second);
				}
				vector<int> lef;
				for (auto i : nodes_left) {
					lef.emplace_back(i);
				}
				int rand_poz = rand() % lef.size();
				int random_el = lef[rand_poz];
				int val = sol2[pos_to_change].second;

				if (sol2[pos_to_change].first == random_el) {
					continue;
				}
				sol2[pos_to_change].second = random_el;
				for (int i = pos_to_change + 1; i < sol2.size(); ++i) {
					if (sol2[i].first == random_el) {
						sol2[i].first = val;
					}
					if (sol2[i].second == random_el) {
						sol2[i].second = val;
					}
				}
				nodes_left = unordered_set<int>();
			}
		}

		if (l % 3 == 0)
		{
			for (int p = 1; p <= 1; ++p) {

				for (auto i : components[cc]) {
					nodes_left.insert(i);
				}
				int pos_to_change = rand() % (max_pos + 1);
				for (int i = 0; i <= pos_to_change; ++i) {
					nodes_left.erase(sol2[i].second);
				}
				vector<int> lef;
				for (auto i : nodes_left) {
					lef.emplace_back(i);
				}
				int random_el = lef[rand() % lef.size()];
				nodes_left = unordered_set<int>();
				sol2[pos_to_change].first = random_el;

				/////////////
				pos_to_change = rand() % (max_pos + 1);
				nodes_left.insert(sol2[pos_to_change].second);
				for (int i = pos_to_change + 1; i < sol2.size(); ++i) {
					nodes_left.insert(sol2[i].first);
					nodes_left.insert(sol2[i].second);
				}
				lef = vector<int>();
				for (auto i : nodes_left) {
					lef.emplace_back(i);
				}
				int rand_poz = rand() % lef.size();
				random_el = lef[rand_poz];
				int val = sol2[pos_to_change].second;

				if (sol2[pos_to_change].first == random_el) {
					continue;
				}
				sol2[pos_to_change].second = random_el;
				for (int i = pos_to_change + 1; i < sol2.size(); ++i) {
					if (sol2[i].first == random_el) {
						sol2[i].first = val;
					}
					if (sol2[i].second == random_el) {
						sol2[i].second = val;
					}
				}
				nodes_left = unordered_set<int>();
			}
		}
		auto sol3 = go_heuristic(cc, sol2, max_pos);

		int new_pos = eval_state(cc, sol2, lim);
		int new_pos2 = eval_state(cc, sol3, lim);
		if (new_pos >= max_pos) {
			sol = sol2;
			max_pos = new_pos;
		}

		if (new_pos2 > max_pos) {
			sol = sol3;
			max_pos = new_pos2;
		}

	}
	//cout  << l << ' ' << max_pos * 100 << '\n';
	if (max_pos == (int) sol.size()) {
		return make_pair(lim, sol);

	}
	
	return make_pair(-1, sol);


}

void solver() {

	testNo = 26;
	readData();
	split_into_components();
	int first_comp_representant = 0;
	for (int cc = 0; cc < components.size(); ++cc) {
		LOWER_BOUND = -1;
		auto sol = reduce_component(cc);
		for (auto it : sol) {
			best_solution_all_graph.emplace_back(it);
		}
		best_red_degree = 1000000000;
		best_solution = vector< pair <int, int> >();
		auto upperbound = get_initial_upper_bound(cc);
		lower_bound_tl = 5 * 60;
		for (int best_degree = upperbound.first - 1; best_degree >= 0 && getElapsed() < lower_bound_tl; best_degree--) {
			auto new_sol = find_heuristic_upper_bound(cc, upperbound.second, upperbound.first - 1);
			if (new_sol.first != -1) {
				upperbound = new_sol;
			}
			else
			{
				break;
			}
		}
		auto best_red_degree_local = upperbound.first;
		auto best_solution_local = upperbound.second;
		best_solution = best_solution_local;
		lower_bound_tl = 10 * 60;
		LOWER_BOUND = get_lowerbound(cc, best_red_degree_local);
		best_solution = best_solution_local;
		best_red_degree = best_red_degree_local;
		if (LOWER_BOUND != best_red_degree) {
			if (LOWER_BOUND != best_red_degree) {
				lower_bound_tl = 99999999;
				max_clusters = components[cc].size();
				branch_and_bound(cc);
			}
		}
		for (auto it : best_solution) {
			best_solution_all_graph.emplace_back(it);
		}
		best_red_degree_all_graph = max(best_red_degree_all_graph, best_red_degree);

		if (cc > 0) {
			if (best_solution.empty()) {
				best_solution_all_graph.emplace_back(first_comp_representant, components[cc].back());

			}
			else
			{
				best_solution_all_graph.emplace_back(first_comp_representant, best_solution.back().first);
			}
		}
		else
		{
			if (best_solution.empty()) {
				first_comp_representant = components[cc].back();
			}
			else
			{
				first_comp_representant = best_solution.back().first;
			}
		}
	}
}




signed main() {
	srand(0);
	//ios_base::sync_with_stdio(false);
	//cin.tie();
	//std::cout.tie();
	begin_ = std::chrono::high_resolution_clock::now();
	solver();
	//cout << "Test no. " << testNo / 2 << '\n';
	//cout << "Max degree: " << best_red_degree_all_graph << '\n';
	//cout << "N: " << n << '\n';
	//cout << best_solution_all_graph.size() << '\n';
	for (auto it : best_solution_all_graph) {
		cout << it.first << ' ' << it.second << '\n';
	}
	//cout << "Time: " << getElapsed() << " s\n";
	return 0;
}