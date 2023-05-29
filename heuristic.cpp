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
#include <set>       

using namespace std;
int testNo, n, m;
const int SECONDS = 300;
const int MOD1 = 100000007;
const int PRIME1 = 10000019;
const int MOD2 = 100000037;
const int PRIME2 = 10000079;

std::chrono::time_point<std::chrono::high_resolution_clock> begin_;

string s;
vector<vector<int>>sorted_edges;
vector<unordered_set<int>>red, black;
vector<pair<int, int>> solution, best_solution,temp1, hashes,initial_solution;
vector<int> paired, paired2, nodes, ordered_nodes,bestDegree, temp_vector, single_nodes, position, position2;
unordered_set<int> neighbors, temp, nodes2, nodes_with_degree2, all_nodes;
vector < pair < int, pair < int, int > > > choices;
vector< vector< pair <int, int>> >node_hashes;
set < long long > nodes_set;
int best_red_degree = (int)1e9, number_of_edges, large_test;

double getElapsed() {
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin_);
	double sec = elapsed.count() * 1e-9;
	return sec;
}

void checkTime() {
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin_);
	double sec = elapsed.count() * 1e-9;

	if (sec >= SECONDS - 40 && initial_solution.size() + best_solution.size() != n - 1) {
		nodes.clear();
		for (auto it : nodes2) {
			nodes.emplace_back(it);
		}
		std::sort(nodes.begin(), nodes.end(), [](int i, int j) {return red[i].size() + black[i].size() < red[j].size() + black[j].size() || (red[i].size() + black[i].size() == red[j].size() + black[j].size() && red[i].size() < red[j].size());  });
		for (int i = 0; i + 1 < nodes.size(); ++i) {
			solution.emplace_back(nodes[i + 1], nodes[i]);
		}
		best_solution = solution;
		ios_base::sync_with_stdio(false);
		cout.tie();
		for (auto it : initial_solution) {
			cout << it.first << ' ' << it.second << '\n';
		}
		for (auto it : best_solution) {
			cout << it.first << ' ' << it.second << '\n';
		}
		exit(0);
	}

	if (sec >= SECONDS - 20) {
		ios_base::sync_with_stdio(false);
		cout.tie();
		for (auto it : initial_solution) {
			cout << it.first << ' ' << it.second << '\n';
		}
		for (auto it : best_solution) {
			cout << it.first << ' ' << it.second << '\n';
		}
		exit(0);
	}
}


void init_data() {
	sorted_edges.resize(n + 5, vector <int>());
	node_hashes.resize(n + 5, vector <pair <int, int >>());
	red.resize(n + 5, unordered_set <int>());
	black.resize(n + 5, unordered_set <int>());
	paired.resize(n + 5, 0);
	paired2.resize(n + 5, 0);
	position.resize(n + 5, 0);
	position2.resize(n + 5, 1e9);
	hashes.resize(n + 5, make_pair(PRIME1, PRIME2));
	for (int i = 1; i <= n; ++i) {
		hashes[i].first = 1ll * hashes[i - 1].first * PRIME1 % MOD1;
		hashes[i].second = 1ll * hashes[i - 1].second * PRIME2 % MOD2;
	}
}

void clear_data() {
	red.clear();
	black.clear();
	paired.clear();
	paired2.clear();
	nodes.clear();
	ordered_nodes.clear();
	neighbors.clear();
	choices.clear();
	nodes_set.clear();
}


inline int eval(int x, int y) {
	if (x == y) {
		return 1000000000;
	}
	if (black[x].size() + red[x].size() > black[y].size() + red[y].size()) {
		swap(x, y);
	}
	int common_both_black = 0, common_mix = 0, common_both_red = 0;
	int cnt_max = black[x].size() + red[x].size(), cnt1 = 0, cnt2 = 0;

	for (auto it : black[x]) {
		checkTime();
		if (cnt1 == cnt_max) {
			break;
		}
		if (black[y].count(it)) {
			++common_both_black;
		}
		if (red[y].count(it)) {
			++common_mix;
		}
		++cnt1;
	}
	for (auto it : red[x]) {
		checkTime();
		if (cnt2 == cnt_max) {
			break;
		}
		if (black[y].count(it)) {
			++common_mix;
		}
		if (red[y].count(it)) {
			++common_both_red;
		}
		++cnt2;
	}
	return red[x].size() + black[x].size() + red[y].size() + black[y].size() - 2 * common_both_black - common_both_red - common_mix;
}
vector<int> singletons;
int solve2(int max_d, int max_edges) {

	int ones = (1 << 27) - 1;
	int max_degree = 0, round = 0;
	int force_stop = 0;
	number_of_edges = m;

	for (auto i : all_nodes) {
		paired[i] = 0;
		paired2[i] = 0;
		nodes2.insert(i);
		nodes_set.insert(1ll*((1ll * black[i].size() << 28) + i));
	}
	while (!nodes_set.empty() && number_of_edges > max_edges) {
		checkTime();
		int x = 0;
		x = 1ll * (*nodes_set.begin()) & ones;

		if (  red[x].size() + black[x].size() == 0 )	{
			nodes_set.erase( 1ll * ( (1ll * (red[x].size() + black[x].size()) << 28) + 1ll * x));
			singletons.emplace_back(x);
			continue;
		}


		int best_degree = (int)1e9, best_node = 0, worst_degree = 0, worst_node = 0;
		for (auto it : black[x]) {
			checkTime();
			if (black[it].size() + red[it].size() < best_degree) {
				best_degree = (int)(black[it].size() + red[it].size());
				best_node = it;
			}
			if (black[it].size() + red[it].size() > worst_degree) {
				worst_degree = (int)(black[it].size() + red[it].size());
				worst_node = it;
			}
		}
		for (auto it : red[x]) {
			checkTime();
			if (black[it].size() + red[it].size() < best_degree) {
				best_degree = (int)(black[it].size() + red[it].size());
				best_node = it;
			}
			if (black[it].size() + red[it].size() > worst_degree) {
				worst_degree = (int)(black[it].size() + red[it].size());
				worst_node = it;
			}
		}
		if (best_node == 0) {
			auto it = nodes_set.begin();
			++it;
			best_node = 1ll * (*it) & ones;
			best_degree = (1ll * (*it) >> 28);
		}

		int best_nbr_score = eval(x, best_node), best_nbr = best_node;

		if (worst_node != 0) {
			int cnt = 0;
			for (auto it : black[worst_node]) {
				if (cnt == 10) {
					break;
				}
				if (it == x) {
					continue;
				}
				int ll = eval(x, it);
				if (ll < best_nbr_score) {
					best_nbr_score = ll;
					best_nbr = it;
				}
				++cnt;
			}
			for (auto it : red[worst_node]) {
				if (cnt == 10) {
					break;
				}
				if (it == x) {
					continue;
				}
				int ll = eval(x, it);
				if (ll < best_nbr_score) {
					best_nbr_score = ll;
					best_nbr = it;
				}
				++cnt;
			}
		}
		best_node = best_nbr;


		int y = best_node;
		nodes_set.erase( 1ll * ( (1ll * (red[best_node].size() + black[best_node].size()) << 28) + 1ll * best_node));
		nodes_set.erase(1ll * (1ll * (red[x].size() + black[x].size()) << 28) +1ll *  x);
		if (red[x].size() < red[y].size()) {
			swap(x, y);
		}
		nodes2.erase(y);
		if (red[x].count(y) || black[x].count(y)) {
			--number_of_edges;
		}
		red[x].erase(y);
		red[y].erase(x);
		black[y].erase(x);
		black[x].erase(y);
		for (auto it : red[y]) {
			if (red[x].count(it) || black[x].count(it)) {
				--number_of_edges;
			}
			unsigned val = red[it].size() + black[it].size();
			red[it].erase(y);
			red[x].insert(it);
			red[it].insert(x);
			if (red[it].size() + black[it].size() != val) {
				nodes_set.erase( 1ll * ( (1ll * val << 28) + it));
				nodes_set.insert( 1ll * ( (1ll * (red[it].size() + black[it].size()) << 28) + 1ll * it));
			}
			max_degree = max(max_degree, (int)red[it].size());
		}
		for (auto it : black[x]) {
			if (black[y].count(it)) {
				--number_of_edges;
				temp.insert(it);
			}
			else
			{
				unsigned val = red[it].size() + black[it].size();
				black[it].erase(x);
				red[x].insert(it);
				red[it].insert(x);
				if (red[it].size() + black[it].size() != val) {
				nodes_set.erase( 1ll * ( (1ll * val << 28) + it));
				nodes_set.insert( 1ll * ( (1ll * (red[it].size() + black[it].size()) << 28) +1ll *  it));
				}
				max_degree = max(max_degree, (int)red[it].size());
			}
		}
		for (auto it : black[y]) {
			if (red[x].count(it)) {
				--number_of_edges;
			}
			unsigned val = red[it].size() + black[it].size();
			if (!temp.count(it)) {
				red[x].insert(it);
				red[it].insert(x);
				max_degree = max(max_degree, (int)red[it].size());
			}
			black[it].erase(y);
			if (red[it].size() + black[it].size() != val) {
				nodes_set.erase( 1ll * ( (1ll * val << 28) + it));
				nodes_set.insert( 1ll * ( (1ll * (red[it].size() + black[it].size()) << 28) + 1ll * it));
			}
		}
		black[x] = temp;
		red[y] = unordered_set<int>();
		black[y] = unordered_set<int>();
		nodes_set.insert((1ll * (red[x].size() + black[x].size()) << 28) + 1ll * x);
		max_degree = max(max_degree, (int)red[x].size());
		solution.emplace_back(x, y);
		temp = unordered_set<int>();
		if (max_degree > max_d) {
			force_stop = 1;
			break;
		}
	}
	for ( int i = 1 ; i < singletons.size(); ++ i )	{
		solution.emplace_back(singletons[0],  singletons[i]);
		nodes2.erase(singletons[i]);
	}
	singletons = vector<int>();
	if (max_degree < best_red_degree && !force_stop && nodes_set.size() == 0) {
		best_red_degree = max_degree;
		best_solution = solution;
	}
	return max_degree;
}

int lgp(int a, int b, int m) {
	if (!b) {
		return 1;
	}
	if (b % 2) {
		return 1ll * a * lgp(1ll * a * a % m, b / 2, m) % m;
	}
	return lgp(1ll * a * a % m, b / 2, m);
}

void readData() {
	ios_base::sync_with_stdio(false);
	cin.tie();
	string problem_id;
	cin >> problem_id >> problem_id >> n >> m;
	init_data();
	for (int j = 1; j <= m; ++j) {
		int x = 0, y = 0;
		cin >> x >> y;
		black[x].insert(y);
		black[y].insert(x);
		sorted_edges[x].emplace_back(y);
		sorted_edges[y].emplace_back(x);
	}

	for (int i = 1; i <= n; ++i) {
		sort(sorted_edges[i].begin(), sorted_edges[i].end());
		node_hashes[i].resize(sorted_edges[i].size());
		pair <int, int> ret = make_pair(0, 0);
		for (int j = 0; j < sorted_edges[i].size(); ++j) {
			ret.first = 1ll * (1ll * ret.first + 1ll * sorted_edges[i][j] * hashes[j].first) % MOD1;
			ret.second = 1ll * (1ll * ret.second + 1ll * sorted_edges[i][j] * hashes[j].second) % MOD2;
			node_hashes[i][j] = ret;
		}
	}
	if (n == 0) {
		exit(0);
	}
}

void solve1(int t, int its) {
	int max_degree = 0, round = 0, mid_value = 0;
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	solution.clear();
	for (auto i : all_nodes) {
		if (large_test) {
			black[i] = unordered_set<int>();
			for (auto it : sorted_edges[i]) {
				black[i].insert(it);
			}
			for (auto it : sorted_edges[i]) {
				black[i].insert(it);
			}	red[i] = unordered_set<int>();
		}
		nodes.push_back(i);
		nodes2.insert(i);
		paired[i] = 0;
		paired2[i] = 0;
	}
	large_test = 1;
	checkTime();
	ordered_nodes = vector<int>();
	int force_stop = 0, iter = 0;

	while (nodes.size() > 1) {
		++iter;
		checkTime();
		++round;
		int size = (int)nodes.size();
		std::sort(nodes.begin(), nodes.end(), [](int i, int j) {return red[i].size() + black[i].size() > red[j].size() + black[j].size(); });
		for (int i = nodes.size() - 1; i - 1 >= 0 && red[nodes[i]].size() + black[nodes[i]].size() == 0 && red[nodes[i - 1]].size() + black[nodes[i - 1]].size() == 0; i -= 2) {
			choices.emplace_back(make_pair(0, make_pair(nodes[i], nodes[i - 1])));
		}
		checkTime();
		mid_value = red[nodes[max(0, ((int)(nodes.size()) - 10000))]].size() + black[nodes[max(0, ((int)(nodes.size()) - 10000))]].size();

		for (auto it : nodes) {
			checkTime();
			paired[it] = 0;
			paired2[it] = 0;
		}
		for (auto node : nodes) {
			checkTime();
			if (paired[node]) {
				continue;
			}
			temp_vector.push_back(node);
			paired[node] = round;
			for (auto it : black[node]) {
				checkTime();
				if (!paired[it]) {
					temp_vector.emplace_back(it);
					paired[it] = round;
				}
			}
			for (auto it : red[node]) {
				checkTime();
				if (paired[it] != 0 && paired[it] != round) {
					exit(1);
				}
				if (!paired[it]) {
					temp_vector.emplace_back(it);
					paired[it] = round;
				}
			}
			sort(temp_vector.begin(), temp_vector.end(), [](const int i, const int j) {return red[i].size() + black[i].size() > red[j].size() + black[j].size() || (red[i].size() + black[i].size() == red[j].size() + black[j].size() && red[i].size() > red[j].size()); });
			for (auto it : temp_vector) {
				checkTime();
				ordered_nodes.emplace_back(it);
			}
			temp_vector = vector<int>();
		}

		for (int i = 0; i < ordered_nodes.size() && choices.size() < 50000000; ++i) {
			checkTime();
			if (red[ordered_nodes[i]].size() + black[ordered_nodes[i]].size() > mid_value) {
				continue;
			}
			for (int l = 1; l <= its && i - l >= 0 && choices.size() < 50000000; ++l) {
				checkTime();
				int x = ordered_nodes[i];
				int y = ordered_nodes[i - l];
				if (x == y) {
					continue;
				}
				if (black[x].size() + red[x].size() > black[y].size() + red[y].size()) {
					swap(x, y);
				}
				int common_both_black = 0, common_mix = 0, common_both_red = 0, edges_between = 0, max_red_degree = 0;
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
					}
				}
				for (auto it : black[y]) {
					if (it != x && !red[x].count(it) && !black[x].count(it)) {
						max_red_degree = max(max_red_degree, (int)(red[it].size() + 1));
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
				choices.emplace_back(make_pair(max_red_degree, make_pair(x, y)));
			}
			if (choices.size() >= 50000000) {
				break;
			}
		}


		checkTime();
		std::sort(choices.begin(), choices.end());
		nodes = vector<int>();
		int made = 0;
		checkTime();
		for (int i = 0; i < choices.size(); ++i) {
			checkTime();
			int x, y;
			x = choices[i].second.first;
			y = choices[i].second.second;

			checkTime();
			temp = unordered_set<int>();
			if (made < max(1, (int)size / (t + 1)) && !paired2[x] && !paired2[y]) {
				++made;
				paired2[x] = 1;
				paired2[y] = 1;
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
				max_degree = max(max_degree, (int)red[x].size());
				solution.emplace_back(x, y);
				nodes2.erase(y);
				temp = unordered_set<int>();
				checkTime();
			}
		}
		for (auto it : nodes2) {
			nodes.emplace_back(it);
		}
		temp.clear();


		ordered_nodes = vector<int>();
		choices = vector < pair < int, pair <int, int> > >();

		checkTime();
		if (bestDegree.size() + 1 == iter) {
			bestDegree.emplace_back(max_degree);
		}
		else
		{
			if (max_degree > bestDegree[iter - 1]) {
				force_stop = 1;
				break;
			}
		}
		if (max_degree >= best_red_degree) {
			force_stop = 1;
			break;
		}

		if (nodes.size() - 1 < max_degree) {
			while (nodes.size() > 1) {
				solution.emplace_back(nodes[0], nodes.back());
				nodes.pop_back();
			}
		}

	}
	for (auto i : all_nodes) {
		black[i] = unordered_set<int>();
		red[i] = unordered_set<int>();
	}
	checkTime();
	nodes = vector<int>();
	nodes2 = unordered_set<int>();
	if (max_degree < best_red_degree && !force_stop) {
		best_red_degree = max_degree;
		best_solution = solution;
	}
	solution.clear();
	nodes_set = set < long long >();
	for (auto i : all_nodes) {
		paired[i] = 0;
		paired2[i] = 0;
	}
}


vector<pair < pair <int, int>, int>> nodes_hashes;
vector<int> sorted_nbr;

void eliminate_duplicates() {

	for (int i = 1; i <= n; ++i) {
		all_nodes.insert(i);
		if ( black[i].empty() )	{
			nodes_hashes.emplace_back( make_pair(0, 0) , i);
		}	else{
			nodes_hashes.emplace_back(node_hashes[i].back(), i);
		}
	}



	sort(nodes_hashes.begin(), nodes_hashes.end());

	int cnt = 0, maxi = 0;
	for (int i = 0; i < nodes_hashes.size();) {
		int j = i + 1;
		for (; j < nodes_hashes.size() && nodes_hashes[i].first == nodes_hashes[j].first; ++j) {
			++cnt;
			auto it = make_pair(nodes_hashes[i].second, nodes_hashes[j].second);
			initial_solution.emplace_back(it);
			for (auto it2 : black[it.second]) {
				black[it2].erase(it.second);
			}
			black[it.second] = unordered_set<int>();
			all_nodes.erase(it.second);
		}
		i = j;
	}

	int cnt2 = 0;
	for (int i = 1; i <= n; ++i) {
		if (!all_nodes.count(i)) {
			continue;
		}
		for (auto it : sorted_edges[i]) {
			if (it < i || !all_nodes.count(it)) {
				continue;
			}
			int x = i;
			int y = it;
			if (black[x].size() != black[y].size()) {
				continue;
			}
			// position of x in vector y
			int st = 1, dr = sorted_edges[y].size(), mid, poz_x_in_y = -1, poz_y_in_x = -1;
			while (st <= dr) {
				mid = (st + dr) / 2;
				if (sorted_edges[y][mid - 1] <= x) {
					st = mid + 1;
					if (sorted_edges[y][mid - 1] == x) {
						poz_x_in_y = mid - 1;
					}
				}
				else
				{
					dr = mid - 1;
				}
			}
			st = 1;
			dr = sorted_edges[x].size();
			// position of y in vector x
			while (st <= dr) {
				mid = (st + dr) / 2;
				if (sorted_edges[x][mid - 1] <= y) {
					st = mid + 1;
					if (sorted_edges[x][mid - 1] == y) {
						poz_y_in_x = mid - 1;
					}
				}
				else
				{
					dr = mid - 1;
				}
			}
			if (poz_y_in_x == -1) {
				exit(1);
			}

			// get hash of x list without y
			pair <int, int > prev;
			prev = make_pair(0, 0);
			if (poz_y_in_x != 0) {
				prev = node_hashes[x][poz_y_in_x - 1];
			}
			pair <int, int> res_x;
			res_x.first = (1ll * (1ll * (1ll * node_hashes[x].back().first - 1ll * node_hashes[x][poz_y_in_x].first + 1ll * MOD1) % MOD1) * lgp(PRIME1, MOD1 - 2, MOD1) % MOD1 + 1ll * prev.first + 1ll * MOD1) % MOD1;
			res_x.second = (1ll * (1ll * (1ll * node_hashes[x].back().second - 1ll * node_hashes[x][poz_y_in_x].second + 1ll * MOD2) % MOD2) * lgp(PRIME2, MOD2 - 2, MOD2) % MOD2 + 1ll * prev.second + 1ll * MOD2) % MOD2;
			// get hash of y list without x
			prev = make_pair(0, 0);
			if (poz_x_in_y != 0) {
				prev = node_hashes[y][poz_x_in_y - 1];
			}
			pair <int, int> res_y;
			res_y.first = (1ll * (1ll * (1ll * node_hashes[y].back().first - 1ll * node_hashes[y][poz_x_in_y].first + 1ll * MOD1) % MOD1) * lgp(PRIME1, MOD1 - 2, MOD1) % MOD1 + 1ll * prev.first + 1ll * MOD1) % MOD1;
			res_y.second = (1ll * (1ll * (1ll * node_hashes[y].back().second - 1ll * node_hashes[y][poz_x_in_y].second + 1ll * MOD2) % MOD2) * lgp(PRIME2, MOD2 - 2, MOD2) % MOD2 + 1ll * prev.second + 1ll * MOD2) % MOD2;

			if (res_x != res_y) {
				continue;
			}

			++cnt2;
			for (auto it2 : black[y]) {
				black[it2].erase(y);
			}
			black[y] = unordered_set<int>();
			all_nodes.erase(y);
			initial_solution.emplace_back(x, y);

		}
	}
	nodes_hashes = vector<pair < pair <int, int>, int>>();

	for (auto it : all_nodes) {
		sorted_edges[it].clear();
		for (auto it2 : black[it]) {
			sorted_edges[it].emplace_back(it2);
		}
	}

}

void solver() {
	readData();
	int max_deg = 0;
	for (int i = 1; i <= n; ++i) {
		max_deg = max(max_deg, (int)black[i].size());
	}
	eliminate_duplicates();
	large_test = 0;
	if (max_deg < 2000) {
		large_test = 1;
		for (int t = 0; t <= 0; ++t) {
			checkTime();
			solve2(m + 1, 0);
			solution.clear();
			nodes_set = set < long long>();
		}
	}	
	for (int t = 7; ; t += 3) {
		checkTime();
		solve1(t, t + 2);
	}
}

signed main() {
	begin_ = std::chrono::high_resolution_clock::now();
	solver();
	return 0;
}