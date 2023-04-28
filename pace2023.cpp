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


std::chrono::time_point<std::chrono::high_resolution_clock> begin_;

vector<unordered_set<int>>edges;
vector<unordered_set<int>>red;
vector<unordered_set<int>>black;
vector<int> paired;
string s;
vector<pair<int, int>> solution, best_solution;
vector<int> nodes, ordered_nodes;
unordered_set<int> neighbors, temp;
vector < pair < int, pair < int, int > > > choices;
set<pair <int, int > > nodes_set;
vector<int> bestDegree;
int best_red_degree = 1e9;
vector<int>temp_vector;


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
	if (sec >= SECONDS - 30) {
		//	cout << best_red_degree << ' ' << best_solution.size() << '\n';
		//	cout << sec << '\n';
		for (auto it : best_solution) {
			cout << it.first << ' ' << it.second << '\n';
		}
		exit(0);
	}
}


void init_data() {
	edges.resize(n + 1, unordered_set <int>());
	red.resize(n + 1, unordered_set <int>());
	black.resize(n + 1, unordered_set <int>());
	paired.resize(n + 1, 0);
}

void clear_data() {
	edges.clear();
	red.clear();
	black.clear();
	paired.clear();
	nodes.clear();
	ordered_nodes.clear();
	neighbors.clear();
	choices.clear();
	nodes_set.clear();
}

void readData() {
	string path_input = R"(C:\Users\andre\OneDrive\Desktop\teste\heuristic-public\heuristic_)";
	for (int l = 1; l <= 3 - to_string(testNo).size(); ++l) {
		path_input.push_back('0');
	}
	path_input += to_string(testNo);
	path_input += ".gr";
	//ifstream in(path_input);
	cin.tie();
	int t;
	string problem_id;
	cin >> problem_id >> problem_id >> n >> m;
	init_data();
	/*
	while (getline(in, s)) {
		deque<char> input;
		for (auto it : s) {
			input.push_back(it);
		}
		while (!input.empty() && input.front() == ' ') {
			input.pop_front();
		}
		s.clear();
		for (char i : input) {
			s.push_back(i);
		}
		if (s.empty() || s[0] != 'p') continue;
		istringstream is(s);
		string problem_id;
		is >> problem_id >> problem_id >> n >> m;
		init_data();
		break;
	}
	*/
	for (int j = 1; j <= m; ++j) {
		int x = 0, y = 0;

		cin >> x >> y;
		//cout << x << ' ' << y << '\n';
		//in.close();
		//exit(0);
		/*
		getline(in, s);
		deque<char> input;
		for (auto it : s) {
			input.push_back(it);
		}
		while (!input.empty() && input.front() == ' ') {
			input.pop_front();
		}
		s.clear();
		if (!input.empty() && input.front() == 'c') {
			j--;
			continue;
		}
		while ('0' <= input.front() && input.front() <= '9') {
			x = x * 10 + input.front() - '0';
			input.pop_front();
		}
		while (!('0' <= input.front() && input.front() <= '9')) {
			input.pop_front();
		}
		while (!input.empty() && '0' <= input.front() && input.front() <= '9') {
			y = y * 10 + input.front() - '0';
			input.pop_front();
		}
		*/
		edges[x].insert(y);
		edges[y].insert(x);
	}
	//in.close();
	if (n == 0) {
		exit(0);
	}
}


int solve2(int max_d) {

	for (int i = 1; i <= n; ++i) {
		black[i] = edges[i];
		red[i] = unordered_set<int>();
		paired[i] = 0;
		nodes_set.insert(make_pair(red[i].size() + black[i].size(), i));
	}
	int max_degree = 0, round = 0;
	int force_stop = 0;
	while (nodes_set.size() > 1 && max_degree < max_d) {
		checkTime();
		//if (nodes_set.size() % 10000 == 1) {
			//cout << nodes_set.size() << ' ' << max_degree << '\n';
		//}
		int x = 0;
		int best_degree, best_node = 0, second_best_node = 0, second_best_degree;
		/*cel mai mic nod cu cel mai mic vecin sau urmatorul cel mai mic nod*/
		x = (*nodes_set.begin()).second;
		best_degree = 1e9, best_node = 0;
		for (auto it : black[x]) {
			if (black[it].size() + red[it].size() < best_degree) {
				best_degree = black[it].size() + red[it].size();
				best_node = it;
			}
		}
		for (auto it : red[x]) {
			if (black[it].size() + red[it].size() < best_degree) {
				best_degree = black[it].size() + red[it].size();
				best_node = it;
			}
		}
		if (best_node == 0) {
			auto it = nodes_set.begin();
			++it;
			best_node = (*it).second;
			best_degree = (*it).first;

		}
		if (second_best_node != 0) {
			x = second_best_node;
		}
		int y = best_node;
		nodes_set.erase(make_pair(red[best_node].size() + black[best_node].size(), best_node));
		nodes_set.erase(make_pair(red[x].size() + black[x].size(), x));
		if (red[x].size() < red[y].size()) {
			swap(x, y);
		}
		red[x].erase(y);
		red[y].erase(x);
		black[y].erase(x);
		black[x].erase(y);
		for (auto it : red[y]) {
			nodes_set.erase(make_pair(red[it].size() + black[it].size(), it));
			red[it].erase(y);
			red[x].insert(it);
			red[it].insert(x);
			nodes_set.insert(make_pair(red[it].size() + black[it].size(), it));
			max_degree = max(max_degree, (int)red[it].size());
		}
		for (auto it : black[x]) {
			if (black[y].count(it)) {
				temp.insert(it);
			}
			else
			{
				nodes_set.erase(make_pair(red[it].size() + black[it].size(), it));
				black[it].erase(x);
				red[x].insert(it);
				red[it].insert(x);
				nodes_set.insert(make_pair(red[it].size() + black[it].size(), it));
				max_degree = max(max_degree, (int)red[it].size());
			}
		}
		for (auto it : black[y]) {
			nodes_set.erase(make_pair(red[it].size() + black[it].size(), it));
			if (!temp.count(it)) {
				red[x].insert(it);
				red[it].insert(x);
				max_degree = max(max_degree, (int)red[it].size());
			}
			black[it].erase(y);
			nodes_set.insert(make_pair(red[it].size() + black[it].size(), it));
		}
		black[x] = temp;
		red[y].clear();
		black[y].clear();
		nodes_set.insert(make_pair(red[x].size() + black[x].size(), x));
		max_degree = max(max_degree, (int)red[x].size());
		solution.emplace_back(x, y);
		temp.clear();
		if (max_degree > max_d) {
			force_stop = 1;
			break;
		}
	}
	if (max_degree < best_red_degree && !force_stop && nodes_set.size() == 1) {
		best_red_degree = max_degree;
		best_solution = solution;
	}
	return max_degree;
}

void solve1(int t) {
	int max_degree = 0, round = 0;
	nodes.clear();
	solution.clear();
	if (t % 2 == 1) {
		for (int i = 1; i <= n; ++i) {
			black[i] = edges[i];
			red[i].clear();
			nodes.push_back(i);
			paired[i] = 0;
		}
	}
	else
	{
		max_degree = max(solve2(50), max_degree);
		for (auto it : nodes_set) {
			nodes.emplace_back(it.second);
		}
	}
	ordered_nodes.clear();
	int force_stop = 0, iter = 0;
	while (nodes.size() > 1) {
		++iter;
		checkTime();
		++round;
		//unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
		sort(nodes.begin(), nodes.end(), [](int i, int j) {return red[i].size() + black[i].size() > red[j].size() + black[j].size(); });
		for (auto it : nodes) {
			paired[it] = 0;
		}
		for (auto node : nodes) {
			if (paired[node]) {
				continue;
			}
			ordered_nodes.push_back(node);
			paired[node] = round;
			for (auto it : black[node]) {
				if (!paired[it]) {
					temp_vector.emplace_back(it);
					paired[it] = round;
				}
			}
			for (auto it : red[node]) {
				if (paired[it] != 0 && paired[it] != round) {
					cout << "WTF\n";
				}
				if (!paired[it]) {
					temp_vector.emplace_back(it);
					paired[it] = round;
				}
			}
			//if (t % 2) {
				//shuffle(temp_vector.begin(), temp_vector.end(), std::default_random_engine(t));
			sort(temp_vector.begin(), temp_vector.end(), [](const int i, const int j) { return black[i].size() + red[i].size() > black[j].size() + red[j].size(); });
			//}
			//else
			//{
				//sort(temp_vector.begin(), temp_vector.end(), [](const int i, const int j) { return black[i].size() + red[i].size() < black[j].size() + red[j].size(); });
			//}
			for (auto it : temp_vector) {
				ordered_nodes.emplace_back(it);
			}
			temp_vector.clear();
		}
		checkTime();
		for (int i = 0; i + 1 < ordered_nodes.size(); i += 2) {
			int x = ordered_nodes[i];
			int y = ordered_nodes[i + 1];
			int black_common = 0;
			for (auto it : black[x]) {
				neighbors.insert(it);
				if (black[y].count(it)) {
					black_common++;
				}
			}
			for (auto it : black[y]) {
				neighbors.insert(it);
			}
			for (auto it : red[x]) {
				neighbors.insert(it);
			}
			for (auto it : red[y]) {
				neighbors.insert(it);
			}
			choices.emplace_back(make_pair(neighbors.size() - black_common, make_pair(x, y)));
			neighbors.clear();
		}
		sort(choices.begin(), choices.end());
		int act = -1;
		if (ordered_nodes.size() % 2) {
			act = ordered_nodes.back();
		}
		nodes.clear();
		if (act != -1) {
			nodes.emplace_back(act);
		}
		for (int i = 0; i < choices.size(); ++i) {
			checkTime();
			int x, y;
			x = choices[i].second.first;
			y = choices[i].second.second;

			temp.clear();
			if (i < max(1, (int)choices.size() / (t + 1))) {
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
				red[y].clear();
				black[y].clear();
				max_degree = max(max_degree, (int)red[x].size());
				solution.emplace_back(x, y);
				nodes.emplace_back(x);
				temp.clear();
			}
			else
			{
				nodes.emplace_back(x);
				nodes.emplace_back(y);
			}
		}
		ordered_nodes.clear();
		choices.clear();
		//cout << "Round: " << round << " Max degree: " << max_degree << '\n';
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
		if (max_degree > best_red_degree) {
			force_stop = 1;
			break;
		}
	}
	for (int i = 1; i <= n; ++i) {
		black[i].clear();
		red[i].clear();
	}
	nodes.clear();
	if (max_degree < best_red_degree && !force_stop) {
		best_red_degree = max_degree;
		best_solution = solution;
	}
	solution.clear();
	nodes_set.clear();
	for (int i = 1; i <= n; ++i) {
		paired[i] = 0;
	}
}


void solver() {
	int  cnt = 0;
	int maxN = 0;
	for (testNo = 176; testNo <= 176; testNo += 2) {
		readData();

		//auto end = std::chrono::high_resolution_clock::now();
		//auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin_);
		//double sec = elapsed.count() * 1e-9;
		//cout << sec << '\n';

		//cout << testNo << ' ' << n << ' ' << m << '\n';
		maxN = max(maxN, n);
		for (int t = 0; t <= 0; ++t) {
			solve2(m + 1);
			solution.clear();
			nodes_set.clear();
			for (int i = 1; i <= n; ++i) {
				black[i].clear();
				red[i].clear();
			}
		}
		for (int t = 1; ; t += 2) {
			solve1(t);
		}
		//clear_data();
	}
}


signed main() {
	ios_base::sync_with_stdio(false);
	cin.tie();
	cout.tie();
	begin_ = std::chrono::high_resolution_clock::now();
	solver();
	return 0;
}

int check() {
	for (int i = 1; i <= n; ++i) {
		black[i] = edges[i];
		red[i].clear();
	}
	int ret = 0;
	vector<int> ok(n + 1, 1);
	for (auto it : solution) {
		//if (!ok[it.first] || !ok[it.second]) {
			//cout << "HEY!\n";
			//exit(0);
		//}
		ok[it.second] = 0;
		int x = it.first;
		int y = it.second;
		red[x].erase(y);
		red[y].erase(x);
		black[y].erase(x);
		black[x].erase(y);
		for (auto it : red[y]) {
			red[it].erase(y);
			red[x].insert(it);
			red[it].insert(x);
			ret = max(ret, (int)red[it].size());
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
				ret = max(ret, (int)red[it].size());
			}
		}
		for (auto it : black[y]) {
			if (!temp.count(it)) {
				red[x].insert(it);
				red[it].insert(x);
				ret = max(ret, (int)red[it].size());
			}
			black[it].erase(y);
		}
		black[x] = temp;
		red[y].clear();
		black[y].clear();
		ret = max(ret, (int)red[x].size());
		temp.clear();
	}
	return ret;
}