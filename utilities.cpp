#include "utilities.h"

double getLargest(vector<double> x) {
	if (x.empty()) return BIGVALUE;
	double temp = x[0];
	for (int i = 1; i < (int)x.size(); i++) {
		if (temp < x[i]) {
			temp = x[i];
		}
	}
	return temp;
}

double getSmallest(vector<double> x) {
	if (x.empty()) return -1;
	double temp = x[0];
	for (int i = 1; i < (int)x.size(); i++) {
		if (temp > x[i]) {
			temp = x[i];
		}
	}
	return temp;
}

vector<int> convexHull(vector<int> nodesemp, Case* instance) {
	vector<int> nodes = nodesemp;
	vector<int> circle;
	int lowestn = 0;
	for (int i = 0; i < (int)nodes.size(); i++) {
		if (instance->positions[nodes[i]].second < instance->positions[nodes[lowestn]].second ||
			(instance->positions[nodes[i]].second == instance->positions[nodes[lowestn]].second &&
				instance->positions[nodes[i]].first < instance->positions[nodes[lowestn]].first)) {
			lowestn = i;
		}
	}
	int startNode = nodes[lowestn];
	circle.push_back(startNode);
	nodes.erase(nodes.begin() + lowestn);
	vector<double> angles;
	for (int i = 0; i < (int)nodes.size(); i++) {
		pair<double, double> vec = make_pair(instance->positions[nodes[i]].first - instance->positions[startNode].first,
			instance->positions[nodes[i]].second - instance->positions[startNode].second);
		double onecos = vec.first / sqrt(vec.first * vec.first + vec.second * vec.second);
		angles.push_back(onecos);
	}

	while (!nodes.empty())
	{
		int thebig = 0;
		for (int i = 0; i < (int)angles.size(); i++) {
			if (angles[i] > angles[thebig]) {
				thebig = i;
			}
		}
		int thebignode = nodes[thebig];
		if (circle.size() <= 2) {
			circle.push_back(thebignode);
			nodes.erase(nodes.begin() + thebig);
			angles.erase(angles.begin() + thebig);
		}
		else {
			while (circle.size() > 2) {
				int lastnode = circle[circle.size() - 1];
				int llastnode = circle[circle.size() - 2];
				pair<double, double> lastvec = make_pair(instance->positions[lastnode].first - instance->positions[llastnode].first,
					instance->positions[lastnode].second - instance->positions[llastnode].second);
				pair<double, double> nowvec = make_pair(instance->positions[thebignode].first - instance->positions[lastnode].first,
					instance->positions[thebignode].second - instance->positions[lastnode].second);
				if (lastvec.first * nowvec.second - lastvec.second * nowvec.first >= 0) {
					circle.push_back(thebignode);
					nodes.erase(nodes.begin() + thebig);
					angles.erase(angles.begin() + thebig);
					break;
				}
				else {
					circle.pop_back();
				}
			}
		}
	}
	//return circle;
	nodes = nodesemp;
	for (int i = 0; i < (int)circle.size(); i++) {
		for (int j = 0; j < (int)nodes.size(); j++) {
			if (circle[i] == nodes[j]) {
				nodes.erase(nodes.begin() + j);
				break;
			}
		}
	}

	while (!nodes.empty())
	{
		int theind = 0;
		int thenode = nodes[theind];
		int thepos = 0;
		double thedis = DBL_MAX;
		for (int i = 0; i < (int)nodes.size(); i++) {
			for (int j = 0; j < (int)circle.size(); j++) {
				int k = j;
				if (j == 0) k = circle.size();
				if (instance->distances[circle[k - 1]][nodes[i]] + instance->distances[nodes[i]][circle[k % circle.size()]]
					- instance->distances[circle[k - 1]][circle[k % circle.size()]] < thedis) {
					theind = i;
					thenode = nodes[theind];
					thepos = k;
					thedis = instance->distances[circle[k - 1]][nodes[i]] + instance->distances[nodes[i]][circle[k % circle.size()]] - instance->distances[circle[k - 1]][circle[k % circle.size()]];
				}
			}
		}
		circle.insert(circle.begin() + thepos, thenode);
		nodes.erase(nodes.begin() + theind);
	}

	bool flag = false;
	int depotpos = 0;
	for (int i = 0; i < (int)circle.size(); i++) {
		if (circle[i] == 0) {
			flag = true;
			depotpos = i;
			break;
		}
	}
	if (flag) {
		vector<int> newcircle;
		for (int i = 0; i < (int)circle.size(); i++) {
			newcircle.push_back(circle[depotpos]);
			depotpos++;
			if (depotpos == (int)circle.size()) {
				depotpos = 0;
			}
		}
		return newcircle;
	}
	else {
		return circle;
	}
}

//insert every time really needs to charge
pair<vector<int>, double> insertStationByNecessaryCharge(vector<int> route, Case* instance) {
	int counter = 1;
	double nowdis = 0;
	double totaldis = 0;
	while (counter < (int)route.size())
	{
		nowdis += instance->distances[route[counter - 1]][route[counter]];
		if (nowdis > instance->maxDis) {
			int thestation = -1;
			do
			{
				nowdis -= instance->distances[route[counter - 1]][route[counter]];
				thestation = instance->findNearestStationFeasible(route[counter - 1], route[counter], instance->maxDis - nowdis);
				counter--;
			} while (thestation == -1 && counter >= 1);
			if (thestation == -1) return make_pair(route, -1);
			route.insert(route.begin() + counter + 1, thestation);
			totaldis += (nowdis + instance->distances[route[counter]][thestation]);
			nowdis = 0;
			counter++;
		}
		counter++;
	}
	totaldis += nowdis;
	return make_pair(route, totaldis);
}

pair<vector<int>, double> insertStationByRemove(vector<int> route, Case* instance) {
	list<pair<int, int>> stationInserted;
	for (int i = 0; i < (int)route.size() - 1; i++) {
		double allowedDis = instance->maxDis;
		if (i != 0) {
			allowedDis = instance->maxDis - instance->distances[stationInserted.back().second][route[i]];
		}
		int onestation = instance->findNearestStationFeasible(route[i], route[i + 1], allowedDis);
		if (onestation == -1) return make_pair(route, -1);
		stationInserted.push_back(make_pair(i, onestation));
	}
	/*for (int i = 0; i < (int)route.size() - 1; i++) {
		stationInserted.push_back(make_pair(i, instance->bestStation[route[i]][route[i + 1]]));
	}*/
	while (!stationInserted.empty())
	{
		bool changed = false;
		list<pair<int, int>>::iterator delone = stationInserted.begin();
		double savedis = 0;
		list<pair<int, int>>::iterator itr = stationInserted.begin();
		list<pair<int, int>>::iterator next = itr;
		next++;
		if (next != stationInserted.end()) {
			int endInd = next->first;
			int endstation = next->second;
			double sumdis = 0;
			for (int i = 0; i < endInd; i++) {
				sumdis += instance->distances[route[i]][route[i + 1]];
			}
			sumdis += instance->distances[route[endInd]][endstation];
			if (sumdis <= instance->maxDis) {
				savedis = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
					- instance->distances[route[itr->first]][route[itr->first + 1]];
			}
		}
		else {
			double sumdis = 0;
			for (int i = 0; i < (int)route.size() - 1; i++) {
				sumdis += instance->distances[route[i]][route[i + 1]];
			}
			if (sumdis <= instance->maxDis) {
				savedis = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
					- instance->distances[route[itr->first]][route[itr->first + 1]];
			}
		}
		itr++;
		while (itr != stationInserted.end())
		{
			int startInd, endInd;
			next = itr;
			next++;
			list<pair<int, int>>::iterator prev = itr;
			prev--;
			double sumdis = 0;
			if (next != stationInserted.end()) {
				startInd = prev->first + 1;
				endInd = next->first;
				sumdis += instance->distances[prev->second][route[startInd]];
				for (int i = startInd; i < endInd; i++) {
					sumdis += instance->distances[route[i]][route[i + 1]];
				}
				sumdis += instance->distances[route[endInd]][next->second];
				if (sumdis <= instance->maxDis) {
					double savedistemp = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
						- instance->distances[route[itr->first]][route[itr->first + 1]];
					if (savedistemp > savedis) {
						savedis = savedistemp;
						delone = itr;
					}
				}
			}
			else {
				startInd = prev->first + 1;
				sumdis += instance->distances[prev->second][route[startInd]];
				for (int i = startInd; i < (int)route.size() - 1; i++) {
					sumdis += instance->distances[route[i]][route[i + 1]];
				}
				if (sumdis <= instance->maxDis) {
					double savedistemp = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
						- instance->distances[route[itr->first]][route[itr->first + 1]];
					if (savedistemp > savedis) {
						savedis = savedistemp;
						delone = itr;
					}
				}
			}
			itr++;
		}
		if (savedis != 0) {
			stationInserted.erase(delone);
			changed = true;
		}
		if (!changed) {
			break;
		}
	}
	while (!stationInserted.empty())
	{
		route.insert(route.begin() + stationInserted.back().first + 1, stationInserted.back().second);
		stationInserted.pop_back();
	}
	double summ = 0;
	for (int i = 0; i < (int)route.size() - 1; i++) {
		summ += instance->distances[route[i]][route[i + 1]];
	}
	return make_pair(route, summ);
}

pair<vector<int>, double> insertStationByRemove2(vector<int> route, Case* instance) {
	list<pair<int, int>> stationInserted;
	list<double> distanceRecord;
	for (int i = 0; i < (int)route.size() - 1; i++) {
		double allowedDis = instance->maxDis;
		if (i != 0) {
			allowedDis = instance->maxDis - instance->distances[stationInserted.back().second][route[i]];
		}
		int onestation = instance->findNearestStationFeasible(route[i], route[i + 1], allowedDis);
		if (onestation == -1) return make_pair(route, -1);
		stationInserted.push_back(make_pair(i, onestation));
		distanceRecord.push_back(0);
	}
	stationInserted.push_back(make_pair(route.size() - 1, 0));
	distanceRecord.push_back(0);
	while (stationInserted.size() > 1)
	{
		bool changed = false;
		double savedis = 0;
		list<pair<int, int>>::iterator sdelone = stationInserted.begin();
		list<pair<int, int>>::iterator sitr = stationInserted.begin();
		list<pair<int, int>>::iterator snext = sitr;
		list<double>::iterator ddelone = distanceRecord.begin();
		list<double>::iterator ditr = distanceRecord.begin();
		list<double>::iterator dnext = ditr;
		snext++;
		dnext++;
		double sumdis = *ditr + instance->distances[route[sitr->first]][route[sitr->first + 1]] + *dnext + instance->distances[route[snext->first]][snext->second];
		if (sumdis <= instance->maxDis) {
			savedis = instance->distances[route[sitr->first]][sitr->second] + instance->distances[sitr->second][route[sitr->first + 1]]
					- instance->distances[route[sitr->first]][route[sitr->first + 1]];
		}
		sitr++;
		snext++;
		ditr++;
		dnext++;
		while (snext != stationInserted.end())
		{
			list<pair<int, int>>::iterator prev = sitr;
			prev--;
			sumdis = instance->distances[prev->second][route[prev->first + 1]] + *ditr + instance->distances[route[sitr->first]][route[sitr->first + 1]] + *dnext + instance->distances[route[snext->first]][snext->second];
			if (sumdis <= instance->maxDis) {
				double savedistemp = instance->distances[route[sitr->first]][sitr->second] + instance->distances[sitr->second][route[sitr->first + 1]]
						- instance->distances[route[sitr->first]][route[sitr->first + 1]];
				if (savedistemp > savedis) {
					savedis = savedistemp;
					sdelone = sitr;
					ddelone = ditr;
				}
			}
			sitr++;
			snext++;
			ditr++;
			dnext++;
		}
		if (savedis != 0) {
			list<double>::iterator ddelonenext = ddelone;
			ddelonenext++;
			(*ddelonenext) =(*ddelonenext) + (*ddelone) + instance->distances[route[sdelone->first]][route[sdelone->first + 1]];
			stationInserted.erase(sdelone);
			distanceRecord.erase(ddelone);
			changed = true;
		}
		if (!changed) break;
	}
	stationInserted.pop_back();
	while (!stationInserted.empty())
	{
		route.insert(route.begin() + stationInserted.back().first + 1, stationInserted.back().second);
		stationInserted.pop_back();
	}
	double summ = 0;
	for (int i = 0; i < (int)route.size() - 1; i++) {
		summ += instance->distances[route[i]][route[i + 1]];
	}
	return make_pair(route, summ);
}

//starts and ends at 0
bool opt2noStation(vector<int>& route, Case* instance) {
	if (route.size() <= 4) return false;
	double minchange = 0;
	bool flag = false;
	do
	{
		minchange = 0;
		int mini = 0, minj = 0;
		for (int i = 0; i < (int)route.size() - 3; i++) {
			for (int j = i + 2; j < (int)route.size() - 1; j++) {
				double xx1 = instance->distances[route[i]][route[j]] + instance->distances[route[i + 1]][route[j + 1]];
				double xx2 = instance->distances[route[i]][route[i + 1]] + instance->distances[route[j]][route[j + 1]];
				double change = xx1 - xx2;
				if (fabs(change) < 0.00000001) change = 0;
				if (minchange > change) {
					minchange = change;
					mini = i;
					minj = j;
					flag = true;
				}
			}
		}
		if (minchange < 0)
			reverse(route.begin() + mini + 1, route.begin() + minj + 1);
	} while (minchange < 0);
	return flag;
}

bool opt2noStation2(int* route, int length, double& fitv, Case* instance) {
	if (length <= 4) return false;
	double minchange = 0;
	bool flag = false;
	do
	{
		minchange = 0;
		int mini = 0, minj = 0;
		for (int i = 0; i < length - 3; i++) {
			for (int j = i + 2; j < length - 1; j++) {
				double xx1 = instance->distances[route[i]][route[j]] + instance->distances[route[i + 1]][route[j + 1]];
				double xx2 = instance->distances[route[i]][route[i + 1]] + instance->distances[route[j]][route[j + 1]];
				double change = xx1 - xx2;
				if (fabs(change) < 0.00000001) change = 0;
				if (minchange > change) {
					minchange = change;
					mini = i;
					minj = j;
					flag = true;
				}
			}
		}
		if (minchange < 0) {
			reverseAnArray(route, mini + 1, minj);
			fitv += minchange;
		}
	} while (minchange < 0);
	return flag;
}

void reverseAnArray(int* route, int start, int end) {
	while (start < end)
	{
		int temp = route[start];
		route[start] = route[end];
		route[end] = temp;
		start++;
		end--;
	}
}

bool orNoStation(int* route, int length, double& fitv, Case* instance) {
	if (length <= 4) return false;
	double minchange = 0;
	bool flag = false;
	do
	{
		minchange = 0;
		int mini = 0, minj = 0;
		for (int i = 1; i < length - 1; i++) {
			for (int j = 1; j < length - 1; j++) {
				if (i < j) {
					double xx1 = instance->distances[route[i - 1]][route[i]] + instance->distances[route[i]][route[i + 1]] + instance->distances[route[j]][route[j + 1]];
					double xx2 = instance->distances[route[i - 1]][route[i + 1]] + instance->distances[route[j]][route[i]] + instance->distances[route[i]][route[j + 1]];
					double change = xx1 - xx2;
					if (fabs(change) < 0.00000001) change = 0;
					if (minchange < change) {
						minchange = change;
						mini = i;
						minj = j;
						flag = true;
					}
				}
				else if (i > j) {
					double xx1 = instance->distances[route[i - 1]][route[i]] + instance->distances[route[i]][route[i + 1]] + instance->distances[route[j - 1]][route[j]];
					double xx2 = instance->distances[route[j - 1]][route[i]] + instance->distances[route[i]][route[j]] + instance->distances[route[i - 1]][route[i + 1]];
					double change = xx1 - xx2;
					if (fabs(change) < 0.00000001) change = 0;
					if (minchange < change) {
						minchange = change;
						mini = i;
						minj = j;
						flag = true;
					}
				}
			}
		}
		if (minchange > 0) {
			moveItoJ(route, mini, minj);
			fitv -= minchange;
		}
	} while (minchange > 0);
	return flag;
}

void moveItoJ(int* route, int a, int b) {
	int x = route[a];
	if (a < b) {
		for (int i = a; i < b; i++) {
			route[i] = route[i + 1];
		}
		route[b] = x;
	}
	else if (a > b) {
		for (int i = a; i > b; i--) {
			route[i] = route[i - 1];
		}
		route[b] = x;
	}
}

pair<vector<int>, double> greedymethod(Case* instance) {
	set<int> remain;
	for (int i = instance->depotNumber; i < instance->depotNumber + instance->customerNumber; i++) {
		remain.insert(i);
	}
	vector<int> solution = {0};
	while (!remain.empty())
	{
		int thenearst = -1;
		int lastnode = solution.back();
		double mindis = DBL_MAX;
		for (auto e : remain) {
			if (instance->distances[e][lastnode] < mindis) {
				thenearst = e;
				mindis = instance->distances[e][lastnode];
			}
		}
		solution.push_back(thenearst);
		remain.erase(thenearst);
	}
	solution.push_back(0);

	double currentc = 0;
	int counter = 1;
	vector<vector<int>> routes;
	do
	{
		vector<int> aroute = {0};
		currentc = 0;
		while (currentc < instance->maxC)
		{
			if (counter < (int)solution.size() && currentc + instance->demand[solution[counter]] < instance->maxC) {
				aroute.push_back(solution[counter]);
				currentc += instance->demand[solution[counter]];
				counter++;
			}
			else {
				if (aroute.back() != 0) aroute.push_back(0);
				pair<vector<int>, double> xx = insertStationByRemove(aroute, instance);
				//cout << instance->calculateRouteDistance(xx.first) << endl;
				routes.push_back(xx.first);
				break;
			}
		}
	} while (counter < (int)solution.size());

	/*for (auto a : routes) {
		for (auto b : a) {
			cout << b << " ";
		}
		cout << endl;
	}*/

	vector<int> oneSolution;
	for (auto a : routes) {
		for (auto b : a) {
			oneSolution.push_back(b);
		}
		oneSolution.pop_back();
	}
	oneSolution.push_back(0);

	/*for (auto e : oneSolution) {
		cout << e << " ";
	}
	cout << endl;*/

	double totaldis = instance->calculateRouteDistance(oneSolution);
	return make_pair(oneSolution, totaldis);
}

bool compareArr(vector<int> x, int* y) {
	for (int i = 0; i < (int)x.size(); i++) {
		if (x[i] != y[i]) {
			return false;
		}
	}
	return true;
}

void prinsSplitAnt(Ant* x, Case* instance) {
	int arrlength = instance->depotNumber + instance->customerNumber;
	int* pp = new int[arrlength];
	double* vv = new double[arrlength];
	memset(pp, 0, sizeof(int) * arrlength);
	vv[0] = 0;
	for (int i = 1; i < arrlength; i++) {
		vv[i] = DBL_MAX;
	}
	for (int i = 1; i < arrlength; i++) {
		double load = 0;
		double cost = 0;
		int j = i;
		do
		{
			load += instance->demand[x->circle[j]];
			if (i == j) {
				cost = instance->distances[0][x->circle[j]] * 2;
			}
			else {
				cost = cost - instance->distances[x->circle[j - 1]][0];
			    cost = cost + instance->distances[x->circle[j - 1]][x->circle[j]];
			    cost = cost + instance->distances[0][x->circle[j]];
			}
			if (load <= (double)instance->maxC) {
				if (vv[i - 1] + cost < vv[j]) {
					vv[j] = vv[i - 1] + cost;
					pp[j] = i - 1;
				}
				j++;
			}
		} while (!(j >= arrlength || load > instance->maxC));
	}
	int j = arrlength - 1;
	int curR = 0;
	while (true)
	{
		int i = pp[j];
		x->route[curR][0] = 0;
		int routelength = 1;
		double loadf = 0;
		for (int k = i + 1; k <= j; k++) {
			x->route[curR][routelength] = x->circle[k];
			loadf += instance->demand[x->circle[k]];
			routelength++;
		}
		x->route[curR][routelength] = 0;
		routelength++;
		x->nodeNum[curR] = routelength;
		x->demsum[curR] = loadf;
		curR++;
		j = i;
		if (i == 0) {
			break;
		}
	}
	x->routeNum = curR;
	x->fit = vv[arrlength - 1];
	delete[] pp;
	delete[] vv;
}

//start from 0, but not end at 0
vector<vector<int>> prinsSplit(vector<int> x, Case* instance) {
	int arrlength = instance->depotNumber + instance->customerNumber;
	int* pp = new int[arrlength];
	double* vv = new double[arrlength];
	memset(pp, 0, sizeof(int) * arrlength);
	vv[0] = 0;
	for (int i = 1; i < arrlength; i++) {
		vv[i] = DBL_MAX;
	}
	for (int i = 1; i < (int)x.size(); i++) {
		double load = 0;
		double cost = 0;
		int j = i;
		do
		{
			load += instance->demand[x[j]];
			if (i == j) {
				cost = instance->distances[0][x[j]] * 2;
			}
			else {
				cost = cost - instance->distances[x[j - 1]][0];
			    cost = cost + instance->distances[x[j - 1]][x[j]];
			    cost = cost + instance->distances[0][x[j]];
			}
			if (load <= (double)instance->maxC) {
				if (vv[i - 1] + cost < vv[j]) {
					vv[j] = vv[i - 1] + cost;
					pp[j] = i - 1;
				}
				j++;
			}
		} while (!(j >= (int)x.size() || load > instance->maxC));
	}

	vector<vector<int>> allroutes;
	int j = x.size() - 1;
	while (true) {
		int i = pp[j];
		vector<int> temp(x.begin() + i + 1, x.begin() + j + 1);
		allroutes.push_back(temp);
		j = i;
		if (i == 0) {
			break;
		}
	}
	delete[] pp;
	delete[] vv;
	return allroutes;
}

pair<double, double> deficitAndDistance(vector<int>& route, Case* instance) {
	vector<double> deficit;
	//deficit.reserve(route.size());
	deficit.push_back(instance->maxDis);
	double totaldis = 0;
	for (int i = 1; i < (int)route.size(); i++) {
		totaldis += instance->distances[route[i - 1]][route[i]];
		if (route[i - 1] < instance->depotNumber + instance->customerNumber) {
			deficit.push_back(deficit.back() - instance->distances[route[i - 1]][route[i]]);
		}
		else {
			deficit.push_back(instance->maxDis - instance->distances[route[i - 1]][route[i]]);
		}
	}
	double totaldificit = 0;
	for (auto e : deficit) {
		if (e < 0) {
			totaldificit += e;
		}
	}
	return make_pair(totaldificit, totaldis);
}

pair<vector<int>, double> insertStationByMinDeficit(vector<int> route, Case* instance) {
	pair<double, double> xx = deficitAndDistance(route, instance);
	vector<int> shunxu;
	//shunxu.reserve(route.size());
	for (int i = 0; i < (int)route.size(); i++) {
		shunxu.push_back(0);
	}
	double dif = xx.first;
	double ttd = xx.second;
	int counter = 0;
	while (dif < 0)
	{
		counter++;
		int thesta = -1;
		int thepos = -1;
		for (int i = instance->depotNumber + instance->customerNumber; i < instance->depotNumber + instance->customerNumber + instance->stationNumber; i++) {
			for (int j = 1; j < (int)route.size(); j++) {
				vector<int> temp = route;
				temp.insert(temp.begin() + j, i);
				xx = deficitAndDistance(temp, instance);
				if (xx.first > dif) {
					dif = xx.first;
					ttd = xx.second;
					thesta = i;
					thepos = j;
				}
				if (xx.first == dif && xx.second < ttd) {
					ttd = xx.second;
					thesta = i;
					thepos = j;
				}
			}
		}
		if (thesta != -1) {
			route.insert(route.begin() + thepos, thesta);
			shunxu.insert(shunxu.begin() + thepos, counter);
		}
	}
	/*int counter2 = 1;
	while (counter >= 3)
	{
		int thepos = 0;
		for (int i = 0; i < (int)shunxu.size(); i++) {
			if (shunxu[i] == counter2) {
				thepos = i;
				break;
			}
		}
		int prepos = 0;
		int nexpos = shunxu.size() - 1;
		for (int i = thepos - 1; i >= 0; i--) {
			if (shunxu[i] != 0 || i == 0) {
				prepos = i;
				break;
			}
		}
		for (int i = thepos + 1; i < (int)shunxu.size(); i++) {
			if (shunxu[i] != 0 || i == (int)shunxu.size() - 1) {
				nexpos = i;
				break;
			}
		}
		double pieced = 0;
		for (int i = prepos; i < nexpos; i++) {
			if (i + 1 != thepos) {
				pieced += instance->distances[route[i]][route[i + 1]];
			}
			else {
				pieced += instance->distances[route[i]][route[i + 2]];
				i++;
			}
		}
		if (pieced <= instance->maxDis) {
			ttd = ttd - instance->distances[route[thepos]][route[thepos + 1]] - instance->distances[route[thepos]][route[thepos - 1]] + instance->distances[route[thepos - 1]][route[thepos + 1]];
			route.erase(route.begin() + thepos);
			shunxu.erase(shunxu.begin() + thepos);
			counter--;
			counter2++;
		}
		else {
			break;
		}
	}*/
	return make_pair(route, ttd);
}

//to be implemented, with the consideration of station during splitting
vector<vector<int>> prinsSplit2(vector<int> x, Case* instance) {
	int arrlength = instance->depotNumber + instance->customerNumber;
	int* pp = new int[arrlength];
	double* vv = new double[arrlength];
	memset(pp, 0, sizeof(int) * arrlength);
	vv[0] = 0;
	for (int i = 1; i < arrlength; i++) {
		vv[i] = DBL_MAX;
	}
	for (int i = 1; i < (int)x.size(); i++) {
		double load = 0;
		double cost = 0;
		int j = i;
		do
		{ 
			load += instance->demand[x[j]];
			if (i == j) {
				cost = instance->distances[0][x[j]] * 2;
				if (cost > instance->maxDis) {
					vector<int> temparr;
					temparr.push_back(0);
					for (int k = i; k <= j; j++) {
						temparr.push_back(x[k]);
					}
					temparr.push_back(0);
					pair<vector<int>, double> xx = insertStationByRemove(temparr, instance);
					cost = xx.second;
				}
			}
			else {
				cost = cost - instance->distances[x[j - 1]][0];
				cost = cost + instance->distances[x[j - 1]][x[j]];
				cost = cost + instance->distances[0][x[j]];
				if (cost > instance->maxDis) {
					vector<int> temparr;
					temparr.push_back(0);
					for (int k = i; k <= j; j++) {
						temparr.push_back(x[k]);
					}
					temparr.push_back(0);
					pair<vector<int>, double> xx = insertStationByRemove(temparr, instance);
					cost = xx.second;
				}
			}
			if (load <= (double)instance->maxC) {
				if (vv[i - 1] + cost < vv[j]) {
					vv[j] = vv[i - 1] + cost;
					pp[j] = i - 1;
				}
				j++;
			}
		} while (!(j >= (int)x.size() || load > instance->maxC));
	}
	vector<vector<int>> allroutes;
	for (int i = 1; i < (int)x.size(); i++) {
		vector<int> onetemp;
		allroutes.push_back(onetemp);
	}
	int t = 0;
	int j = x.size() - 1;
	while (true) {
		t++;
		int i = pp[j];
		for (int k = i + 1; k <= j; k++) {
			allroutes[t].push_back(x[k]);
		}
		j = i;
		if (i == 0) {
			break;
		}
	}
	vector<vector<int>> newroutes;
	for (auto e : allroutes) {
		if (!e.empty()) {
			newroutes.push_back(e);
		}
	}
	delete[] pp;
	delete[] vv;
	return newroutes;
}

//starts and ends at 0
bool exchangeNoStation(vector<int>& route, Case* instance) {
	if (route.size() <= 4) return false;
	double minchange = 0;
	bool flag = false;
	do
	{
		minchange = 0;
		int mini = 0, minj = 0;
		for (int i = 1; i < (int)route.size() - 2; i++) {
			for (int j = i + 1; j < (int)route.size() - 1; j++) {
				double orid = instance->distances[route[i - 1]][route[i]] + instance->distances[route[i + 1]][route[i]] + instance->distances[route[j - 1]][route[j]] + instance->distances[route[j + 1]][route[j]];
				double aftd = instance->distances[route[i - 1]][route[j]] + instance->distances[route[i + 1]][route[j]] + instance->distances[route[j - 1]][route[i]] + instance->distances[route[j + 1]][route[i]];
				double change = aftd - orid;
				if (fabs(change) < 0.00000001) change = 0;
				if (minchange > change) {
					minchange = change;
					mini = i;
					minj = j;
					flag = true;
				}
			}
		}
		int temp = route[mini];
		route[mini] = route[minj];
		route[minj] = temp;
	} while (minchange < 0);
	return flag;
}

pair<vector<int>, double> insertStationByRemove3(vector<int> route, Case* instance) {
	vector<pair<int, int>> stationsInserted;
	vector<double> distInfor(route.size(), 0);
	for (int i = 0; i < (int)route.size() - 1; i++) {
		int thestation = instance->findNearestStationFeasible2(route[i], route[i + 1], instance->maxDis);
		stationsInserted.push_back(make_pair(i, thestation));
		distInfor[i + 1] = distInfor[i] + instance->distances[route[i + 1]][route[i]];
	}
	//bool suc = false;
	//insert one station
	int theone = -1;
	double theonesavedis = DBL_MAX;
	for (int i = 0; i < (int)stationsInserted.size(); i++) {
		double formerd = distInfor[stationsInserted[i].first] + instance->distances[route[stationsInserted[i].first]][stationsInserted[i].second];
		if (formerd <= instance->maxDis) {
			double latterd = distInfor.back() - distInfor[stationsInserted[i].first + 1] + instance->distances[route[stationsInserted[i].first + 1]][stationsInserted[i].second];
			if (latterd <= instance->maxDis) {
				double extradis = instance->distances[route[stationsInserted[i].first]][stationsInserted[i].second] + instance->distances[route[stationsInserted[i].first + 1]][stationsInserted[i].second];
				extradis -= instance->distances[route[stationsInserted[i].first]][route[stationsInserted[i].first + 1]];
				if (extradis < theonesavedis) {
					theone = i;
					theonesavedis = extradis;
				}
			}
		}
		else {
			break;
		}
	}
	//insert two stations
	int thetwo1 = -1;
	int thetwo2 = -1;
	double thetwosavedis = DBL_MAX;
	for (int i = 0; i < (int)stationsInserted.size() - 1; i++) {
		double formerd = distInfor[stationsInserted[i].first] + instance->distances[route[stationsInserted[i].first]][stationsInserted[i].second];
		if (formerd <= instance->maxDis) {
			for (int j = i + 1; j < (int)stationsInserted.size(); j++) {
				double midis = distInfor[stationsInserted[j].first] - distInfor[stationsInserted[i].first + 1];
				midis += instance->distances[stationsInserted[i].second][route[stationsInserted[i].first + 1]];
				midis += instance->distances[route[stationsInserted[j].first]][stationsInserted[j].second];
				if (midis <= instance->maxDis) {
					double latterd = distInfor.back() - distInfor[stationsInserted[j].first + 1] + instance->distances[route[stationsInserted[j].first + 1]][stationsInserted[j].second];
					if (latterd <= instance->maxDis) {
						double extradis1 = instance->distances[route[stationsInserted[i].first]][stationsInserted[i].second] + instance->distances[route[stationsInserted[i].first + 1]][stationsInserted[i].second];
						extradis1 -= instance->distances[route[stationsInserted[i].first]][route[stationsInserted[i].first + 1]];
						double extradis2 = instance->distances[route[stationsInserted[j].first]][stationsInserted[j].second] + instance->distances[route[stationsInserted[j].first + 1]][stationsInserted[j].second];
						extradis2 -= instance->distances[route[stationsInserted[j].first]][route[stationsInserted[j].first + 1]];
						if (extradis1 + extradis2 < thetwosavedis) {
							thetwo1 = i;
							thetwo2 = j;
							thetwosavedis = extradis1 + extradis2;
						}
					}
				}
				else {
					break;
				}
			}
		}
		else {
			break;
		}
	}
	if (theone == -1 && thetwo1 == -1) {
		pair<vector<int>, double> xx = insertStationByRemove(route, instance);
		return xx;
	}
	else {
		double allcost = distInfor.back();
		if (theonesavedis < thetwosavedis) {
			route.insert(route.begin() + stationsInserted[theone].first + 1, stationsInserted[theone].second);
			allcost += theonesavedis;
		}
		else {
			route.insert(route.begin() + stationsInserted[thetwo2].first + 1, stationsInserted[thetwo2].second);
			route.insert(route.begin() + stationsInserted[thetwo1].first + 1, stationsInserted[thetwo1].second);
			allcost += thetwosavedis;
		}
		return make_pair(route, allcost);
	}
}

//pair<vector<int>, double> insertStationByEnumerate(vector<int> route, Case* instance) {
//	//get the accumulated distance vector
//	vector<double> accumulateDistance(route.size());
//	accumulateDistance[0] = 0;
//	for (int i = 1; i < (int)route.size(); i++) {
//		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
//	}
//	//get the best station list vector
//	//did not consider the extreme situation that there is no feasible station between two customers.
//	vector<int> bestStationList(route.size() - 1);
//	for (int i = 0; i < (int)route.size() - 1; i++) {
//		int astation = instance->findNearestStationFeasible(route[i], route[i + 1], instance->maxDis);
//		bestStationList[i] = astation;
//	}
//	//get the candidate list vector<vector>
//	//did not consider the extreme situation that there is no feasible station between two customers
//	vector<vector<int>> candidateStationList;
//	for (int i = 0; i < (int)route.size() - 1; i++) {
//		vector<int> acandilist;
//		for (int j = instance->customerNumber + instance->depotNumber; j < instance->customerNumber + instance->depotNumber + instance->stationNumber; j++) {
//			bool bedominated = false;
//			for (auto e : acandilist) {
//				if (instance->distances[route[i]][e] <= instance->distances[route[i]][j] &&
//					instance->distances[e][route[i + 1]] <= instance->distances[j][route[i + 1]]) {
//					bedominated = true;
//					break;
//				}
//			}
//			if (bedominated == false) {
//				for (int k = 0; k < (int)acandilist.size(); k++) {
//					if (instance->distances[route[i]][j] <= instance->distances[route[i]][acandilist[k]] &&
//						instance->distances[j][route[i + 1]] <= instance->distances[acandilist[k]][route[j + 1]]) {
//						acandilist.erase(acandilist.begin() + k);
//						k--;
//					}
//				}
//				acandilist.push_back(j);
//			}
//		}
//		candidateStationList.push_back(acandilist);
//	}
//	//get the upper bound of number of stations by remove heuristic
//	pair<vector<int>, double> xx = insertStationByRemove(route, instance);
//	int ub = 0;
//	for (auto e : xx.first) {
//		if (e >= instance->depotNumber + instance->customerNumber) {
//			ub++;
//		}
//	}
//
//}

pair<vector<int>, double> insertStationByEnumerationWithStation(vector<int> route, Case* instance) {
	vector<int> newroute;
	int lb = 0;
	for (int i = 0; i < (int)route.size(); i++) {
		if (route[i] >= instance->depotNumber + instance->customerNumber) {
			lb++;
		}
		else {
			newroute.push_back(route[i]);
		}
	}
	int ub = lb + 1;

	vector<double> accumulateDistance(newroute.size(), 0);
	for (int i = 1; i < (int)newroute.size(); i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[newroute[i]][newroute[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return make_pair(newroute, accumulateDistance.back());
	}

	int* chosenPos = new int[newroute.size()];
	int* chosenSta = new int[newroute.size()];
	vector<int> finalRoute;
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainN(0, i, chosenSta, chosenPos, finalRoute, finalfit, i, newroute, accumulateDistance, instance);
	}
	delete[] chosenPos;
	delete[] chosenSta;
	if (finalfit != DBL_MAX) {
		return make_pair(finalRoute, finalfit);
	}
	else {
		return make_pair(newroute, -1);
	}
}

pair<vector<int>, double> insertStationByEnumeration(vector<int> route, Case* instance) {
	vector<double> accumulateDistance(route.size(), 0);
	for (int i = 1; i < (int)route.size(); i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return make_pair(route, accumulateDistance.back());
	}

	int ub = (int)(accumulateDistance.back() / instance->maxDis + 3);
	int lb = (int)(accumulateDistance.back() / instance->maxDis);
	int* chosenPos = new int[route.size()];
	int* chosenSta = new int[route.size()];
	vector<int> finalRoute;
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainN(0, i, chosenSta, chosenPos, finalRoute, finalfit, i, route, accumulateDistance, instance);
	}
	delete[] chosenPos;
	delete[] chosenSta;
	if (finalfit != DBL_MAX) {
		return make_pair(finalRoute, finalfit);
	}
	else {
		return make_pair(route, -1);
	}
}

pair<vector<int>, double> insertStationByEnumeration(vector<int> route, Case* instance, int ub) {
	vector<double> accumulateDistance(route.size(), 0);
	for (int i = 1; i < (int)route.size(); i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return make_pair(route, accumulateDistance.back());
	}

	//int ub = (int)(accumulateDistance.back() / instance->maxDis + 3);
	int lb = (int)(accumulateDistance.back() / instance->maxDis);
	int* chosenPos = new int[route.size()];
	int* chosenSta = new int[route.size()];
	vector<int> finalRoute;
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainN(0, i, chosenSta, chosenPos, finalRoute, finalfit, i, route, accumulateDistance, instance);
	}
	delete[] chosenPos;
	delete[] chosenSta;
	if (finalfit != DBL_MAX) {
		return make_pair(finalRoute, finalfit);
	}
	else {
		return make_pair(route, -1);
	}
}

void tryACertainN(int mlen, int nlen, int* chosenSta, int* chosenPos, vector<int>& finalRoute, double& finalfit, int curub, vector<int>& route, vector<double>& accumulateDis, Case* instance) {
	for (int i = mlen; i <= (int)route.size() - 1 - nlen; i++) {
		if (curub == nlen) {
			if (accumulateDis[i] >= instance->maxDis) {
				break;
			}
		}
		else {
			if (accumulateDis[i] - accumulateDis[chosenPos[curub - nlen - 1] + 1] >= instance->maxDis) {
				break;
			}
		}
		if (nlen == 1) {
			if (accumulateDis.back() - accumulateDis[i + 1] >= instance->maxDis) {
				continue;
			}
		}
		for (int j = 0; j < (int)instance->bestStations[route[i]][route[i + 1]].size(); j++) {
			chosenSta[curub - nlen] = instance->bestStations[route[i]][route[i + 1]][j];
			chosenPos[curub - nlen] = i;
			if (nlen > 1) {
				tryACertainN(i + 1, nlen - 1, chosenSta, chosenPos, finalRoute, finalfit, curub, route, accumulateDis, instance);
			}
			else {
				bool feasible = true;
				double piecedis = accumulateDis[chosenPos[0]] + instance->distances[route[chosenPos[0]]][chosenSta[0]];
				if (piecedis > instance->maxDis) feasible = false;
				for (int k = 1; feasible && k < curub; k++) {
					piecedis = accumulateDis[chosenPos[k]] - accumulateDis[chosenPos[k - 1] + 1];
					piecedis += instance->distances[chosenSta[k - 1]][route[chosenPos[k - 1] + 1]];
					piecedis += instance->distances[chosenSta[k]][route[chosenPos[k]]];
					if (piecedis > instance->maxDis) feasible = false;
				}
				piecedis = accumulateDis.back() - accumulateDis[chosenPos[curub - 1] + 1];
				piecedis += instance->distances[route[chosenPos[curub - 1] + 1]][chosenSta[curub - 1]];
				if (piecedis > instance->maxDis) feasible = false;
				if (feasible) {
					double totaldis = accumulateDis.back();
					for (int k = 0; k < curub; k++) {
						int firstnode = route[chosenPos[k]];
						int secondnode = route[chosenPos[k] + 1];
						totaldis -= instance->distances[firstnode][secondnode];
						totaldis += instance->distances[firstnode][chosenSta[k]];
						totaldis += instance->distances[chosenSta[k]][secondnode];
					}
					if (totaldis < finalfit) {
						finalfit = totaldis;
						finalRoute = route;
						for (int k = curub - 1; k >= 0; k--) {
							finalRoute.insert(finalRoute.begin() + chosenPos[k] + 1, chosenSta[k]);
						}
					}
				}
			}
		}
	}
}

pair<vector<int>, double> insertStationBySimpleEnumeration(vector<int> route, Case* instance) {
	vector<double> accumulateDistance(route.size(), 0);
	for (int i = 1; i < (int)route.size(); i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return make_pair(route, accumulateDistance.back());
	}

	int ub = (int)(accumulateDistance.back() / instance->maxDis + 1);
	int lb = (int)(accumulateDistance.back() / instance->maxDis);
	int* chosenPos = new int[route.size()];
	vector<int> finalRoute;
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainN(0, i, chosenPos, finalRoute, finalfit, i, route, accumulateDistance, instance);
	}
	delete[] chosenPos;
	if (finalfit != DBL_MAX) {
		return make_pair(finalRoute, finalfit);
	}
	else {
		return make_pair(route, -1);
	}
}

double insertStationBySimpleEnumeration2(vector<int> route, Case* instance) {
	vector<double> accumulateDistance(route.size(), 0);
	for (int i = 1; i < (int)route.size(); i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return accumulateDistance.back();
	}

	int ub = (int)(accumulateDistance.back() / instance->maxDis + 1);
	int lb = (int)(accumulateDistance.back() / instance->maxDis);
	int* chosenPos = new int[route.size()];
	vector<int> finalRoute;
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainN(0, i, chosenPos, finalRoute, finalfit, i, route, accumulateDistance, instance);
	}
	delete[] chosenPos;
	if (finalfit != DBL_MAX) {
		return finalfit;
	}
	else {
		return -1;
	}
}

void tryACertainN(int mlen, int nlen, int* chosenPos, vector<int>& finalRoute, double& finalfit, int curub, vector<int>& route, vector<double>& accumulateDis, Case* instance) {
	for (int i = mlen; i <= (int)route.size() - 1 - nlen; i++) {
		if (curub == nlen) {
			double onedis = instance->distances[route[i]][instance->bestStation[route[i]][route[i + 1]]];
			if (accumulateDis[i] + onedis > instance->maxDis) {
				break;
			}
		}
		else {
			int lastpos = chosenPos[curub - nlen - 1];
			double onedis = instance->distances[route[lastpos + 1]][instance->bestStation[route[lastpos]][route[lastpos + 1]]];
			double twodis = instance->distances[route[i]][instance->bestStation[route[i]][route[i + 1]]];
			if (accumulateDis[i] - accumulateDis[lastpos + 1] + onedis + twodis > instance->maxDis) {
				break;
			}
		}
		if (nlen == 1) {
			double onedis = accumulateDis.back() - accumulateDis[i + 1] + instance->distances[instance->bestStation[route[i]][route[i + 1]]][route[i + 1]];
			if (onedis > instance->maxDis) {
				continue;
			}
		}

		chosenPos[curub - nlen] = i;
		if (nlen > 1) {
			tryACertainN(i + 1, nlen - 1, chosenPos, finalRoute, finalfit, curub, route, accumulateDis, instance);
		}
		else {
			double disum = accumulateDis.back();
			for (int j = 0; j < curub; j++) {
				int firstnode = route[chosenPos[j]];
				int secondnode = route[chosenPos[j] + 1];
				int thestation = instance->bestStation[firstnode][secondnode];
				disum -= instance->distances[firstnode][secondnode];
				disum += instance->distances[firstnode][thestation];
				disum += instance->distances[secondnode][thestation];
			}
			if (disum < finalfit) {
				// finalRoute = route;
				// for (int j = curub - 1; j >= 0; j--) {
				// 	finalRoute.insert(finalRoute.begin() + chosenPos[j] + 1, instance->bestStation[finalRoute[chosenPos[j]]][finalRoute[chosenPos[j] + 1]]);
				// }
				finalfit = disum;
			}
		}
	}
}

//void evaluateAnInsertion(int* chosenPos, int curub, vector<int>& finalRoute, double& finalfit, vector<int>& route, Case* instance) {
//	vector<int> fixRoute = route;
//	for (int i = curub - 1; i >= 0; i--) {
//		fixRoute.insert(fixRoute.begin() + chosenPos[i] + 1, instance->bestStations[fixRoute[chosenPos[i]]][fixRoute[chosenPos[i] + 1]]);
//	}
//	double sumdis = 0;
//	double piecedis = 0;
//	for (int i = 1; i < (int)fixRoute.size(); i++) {
//		piecedis += instance->distances[fixRoute[i]][fixRoute[i - 1]];
//		sumdis += instance->distances[fixRoute[i]][fixRoute[i - 1]];
//		if (piecedis > instance->maxDis) {
//			return;
//		}
//		if (fixRoute[i] == 0 || fixRoute[i] >= instance->depotNumber + instance->customerNumber) {
//			piecedis = 0;
//		}
//	}
//	if (sumdis < finalfit) {
//		finalRoute = fixRoute;
//		finalfit = sumdis;
//	}
//}

bool opt2StarNoStation(Ant* ant, Case* instance) {
	if (ant->routeNum == 1) {
		return false;
	}
	int* tempr = new int[ant->nodeCap];
	int* tempr2 = new int[ant->nodeCap];
	bool updated = false;
	bool updated2 = false;
	do
	{
		updated2 = false;
		//iterate route first
		for (int r1 = 0; r1 < ant->routeNum - 1; r1++) {
			//iterate route second
			for (int r2 = r1 + 1; r2 < ant->routeNum; r2++) {
				//iterate all links in the first rotue
				double frdem = 0;
				for (int n1 = 0; n1 < ant->nodeNum[r1] - 1; n1++) {
					frdem += instance->demand[ant->route[r1][n1]];
					//iterate all links in the second route
					double srdem = 0;
					for (int n2 = 0; n2 < ant->nodeNum[r2] - 1; n2++) {
						srdem += instance->demand[ant->route[r2][n2]];
						//judge two situations
						//swith tails
						if (frdem + ant->demsum[r2] - srdem <= instance->maxC && srdem + ant->demsum[r1] - frdem <= instance->maxC) {
							double xx1 = instance->distances[ant->route[r1][n1]][ant->route[r1][n1 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r2][n2 + 1]];
							double xx2 = instance->distances[ant->route[r1][n1]][ant->route[r2][n2 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r1][n1 + 1]];
							double change = xx1 - xx2;
							//if can be improved
							if (change > 0.00000001) {
								memcpy(tempr, ant->route[r1], sizeof(int) * ant->nodeCap);
								int counter1 = n1 + 1;
								for (int i = n2 + 1; i < ant->nodeNum[r2]; i++) {
									ant->route[r1][counter1] = ant->route[r2][i];
									counter1++;
								}
								int counter2 = n2 + 1;
								for (int i = n1 + 1; i < ant->nodeNum[r1]; i++) {
									ant->route[r2][counter2] = tempr[i];
									counter2++;
								}
								ant->nodeNum[r1] = counter1;
								ant->nodeNum[r2] = counter2;
								double newdemsum1 = frdem + ant->demsum[r2] - srdem;
								double newdemsum2 = srdem + ant->demsum[r1] - frdem;
								ant->demsum[r1] = newdemsum1;
								ant->demsum[r2] = newdemsum2;
								updated = true;
								updated2 = true;
								//judge whether a route already empty
								if (ant->demsum[r1] == 0) {
									int* tempp = ant->route[r1];
									ant->route[r1] = ant->route[ant->routeNum - 1];
									ant->route[ant->routeNum - 1] = tempp;
									ant->demsum[r1] = ant->demsum[ant->routeNum - 1];
									ant->nodeNum[r1] = ant->nodeNum[ant->routeNum - 1];
									ant->routeNum--;
								}
								if (ant->demsum[r2] == 0) {
									int* tempp = ant->route[r2];
									ant->route[r2] = ant->route[ant->routeNum - 1];
									ant->route[ant->routeNum - 1] = tempp;
									ant->demsum[r2] = ant->demsum[ant->routeNum - 1];
									ant->nodeNum[r2] = ant->nodeNum[ant->routeNum - 1];
									ant->routeNum--;
								}
								break;
							}
						}
						else if (frdem + srdem <= instance->maxC && ant->demsum[r1] - frdem + ant->demsum[r2] - srdem <= instance->maxC) {
							double xx1 = instance->distances[ant->route[r1][n1]][ant->route[r1][n1 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r2][n2 + 1]];
							double xx2 = instance->distances[ant->route[r1][n1]][ant->route[r2][n2]] + instance->distances[ant->route[r1][n1 + 1]][ant->route[r2][n2 + 1]];
							double change = xx1 - xx2;
							if (change > 0.00000001) {
								memcpy(tempr, ant->route[r1], sizeof(int) * ant->nodeCap);
								int counter1 = n1 + 1;
								for (int i = n2; i >= 0; i--) {
									ant->route[r1][counter1] = ant->route[r2][i];
									counter1++;
								}
								int counter2 = 0;
								for (int i = ant->nodeNum[r1] - 1; i >= n1 + 1; i--) {
									tempr2[counter2] = tempr[i];
									counter2++;
								}
								for (int i = n2 + 1; i < ant->nodeNum[r2]; i++) {
									tempr2[counter2] = ant->route[r2][i];
									counter2++;
								}
								memcpy(ant->route[r2], tempr2, sizeof(int) * ant->nodeCap);
								ant->nodeNum[r1] = counter1;
								ant->nodeNum[r2] = counter2;
								double newdemsum1 = frdem + srdem;
								double newdemsum2 = ant->demsum[r1] + ant->demsum[r2] - frdem - srdem;
								ant->demsum[r1] = newdemsum1;
								ant->demsum[r2] = newdemsum2;
								updated = true;
								updated2 = true;
								if (ant->demsum[r1] == 0) {
									int* tempp = ant->route[r1];
									ant->route[r1] = ant->route[ant->routeNum - 1];
									ant->route[ant->routeNum - 1] = tempp;
									ant->demsum[r1] = ant->demsum[ant->routeNum - 1];
									ant->nodeNum[r1] = ant->nodeNum[ant->routeNum - 1];
									ant->routeNum--;
								}
								if (ant->demsum[r2] == 0) {
									int* tempp = ant->route[r2];
									ant->route[r2] = ant->route[ant->routeNum - 1];
									ant->route[ant->routeNum - 1] = tempp;
									ant->demsum[r2] = ant->demsum[ant->routeNum - 1];
									ant->nodeNum[r2] = ant->nodeNum[ant->routeNum - 1];
									ant->routeNum--;
								}
								break;
							}
						}
					}
					if (updated2) break;
				}
				if (updated2) break;
			}
			if (updated2) break;
		}
	} while (updated2);
	delete[] tempr;
	delete[] tempr2;
	return updated;
}

bool opt2starNoStation2(Ant* ant, Case* instance) {
	if (ant->routeNum == 1) {
		return false;
	}
	unordered_set<pair<int, int>, pair_hash> routepairs;
	for (int i = 0; i < ant->routeNum - 1; i++) {
		for (int j = i + 1; j < ant->routeNum; j++) {
			routepairs.insert(make_pair(i, j));
		}
	}
	int* tempr = new int[ant->nodeCap];
	int* tempr2 = new int[ant->nodeCap];
	bool updated = false;
	bool updated2 = false;
	while (!routepairs.empty())
	{
		updated2 = false;
		int r1 = routepairs.begin()->first;
		int r2 = routepairs.begin()->second;
		routepairs.erase(routepairs.begin());
		double frdem = 0;
		for (int n1 = 0; n1 < ant->nodeNum[r1] - 1; n1++) {
			frdem += instance->demand[ant->route[r1][n1]];
			double srdem = 0;
			for (int n2 = 0; n2 < ant->nodeNum[r2] - 1; n2++) {
				srdem += instance->demand[ant->route[r2][n2]];
				if (frdem + ant->demsum[r2] - srdem <= instance->maxC && srdem + ant->demsum[r1] - frdem <= instance->maxC) {
					double xx1 = instance->distances[ant->route[r1][n1]][ant->route[r1][n1 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r2][n2 + 1]];
					double xx2 = instance->distances[ant->route[r1][n1]][ant->route[r2][n2 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r1][n1 + 1]];
					double change = xx1 - xx2;
					if (change > 0.00000001) {
						ant->fit -= change;
						memcpy(tempr, ant->route[r1], sizeof(int) * ant->nodeCap);
						int counter1 = n1 + 1;
						for (int i = n2 + 1; i < ant->nodeNum[r2]; i++) {
							ant->route[r1][counter1] = ant->route[r2][i];
							counter1++;
						}
						int counter2 = n2 + 1;
						for (int i = n1 + 1; i < ant->nodeNum[r1]; i++) {
							ant->route[r2][counter2] = tempr[i];
							counter2++;
						}
						ant->nodeNum[r1] = counter1;
						ant->nodeNum[r2] = counter2;
						double newdemsum1 = frdem + ant->demsum[r2] - srdem;
						double newdemsum2 = srdem + ant->demsum[r1] - frdem;
						ant->demsum[r1] = newdemsum1;
						ant->demsum[r2] = newdemsum2;
						updated = true;
						updated2 = true;
						for (int i = 0; i < r1; i++) {
							routepairs.insert({i, r1});
						}
						for (int i = 0; i < r2; i++) {
							routepairs.insert({i, r2});
						}
						if (ant->demsum[r1] == 0) {
							int* tempp = ant->route[r1];
							ant->route[r1] = ant->route[ant->routeNum - 1];
							ant->route[ant->routeNum - 1] = tempp;
							ant->demsum[r1] = ant->demsum[ant->routeNum - 1];
							ant->nodeNum[r1] = ant->nodeNum[ant->routeNum - 1];
							ant->routeNum--;
							for (int i = 0; i < ant->routeNum; i++) {
								routepairs.erase({i, ant->routeNum});
							}
						}
						if (ant->demsum[r2] == 0) {
							int* tempp = ant->route[r2];
							ant->route[r2] = ant->route[ant->routeNum - 1];
							ant->route[ant->routeNum - 1] = tempp;
							ant->demsum[r2] = ant->demsum[ant->routeNum - 1];
							ant->nodeNum[r2] = ant->nodeNum[ant->routeNum - 1];
							ant->routeNum--;
							for (int i = 0; i < ant->routeNum; i++) {
								routepairs.erase({i, ant->routeNum});
							}
						}
						break;
					}
				}
				else if (frdem + srdem <= instance->maxC && ant->demsum[r1] - frdem + ant->demsum[r2] - srdem <= instance->maxC) {
					double xx1 = instance->distances[ant->route[r1][n1]][ant->route[r1][n1 + 1]] + instance->distances[ant->route[r2][n2]][ant->route[r2][n2 + 1]];
					double xx2 = instance->distances[ant->route[r1][n1]][ant->route[r2][n2]] + instance->distances[ant->route[r1][n1 + 1]][ant->route[r2][n2 + 1]];
					double change = xx1 - xx2;
					if (change > 0.00000001) {
						ant->fit -= change;
						memcpy(tempr, ant->route[r1], sizeof(int) * ant->nodeCap);
						int counter1 = n1 + 1;
						for (int i = n2; i >= 0; i--) {
							ant->route[r1][counter1] = ant->route[r2][i];
							counter1++;
						}
						int counter2 = 0;
						for (int i = ant->nodeNum[r1] - 1; i >= n1 + 1; i--) {
							tempr2[counter2] = tempr[i];
							counter2++;
						}
						for (int i = n2 + 1; i < ant->nodeNum[r2]; i++) {
							tempr2[counter2] = ant->route[r2][i];
							counter2++;
						}
						memcpy(ant->route[r2], tempr2, sizeof(int) * ant->nodeCap);
						ant->nodeNum[r1] = counter1;
						ant->nodeNum[r2] = counter2;
						double newdemsum1 = frdem + srdem;
						double newdemsum2 = ant->demsum[r1] + ant->demsum[r2] - frdem - srdem;
						ant->demsum[r1] = newdemsum1;
						ant->demsum[r2] = newdemsum2;
						updated = true;
						updated2 = true;
						for (int i = 0; i < r1; i++) {
							routepairs.insert({i, r1});
						}
						for (int i = 0; i < r2; i++) {
							routepairs.insert({i, r2});
						}
						if (ant->demsum[r1] == 0) {
							int* tempp = ant->route[r1];
							ant->route[r1] = ant->route[ant->routeNum - 1];
							ant->route[ant->routeNum - 1] = tempp;
							ant->demsum[r1] = ant->demsum[ant->routeNum - 1];
							ant->nodeNum[r1] = ant->nodeNum[ant->routeNum - 1];
							ant->routeNum--;
							for (int i = 0; i < ant->routeNum; i++) {
								routepairs.erase({i, ant->routeNum});
							}
						}
						if (ant->demsum[r2] == 0) {
							int* tempp = ant->route[r2];
							ant->route[r2] = ant->route[ant->routeNum - 1];
							ant->route[ant->routeNum - 1] = tempp;
							ant->demsum[r2] = ant->demsum[ant->routeNum - 1];
							ant->nodeNum[r2] = ant->nodeNum[ant->routeNum - 1];
							ant->routeNum--;
							for (int i = 0; i < ant->routeNum; i++) {
								routepairs.erase({i, ant->routeNum});
							}
						}
						break;
					}
				}
			}
			if (updated2) break;
		}
	}
	delete[] tempr;
	delete[] tempr2;
	return updated;
}

void opt2ToAnt(Ant* ant, Case* instance) {
	for (int i = 0; i < ant->routeNum; i++) {
		opt2noStation2(ant->route[i], ant->nodeNum[i], ant->fit, instance);
	}
}

void orToAnt(Ant* ant, Case* instance) {
	for (int i = 0; i < ant->routeNum; i++) {
		orNoStation(ant->route[i], ant->nodeNum[i], ant->fit, instance);
	}
}

//have not got a clear idea
// bool opt2starNeighborNoStation(Solution route, Case* instance) {
// 	if (route.leftNodes.size() == 1) {
// 		return false;
// 	}
// 	bool updated1 = false;
// 	bool updated2 = false;
// 	int startRoute = 0;
// 	do
// 	{
// 		updated2 = false;
// 		for (int r1 = startRoute; r1 < (int)route.leftNodes.size() - 1; r1++) {
// 			int s = 0;
// 			int e = route.rightNodes[r1];
// 			do
// 			{
				
// 			} while (!updated2 && e != 0);
// 		}
// 	} while (updated2);
	
// }

double insertStationBySimpleEnumerationArray(int* route, int length, Case* instance) {
	vector<double> accumulateDistance(length, 0);
	for (int i = 1; i < length; i++) {
		accumulateDistance[i] = accumulateDistance[i - 1] + instance->distances[route[i]][route[i - 1]];
	}
	if (accumulateDistance.back() <= instance->maxDis) {
		return accumulateDistance.back();
	}

	int ub = (int)(accumulateDistance.back() / instance->maxDis + 1);
	int lb = (int)(accumulateDistance.back() / instance->maxDis);
	int* chosenPos = new int[length];
	double finalfit = DBL_MAX;
	for (int i = lb; i <= ub; i++) {
		tryACertainNArray(0, i, chosenPos, finalfit, i, route, length, accumulateDistance, instance);
	}
	delete[] chosenPos;
	if (finalfit != DBL_MAX) {
		return finalfit;
	}
	else {
		return -1;
	}
}

void tryACertainNArray(int mlen, int nlen, int* chosenPos, double& finalfit, int curub, int* route, int length, vector<double>& accumulateDis, Case* instance) {
	for (int i = mlen; i <= length - 1 - nlen; i++) {
		if (curub == nlen) {
			double onedis = instance->distances[route[i]][instance->bestStation[route[i]][route[i + 1]]];
			if (accumulateDis[i] + onedis > instance->maxDis) {
				break;
			}
		}
		else {
			int lastpos = chosenPos[curub - nlen - 1];
			double onedis = instance->distances[route[lastpos + 1]][instance->bestStation[route[lastpos]][route[lastpos + 1]]];
			double twodis = instance->distances[route[i]][instance->bestStation[route[i]][route[i + 1]]];
			if (accumulateDis[i] - accumulateDis[lastpos + 1] + onedis + twodis > instance->maxDis) {
				break;
			}
		}
		if (nlen == 1) {
			double onedis = accumulateDis.back() - accumulateDis[i + 1] + instance->distances[instance->bestStation[route[i]][route[i + 1]]][route[i + 1]];
			if (onedis > instance->maxDis) {
				continue;
			}
		}

		chosenPos[curub - nlen] = i;
		if (nlen > 1) {
			tryACertainNArray(i + 1, nlen - 1, chosenPos, finalfit, curub, route, length, accumulateDis, instance);
		}
		else {
			double disum = accumulateDis.back();
			for (int j = 0; j < curub; j++) {
				int firstnode = route[chosenPos[j]];
				int secondnode = route[chosenPos[j] + 1];
				int thestation = instance->bestStation[firstnode][secondnode];
				disum -= instance->distances[firstnode][secondnode];
				disum += instance->distances[firstnode][thestation];
				disum += instance->distances[secondnode][thestation];
			}
			if (disum < finalfit) {
				finalfit = disum;
			}
		}
	}
}

double insertStationByRemoveArray(int* route, int length, Case* instance) {
	list<pair<int, int>> stationInserted;
	for (int i = 0; i < length - 1; i++) {
		double allowedDis = instance->maxDis;
		if (i != 0) {
			allowedDis = instance->maxDis - instance->distances[stationInserted.back().second][route[i]];
		}
		int onestation = instance->findNearestStationFeasible(route[i], route[i + 1], allowedDis);
		if (onestation == -1) return -1;
		stationInserted.push_back(make_pair(i, onestation));
	}
	while (!stationInserted.empty())
	{
		bool change = false;
		list<pair<int, int>>::iterator delone = stationInserted.begin();
		double savedis = 0;
		list<pair<int, int>>::iterator itr = stationInserted.begin();
		list<pair<int, int>>::iterator next = itr;
		next++;
		if (next != stationInserted.end()) {
			int endInd = next->first;
			int endstation = next->second;
			double sumdis = 0;
			for (int i = 0; i < endInd; i++) {
				sumdis += instance->distances[route[i]][route[i + 1]];
			}
			sumdis += instance->distances[route[endInd]][endstation];
			if (sumdis <= instance->maxDis) {
				savedis = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
					- instance->distances[route[itr->first]][route[itr->first + 1]];
			}
		}
		else {
			double sumdis = 0;
			for (int i = 0; i < length - 1; i++) {
				sumdis += instance->distances[route[i]][route[i + 1]];
			}
			if (sumdis <= instance->maxDis) {
				savedis = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
					- instance->distances[route[itr->first]][route[itr->first + 1]];
			}
		}
		itr++;
		while (itr != stationInserted.end())
		{
			int startInd, endInd;
			next = itr;
			next++;
			list<pair<int, int>>::iterator prev = itr;
			prev--;
			double sumdis = 0;
			if (next != stationInserted.end()) {
				startInd = prev->first + 1;
				endInd = next->first;
				sumdis += instance->distances[prev->second][route[startInd]];
				for (int i = startInd; i < endInd; i++) {
					sumdis += instance->distances[route[i]][route[i + 1]];
				}
				sumdis += instance->distances[route[endInd]][next->second];
				if (sumdis <= instance->maxDis) {
					double savedistemp = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
						- instance->distances[route[itr->first]][route[itr->first + 1]];
					if (savedistemp > savedis) {
						savedis = savedistemp;
						delone = itr;
					}
				}
			}
			else {
				startInd = prev->first + 1;
				sumdis += instance->distances[prev->second][route[startInd]];
				for (int i = startInd; i < length - 1; i++) {
					sumdis += instance->distances[route[i]][route[i + 1]];
				}
				if (sumdis <= instance->maxDis) {
					double savedistemp = instance->distances[route[itr->first]][itr->second] + instance->distances[itr->second][route[itr->first + 1]]
						- instance->distances[route[itr->first]][route[itr->first + 1]];
					if (savedistemp > savedis) {
						savedis = savedistemp;
						delone = itr;
					}
				}
			}
			itr++;
		}
		if (savedis != 0) {
			stationInserted.erase(delone);
			change = true;
		}
		if (!change) {
			break;
		}
	}
	double sum = 0;
	for (int i = 0; i < length - 1; i++) {
		sum += instance->distances[route[i]][route[i + 1]];
	}
	for (auto& e : stationInserted) {
		int pos = e.first;
		int stat = e.second;
		sum -= instance->distances[route[pos]][route[pos + 1]];
		sum += instance->distances[route[pos]][stat];
		sum += instance->distances[stat][route[pos + 1]];
	}
	return sum;
}

bool opt2withStation(int* route, int length, Case* instance) {
	if (length <= 4) return false;

	vector<int> tt(route, route + length);
	pair<vector<int>, double> xx = insertStationByEnumeration(tt, instance);
	double minchange = 0;
	bool flag = false;
	int* tempRoute = new int[length];
	do {
		minchange = 0;
		int mini = 0, minj = 0;
		for (int i = 0; i < length - 3; i++) {
			for (int j = i + 2; j < length - 1; j++) {
				
			}
		}
	} while (minchange < 0);
	delete[] tempRoute;
	return flag;
}