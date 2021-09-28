#include "utilities.h"

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

bool opt2noStation2(int* route, int length, Case* instance) {
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

bool orNoStation(int* route, int length, Case* instance) {
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
	/*for (int i = 1; i < (int)x.size(); i++) {
		vector<int> onetemp;
		allroutes.push_back(onetemp);
	}*/
	//int t = 0;
	int j = x.size() - 1;
	while (true) {
		int i = pp[j];
		vector<int> temp(x.begin() + i + 1, x.begin() + j + 1);
		// for (int k = i + 1; k <= j; k++) {
		// 	allroutes[t].push_back(x[k]);
		// }
		allroutes.push_back(temp);
		j = i;
		//t++;
		if (i == 0) {
			break;
		}
	}
	// vector<vector<int>> newroutes;
	// for (auto e : allroutes) {
	// 	if (!e.empty()) {
	// 		newroutes.push_back(e);
	// 	}
	// }
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
