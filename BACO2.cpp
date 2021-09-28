#include "BACO2.h"

BACO2::BACO2(Case* instance, int seed) {
	this->instance = instance;
	this->cdnumber = instance->depotNumber + instance->customerNumber;
	this->antno = this->cdnumber;
	this->ants = new int* [antno];
	for (int i = 0; i < antno; i++) {
		this->ants[i] = new int[cdnumber];
	}
	this->fit = new double[antno];
	this->gbest = new int[cdnumber];
	this->pher = new double* [cdnumber];
	for (int i = 0; i < cdnumber; i++) {
		this->pher[i] = new double[cdnumber];
	}
	vector<int> remain;
	for (int i = 0; i < cdnumber; i++) {
		remain.push_back(i);
	}
	remain = convexHull(remain, instance);
	for (int i = 0; i < cdnumber; i++) {
		gbest[i] = remain[i];
	}
	vector<vector<int>> solution = prinsSplit(remain, instance);
	for (int i = 0; i < (int)solution.size(); i++) {
		if (solution[i][0] != 0) {
			solution[i].insert(solution[i].begin(), 0);
		}
		if (solution[i].back() != 0) {
			solution[i].push_back(0);
		}
		opt2noStation(solution[i], instance);
	}
	this->gbestf = 0;
	for (int i = 0; i < (int)solution.size(); i++) {
		pair<vector<int>, double> xx = insertStationByRemove(solution[i], instance);
		solution[i] = xx.first;
		this->gbestf += xx.second;
	}
	for (auto a : solution) {
		for (auto b : a) {
			bestSolution.push_back(b);
		}
		bestSolution.pop_back();
	}
	bestSolution.push_back(0);
	this->t0 = 1.0 / (0.02 * gbestf);
	this->maxt = t0;
	double px = exp(log(0.05) / cdnumber);
	this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
	this->mint = this->maxt * this->mint;
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] = t0;
			if (i == j) pher[i][j] = 0;
		}
	}

	this->ibest = 0;
	stringstream ss;
	ss << instance->ID << "." << seed << ".BACO2.result.txt";
	string filename;
	ss >> filename;
	ss.clear();
	result.open(filename, ios::app);
	ss << instance->ID << "." << seed << ".BACO2.solution.txt";
	string sofilename;
	ss >> sofilename;
	ss.clear();
	sofile.open(sofilename, ios::app);
	default_random_engine gent(seed);
	gen = gent;
	uniform_real_distribution<double> udist(0.0, 1.0);
	udis = udist;

	this->usedFes = 0;
	this->candinumber = 20;
	this->candidatelist = instance->candidatelist;
	t1 = clock();
}

BACO2::~BACO2() {
	for (int i = 0; i < antno; i++) {
		delete[] ants[i];
	}
	delete[] ants;;
	for (int i = 0; i < cdnumber; i++) {
		delete[] pher[i];
	}
	delete[] pher;
	delete[] fit;
	delete[] gbest;
	result.close();
	sofile.close();
}

void BACO2::run() {
	long timelimited = (cdnumber + instance->stationNumber) * 36;
	long timeused = 0;
	while (timeused < timelimited)//(usedFes < MAXFES)
	{
		buildSolutionsByCL();
		evaluateAndUpdatePher();
		result << usedFes << ',' << gbestf << endl;
		//gettimeofday(&t2, NULL);
		t2 = clock();
		//timeused = t2.tv_sec - t1.tv_sec;
		timeused = (t2 - t1) / CLOCKS_PER_SEC;
	}
	sofile << fixed << setprecision(8) << gbestf << endl;
	for (auto e : bestSolution) {
		sofile << e << ',';
	}
	sofile << endl;
}

pair<vector<int>, double> BACO2::interpretACircleVec(int* circle) {
	vector<int> newdcirle;
	int zeropos = 0;
	for (int i = 0; i < cdnumber; i++) {
		if (circle[i] == 0) {
			zeropos = i;
			break;
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		if (zeropos == cdnumber)
			zeropos = 0;
		newdcirle.push_back(circle[zeropos]);
		zeropos++;
	}
	vector<vector<int>> solution = prinsSplit(newdcirle, instance);
	for (int i = 0; i < (int)solution.size(); i++) {
		if (solution[i][0] != 0) {
			solution[i].insert(solution[i].begin(), 0);
		}
		if (solution[i].back() != 0) {
			solution[i].push_back(0);
		}
		opt2noStation(solution[i], instance);
	}
	double onef = 0;
	for (int j = 0; j < (int)solution.size(); j++) {
		pair<vector<int>, double> xx = insertStationByRemove(solution[j], instance);
		/*if (xx.second < 0) {
			xx = insertStationByMinDeficit(solution[j], instance);
		}*/
		solution[j] = xx.first;
		if (xx.second < 0) {
			onef += 9999999;
		}
		else {
			onef += xx.second;
		}
	}
	newdcirle.clear();
	for (auto a : solution) {
		for (auto b : a) {
			newdcirle.push_back(b);
		}
		newdcirle.pop_back();
	}
	newdcirle.push_back(0);
	return make_pair(newdcirle, onef);
}

void BACO2::evaluateAndUpdatePher() {
	ibest = 0;
	for (int i = 0; i < antno; i++) {
		pair<vector<int>, double> xx = interpretACircleVec(ants[i]);
		usedFes++;
		fit[i] = xx.second;
		if (fit[i] < fit[ibest]) {
			ibest = i;
			if (fit[ibest] < gbestf) {
				bestSolution = xx.first;
				memcpy(gbest, ants[ibest], sizeof(int) * cdnumber);
				gbestf = fit[ibest];
				this->maxt = 1.0 / (0.02 * gbestf);
				double px = exp(log(0.05) / cdnumber);
				this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
				this->mint = this->maxt * this->mint;
			}
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 0; i < cdnumber - 1; i++) {
		int cur = ants[ibest][i];
		int nex = ants[ibest][i + 1];
		pher[cur][nex] += (1.0 / fit[ibest]);
		pher[nex][cur] = pher[cur][nex];
	}
	int cur = ants[ibest][cdnumber - 1];
	int nex = ants[ibest][0];
	pher[cur][nex] += (1.0 / fit[ibest]);
	pher[nex][cur] = pher[cur][nex];
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			if (i == j) continue;
			if (pher[i][j] < mint) pher[i][j] = mint;
			if (pher[i][j] > maxt) pher[i][j] = maxt;
		}
	}
}

void BACO2::buildSolutionsByCL() {
	int* allpoints = new int[cdnumber];
	for (int i = 0; i < cdnumber; i++) {
		allpoints[i] = i;
	}
	double* prob = new double[cdnumber];
	int* havebeenchoosen = new int[cdnumber];
	int* alltemp = new int[cdnumber];
	for (int i = 0; i < antno; i++) {
		memcpy(alltemp, allpoints, sizeof(int) * cdnumber);
		int alllen = cdnumber;
		int curone = 0;
		int counter = 0;
		ants[i][counter] = alltemp[curone];
		memcpy(&alltemp[curone], &alltemp[curone + 1], sizeof(int) * (alllen - curone - 1));
		alllen--;
		memset(havebeenchoosen, 0, sizeof(int) * cdnumber);
		havebeenchoosen[ants[i][counter]] = 1;
		while (alllen > 0)
		{
			counter++;
			int lastone = ants[i][counter - 1];
			vector<int> tobechosen;
			for (auto e : candidatelist[lastone]) {
				if (havebeenchoosen[e] == 0) {
					tobechosen.push_back(e);
				}
			}
			if (tobechosen.empty()) {
				//tobechosen = alltemp;
				tobechosen.insert(tobechosen.begin(), alltemp, alltemp + alllen);
			}
			memset(prob, 0, sizeof(double) * cdnumber);
			for (int j = 0; j < (int)tobechosen.size(); j++) {
				prob[j] = pher[lastone][tobechosen[j]] * (1.0 / instance->distances[lastone][tobechosen[j]]) * (1.0 / instance->distances[lastone][tobechosen[j]]);
			}
			for (int j = 1; j < (int)tobechosen.size(); j++) {
				prob[j] += prob[j - 1];
			}
			for (int j = 0; j < (int)tobechosen.size(); j++) {
				prob[j] /= prob[tobechosen.size() - 1];
			}
			double threshold = udis(gen);
			int theind = 0;
			for (int j = 0; j < (int)tobechosen.size(); j++) {
				if (prob[j] >= threshold || j == (int)tobechosen.size() - 1) {
					theind = j;
					break;
				}
			}
			ants[i][counter] = tobechosen[theind];
			havebeenchoosen[ants[i][counter]] = 1;
			int erapos = 0;
			for (int j = 0; j < alllen; j++) {
				if (alltemp[j] == ants[i][counter]) {
					erapos = j;
					break;
				}
			}
			memcpy(&alltemp[erapos], &alltemp[erapos + 1], sizeof(int) * (alllen - erapos - 1));
			alllen--;
		}
	}
	delete[] allpoints;
	delete[] alltemp;
	delete[] havebeenchoosen;
	delete[] prob;
}

void BACO2::buildSolutionsByAll() {
	vector<int> allpoints;
	for (int i = 0; i < cdnumber; i++) {
		allpoints.push_back(i);
	}
	double* prob = new double[cdnumber];
	for (int i = 0; i < antno; i++) {
		vector<int> alltemp = allpoints;
		int curone = 0;
		int counter = 0;
		ants[i][counter] = alltemp[curone];
		alltemp.erase(alltemp.begin() + curone);
		while (!alltemp.empty())
		{
			counter++;
			int lastone = ants[i][counter - 1];
			memset(prob, 0, sizeof(double) * cdnumber);
			for (int j = 0; j < (int)alltemp.size(); j++) {
				prob[j] = pher[lastone][alltemp[j]] * (1.0 / instance->distances[lastone][alltemp[j]]) * (1.0 / instance->distances[lastone][alltemp[j]]);
			}
			for (int j = 1; j < (int)alltemp.size(); j++) {
				prob[j] += prob[j - 1];
			}
			for (int j = 0; j < (int)alltemp.size(); j++) {
				prob[j] /= prob[alltemp.size() - 1];
			}
			double threshold = udis(gen);
			int theind = 0;
			for (int j = 0; j < (int)alltemp.size(); j++) {
				if (prob[j] >= threshold || j == (int)alltemp.size() - 1) {
					theind = j;
					break;
				}
			}
			ants[i][counter] = alltemp[theind];
			alltemp.erase(alltemp.begin() + theind);
		}
	}
	delete[] prob;
}