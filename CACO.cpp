#include "CACO.h"


CACO::CACO(Case* instance, int seed, int isCan, int isRA, int representation, double timer, double afr) {
    this->instance = instance;
    this->isCan = isCan;
    this->isRA = isRA;
	this->representation = representation;
    this->usedFes = 0;
    this->refined = 0;
    this->repaired = 0;
    this->timerate = timer;
    this->cdnumber = instance->customerNumber + instance->depotNumber;
	this->antno = this->cdnumber;
	ants.reserve(this->antno);
	for (int i = 0; i < this->antno; i++) {
		this->ants.push_back(new Ant(instance->vehicleNumber * 2, cdnumber));
	}
    bestSolution = new Ant(instance->vehicleNumber * 2, cdnumber);
	this->pher = new double*[cdnumber];
	for (int i = 0; i < cdnumber; i++) {
		this->pher[i] = new double[cdnumber];
	}
	if (representation == 1) {
		vector<int> remain;
		for (int i = 0; i < cdnumber; i++) {
			remain.push_back(i);
		}
		remain = convexHull(remain, instance);
		for (int i = 0; i < cdnumber; i++) {
			bestSolution->circle[i] = remain[i];
		}
		prinsSplitAnt(bestSolution, instance);
		opt2ToAnt(bestSolution, instance);
		bestSolution->fit = fixOneSolution(bestSolution);
	}
	else {
		generateASolutionGreedy(bestSolution);
		opt2ToAnt(bestSolution, instance);
		opt2starNoStation2(bestSolution, instance);
		orToAnt(bestSolution, instance);
		bestSolution->fit = fixOneSolution(bestSolution);
	}
    this->t0 = 1.0 / (0.02 * bestSolution->fit);
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

    string instanceName = instance->filename.substr(0, instance->filename.find_last_of('.'));
    string directoryPath = STATS_PATH + "/" + instanceName + "/" + to_string(seed);
    create_directories_if_not_exists(directoryPath);
    string filename = "evols." + to_string(isCan) + "." + to_string(isRA) + "." + to_string(representation) + "." + instanceName + ".csv";
    result.open(directoryPath + "/" + filename);
    result << "accumulated_ant_num" << "," << "upper_size" << "," << "lower_size" << "," << "time_used" << "," << "evals" << "," << "progress" << "," << "min_fit" << endl;

//    string solDirectoryPath = STATS_PATH + "/" + instanceName;
//    create_directories_if_not_exists(solDirectoryPath);
    string sofilename = "solution." + to_string(isCan) + "." + to_string(isRA) + "." + to_string(representation) + "." + instanceName + ".txt";
    sofile.open(directoryPath + "/" + sofilename);

    default_random_engine gent(seed);
	this->gen = gent;
	uniform_real_distribution<double> udist(0.0, 1.0);
	this->udis = udist;
	normal_distribution<double> ndist(0, 1);
	this->ndis = ndist;

    if (isCan == 1) {
		this->candidatelist = instance->candidatelist;
	}
    this->maxRefine = BIGVALUE;
    this->minRepair = 0;
	this->staTime = clock();
	this->afr = afr;
}

void CACO::generateASolutionGreedy(Ant* anant) {
    vector<int> remain(cdnumber - 1);
	for (int i = 1; i < cdnumber; i++) {
		remain[i - 1] = i;
	}
	anant->reset();
	double loadofroute = 0;
	while (!remain.empty())
	{
		int routeindex = anant->routeNum;
		if (anant->nodeNum[routeindex] == 0) {
			anant->route[routeindex][0] = 0;
			anant->nodeNum[routeindex]++;
		}
		else {
			int lastone = anant->route[routeindex][anant->nodeNum[routeindex] - 1];
			int chosenIndex = -1;
			double mindis = 999999999;
			for (int i = 0; i < (int)remain.size(); i++) {
				if (instance->demand[remain[i]] + loadofroute <= instance->maxC &&
					instance->getDistance(lastone, remain[i]) < mindis) {
					mindis = instance->getDistance(lastone, remain[i]);
					chosenIndex = i;
				}
			}
			int nextnode = 0;
			if (chosenIndex != -1) {
				nextnode = remain[chosenIndex];
				remain.erase(remain.begin() + chosenIndex);
			}
			loadofroute += instance->demand[nextnode];
			anant->route[routeindex][anant->nodeNum[routeindex]] = nextnode;
			anant->nodeNum[routeindex]++;
			if (nextnode == 0) {
				anant->routeNum++;
				anant->demsum[routeindex] = loadofroute;
				loadofroute = 0;
			}
		}
	}
	anant->route[anant->routeNum][anant->nodeNum[anant->routeNum]] = 0;
	anant->nodeNum[anant->routeNum]++;
	anant->demsum[anant->routeNum] = loadofroute;
	anant->routeNum++;

    anant->fit = 0;
    for (int i = 0; i < anant->routeNum; i++) {
        for (int j = 0; j < anant->nodeNum[i] - 1; j++) {
            anant->fit += instance->getDistance(anant->route[i][j], anant->route[i][j + 1]);
        }
    }
}

CACO::~CACO() {
    delete bestSolution;
    for (int i = (int)ants.size() - 1; i >= 0; i--) {
        delete ants[i];
    }
    for (int i = 0; i < cdnumber; i++) {
        delete[] pher[i];
    }
    delete[] pher;
    result.close();
    sofile.close();
}

void CACO::run() {
    double timelimited = (cdnumber + instance->stationNumber) * this->timerate * 60 * 60; // seconds
    double timeused = 0;
	int generationNum = 0;
    while (timeused < timelimited)
//    while (generationNum < 10000)
//    while (true)
	{
		generationNum++;
		if (representation != 1) {
			if (isCan == 1) {
				buildSolutionsFromCandi();
			}
			else {
				buildSolutions();
			}
			if (isRA == 1) {
				//evaluateSome();
				evaluateSomeForOnlyFixing0();
			}
			else {
				evaluateAll();
			}
		}
		else {
			if (isCan == 1) {
				buildSolutionsFromCandi2();
			}
			else {
				buildSolutions2();
			}
			if (isRA == 1) {
				//evaluateSome2();
				evaluateSomeForOnlyFixing1();
			}
			else {
				evaluateAll2();
			}
		}
        
        usedFes += antno;
        endTime = clock();
        timeused = static_cast<double>(endTime - staTime) / CLOCKS_PER_SEC; // seconds
        double evals = instance->getEvals();
        result << usedFes << ',' << refined << ',' << repaired << ',' << timeused << ',' << evals << "," << evals/instance->maxEvals << "," << bestSolution->fit << endl;
        // usedFes: 使用过的蚂蚁数
        // refined: 使用Confidence-based selection挑选出来做local search的蚂蚁数
        // repaired: 使用Confidence-based selection挑选出来做recharging的蚂蚁数
        // timeused: 已使用的时间
        // bestSolution->fit: 目前最好的解
//        if (instance->getEvals() > instance->maxEvals) break; //TODO:
    }
    sofile << fixed << setprecision(8) << bestSolution->fit << endl;
	for (int i = 0; i < bestSolution->routeNum; i++) {
		for (int j = 0; j < bestSolution->nodeNum[i]; j++) {
			sofile << bestSolution->route[i][j] << ',';
		}
		sofile << endl;
	}
	sofile << endl;
}

void CACO::buildSolutions() {
    vector<int> allpoints(cdnumber - 1);
	for (int i = 1; i < cdnumber; i++) {
		allpoints[i - 1] = i;
	}
	double* prob = new double[cdnumber];
	int* choices = new int[cdnumber];
	int choinum = 0;
	for (int i = 0; i < antno; i++) {
		unordered_set<int> alltemp(allpoints.begin(), allpoints.end());
		ants[i]->reset();
		double loadofOneRoute = 0;
		double remainDemand = instance->totalDem;
		while (!alltemp.empty())
		{
			int routeindex = ants[i]->routeNum;
			if (ants[i]->nodeNum[routeindex] == 0) {
				ants[i]->route[routeindex][0] = 0;
				ants[i]->nodeNum[routeindex]++;
			}
			else {
				choinum = 0;
				for (auto& e : alltemp) {
					if (instance->demand[e] <= instance->maxC - loadofOneRoute) {
						choices[choinum] = e;
						choinum++;
					}
				}
				if ((choinum == 0 || remainDemand <= instance->maxC * (instance->vehicleNumber - routeindex - 1)) && ants[i]->nodeNum[routeindex] > 1) {
					choices[choinum] = 0;
					choinum++;
				}
				int lastone = ants[i]->route[routeindex][ants[i]->nodeNum[routeindex] - 1];
				for (int j = 0; j < choinum; j++) {
					double heuinfor = (1.0 / instance->getDistance(lastone, choices[j])) * (1.0 / instance->getDistance(lastone, choices[j]));
					prob[j] = pher[lastone][choices[j]] * heuinfor;
				}
				for (int j = 1; j < choinum; j++) {
					prob[j] += prob[j - 1];
				}
				double threshold = udis(gen) * prob[choinum - 1];
				int theind = 0;
				for (int j = 0; j < choinum; j++) {
					if (prob[j] >= threshold || j == choinum - 1) {
						theind = j;
						break;
					}
				}
				int nextnode = choices[theind];
				loadofOneRoute += instance->demand[nextnode];
				ants[i]->route[routeindex][ants[i]->nodeNum[routeindex]] = nextnode;
				ants[i]->nodeNum[routeindex]++;
				remainDemand -= instance->demand[nextnode];
				if (nextnode == 0) {
					ants[i]->routeNum++;
					ants[i]->demsum[routeindex] = loadofOneRoute;
					loadofOneRoute = 0;
				}
				else {
					alltemp.erase(nextnode);
				}
			}
		}
		ants[i]->route[ants[i]->routeNum][ants[i]->nodeNum[ants[i]->routeNum]] = 0;
		ants[i]->nodeNum[ants[i]->routeNum]++;
		ants[i]->demsum[ants[i]->routeNum] = loadofOneRoute;
		ants[i]->routeNum++;
	}
	delete[] prob;
	delete[] choices;
}

void CACO::buildSolutionsFromCandi() {
    bool* checked = new bool[cdnumber];
	double* prob = new double[cdnumber];
	int* choices = new int[cdnumber];
	int choinum = 0;
	for (int i = 0; i < antno; i++) {
		ants[i]->reset();
		int remain = cdnumber - 1;
		memset(checked, 0, sizeof(bool) * cdnumber);
		double loadofOneRoute = 0;
		double remainDemand = instance->totalDem;
		while (remain > 0)
		{
			int routeindex = ants[i]->routeNum;
			if (ants[i]->nodeNum[routeindex] == 0) {
				ants[i]->route[routeindex][0] = 0;
				ants[i]->nodeNum[routeindex]++;
			}
			else {
				int lastone = ants[i]->route[routeindex][ants[i]->nodeNum[routeindex] - 1];
				choinum = 0;
				if (lastone == 0) {
					for (int j = 1; j < cdnumber; j++) {
						if (checked[j] == false) {
							choices[choinum] = j;
							choinum++;
						}
					}
				}
				else {
					//check candidate list
					for (int j = 0; j < (int)candidatelist[lastone].size(); j++) {
						if (checked[candidatelist[lastone][j]] == false && instance->demand[candidatelist[lastone][j]] <= instance->maxC - loadofOneRoute) {
							choices[choinum] = candidatelist[lastone][j];
							choinum++;
						}
					}
					//if none, check all remaining
					if (choinum == 0) {
						for (int j = 1; j < cdnumber; j++) {
							if (checked[j] == false && instance->demand[j] <= instance->maxC - loadofOneRoute) {
								choices[choinum] = j;
								choinum++;
							}
						}
					}
					//whether allow to return
					if (choinum == 0 || remainDemand <= instance->maxC * (instance->vehicleNumber - routeindex - 1)) {
						choices[choinum] = 0;
						choinum++;
					}
				}
				//roullet wheel selection
				for (int j = 0; j < choinum; j++) {
					prob[j] = pher[lastone][choices[j]] * (1.0 / instance->getDistance(lastone, choices[j])) * (1.0 / instance->getDistance(lastone, choices[j]));
				}
				for (int j = 1; j < choinum; j++) {
					prob[j] += prob[j - 1];
				}
				double threshold = udis(gen) * prob[choinum - 1];
				int theind = 0;
				for (int j = 0; j < choinum; j++) {
					if (prob[j] >= threshold || j == choinum - 1) {
						theind = j;
						break;
					}
				}
				int nextnode = choices[theind];
				loadofOneRoute += instance->demand[nextnode];
				ants[i]->route[routeindex][ants[i]->nodeNum[routeindex]] = nextnode;
				ants[i]->nodeNum[routeindex]++;
				remainDemand -= instance->demand[nextnode];
				checked[nextnode] = true;
				if (nextnode == 0) {
					ants[i]->routeNum++;
					ants[i]->demsum[routeindex] = loadofOneRoute;
					loadofOneRoute = 0;
				}
				else {
					remain--;
				}
			}
		}
		ants[i]->route[ants[i]->routeNum][ants[i]->nodeNum[ants[i]->routeNum]] = 0;
		ants[i]->nodeNum[ants[i]->routeNum]++;
		ants[i]->demsum[ants[i]->routeNum] = loadofOneRoute;
		ants[i]->routeNum++;
	}
	delete[] checked;
	delete[] prob;
	delete[] choices;
}

void CACO::evaluateAll() {
	for (int i = 0; i < antno; i++) {
        ants[i]->fit = 0;
        for (int j = 0; j < ants[i]->routeNum; j++) {
            for (int k = 0; k < ants[i]->nodeNum[j] - 1; k++) {
                ants[i]->fit += instance->getDistance(ants[i]->route[j][k], ants[i]->route[j][k + 1]);
            }
        }
        ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
	}
    for (int i = 0; i < antno; i++) {
        opt2ToAnt(ants[i], instance);
    }
    for (int i = 0; i < antno; i++) {
		opt2starNoStation2(ants[i], instance);
	}
	for (int i = 0; i < antno; i++) {
		orToAnt(ants[i], instance);
	}
	for (int i = 0; i < antno; i++) {
		ants[i]->fit = fixOneSolution(ants[i]);
	}
    ibest = 0;
	for (int i = 1; i < antno; i++) {
		if (ants[i]->fit < ants[ibest]->fit) {
			ibest = i;
		}
	}
	if (ants[ibest]->fit < bestSolution->fit) {
		bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 0; i < ants[ibest]->routeNum; i++) {
		if (ants[ibest]->nodeNum[i] <= 3) {
			int cur = ants[ibest]->route[i][0];
			int nex = ants[ibest]->route[i][1];
			pher[cur][nex] += (1.0 / ants[ibest]->fit);
			pher[nex][cur] = pher[cur][nex];
		}
		else {
			for (int j = 1; j < ants[ibest]->nodeNum[i]; j++) {
				int cur = ants[ibest]->route[i][j];
				int nex = ants[ibest]->route[i][j - 1];
				pher[cur][nex] += (1.0 / ants[ibest]->fit);
				pher[nex][cur] = pher[cur][nex];
			}
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
	this->refined += antno;
	this->repaired += antno;
}

void CACO::evaluateSome() {
    int examplar1 = 0;
    double minFitBeforeLS = DBL_MAX;
    for (int i = 0; i < antno; i++) {
        ants[i]->fit = 0;
        for (int j = 0; j < ants[i]->routeNum; j++) {
            for (int k = 0; k < ants[i]->nodeNum[j] - 1; k++) {
                ants[i]->fit += instance->getDistance(ants[i]->route[j][k], ants[i]->route[j][k + 1]);
            }
        }
        ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
        if (minFitBeforeLS > ants[i]->fit) {
            examplar1 = i;
            minFitBeforeLS = ants[i]->fit;
        }
    }
    //screen the promising solutions to do local search
    unordered_set<double> uniq1;
    uniq1.insert(ants[examplar1]->fit);
    opt2ToAnt(ants[examplar1], instance);
    opt2starNoStation2(ants[examplar1], instance);
    orToAnt(ants[examplar1], instance);
    ants[examplar1]->fit = round(ants[examplar1]->fit * 1000000.0) / 1000000.0;
    if (refiningImprovement.size() < 1) {
        this->maxRefine = BIGVALUE;
    }
    else {
        this->maxRefine = getLargest(refiningImprovement);
    }
    if (this->maxRefine < (minFitBeforeLS - ants[examplar1]->fit)) {
        this->maxRefine = (minFitBeforeLS - ants[examplar1]->fit);
    }
    vector<int> screenlist1;
    screenlist1.reserve(antno);
    for (int i = 0; i < antno; i++) {
        if (i != examplar1 && ants[i]->fit - this->maxRefine <= ants[examplar1]->fit && uniq1.find(ants[i]->fit) == uniq1.end()) {
            screenlist1.push_back(i);
            uniq1.insert(ants[i]->fit);
        }
    }
    //do local search on the promising solutions
    double iterMaxRef = minFitBeforeLS - ants[examplar1]->fit;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        double orifit = ants[ind]->fit;
        opt2ToAnt(ants[ind], instance);
        opt2starNoStation2(ants[ind], instance);
        orToAnt(ants[ind], instance);
        if (iterMaxRef < orifit - ants[ind]->fit) {
            iterMaxRef = orifit - ants[ind]->fit;
        }
        ants[ind]->fit = round(ants[ind]->fit * 1000000.0) / 1000000.0;
    }
    refiningImprovement.push_back(iterMaxRef);
    if (refiningImprovement.size() > 1) refiningImprovement.erase(refiningImprovement.begin());
    screenlist1.push_back(examplar1);
    //screen the promising solutions to fix
    int examplar2 = 0;
    double minFitBeforeFix = DBL_MAX;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        if (ants[ind]->fit < minFitBeforeFix) {
            minFitBeforeFix = ants[ind]->fit;
            examplar2 = ind;
        }
    }
    unordered_set<double> uniq2;
    uniq2.insert(ants[examplar2]->fit);
    ants[examplar2]->fit = fixOneSolution(ants[examplar2]);
    if (this->minRepair >= ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = (ants[examplar2]->fit - minFitBeforeFix) * 0.8;
    }
    vector<int> screenlist2;
    screenlist2.reserve(antno);
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        if (ind != examplar2 && ants[ind]->fit + this->minRepair <= ants[examplar2]->fit && uniq2.find(ants[ind]->fit) == uniq2.end()) {
            screenlist2.push_back(ind);
            uniq2.insert(ants[ind]->fit);
        }
    }
    //fix the selected solutions
    for (int i = 0; i < (int)screenlist2.size(); i++) {
        int ind = screenlist2[i];
        double orifit = ants[ind]->fit;
        ants[ind]->fit = fixOneSolution(ants[ind]);
        if (ants[ind]->fit - orifit < this->minRepair || this->minRepair == 0) {
            this->minRepair = ants[ind]->fit - orifit;
        }
    }
    if (this->minRepair > ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = ants[examplar2]->fit - minFitBeforeFix;
    }
    screenlist2.push_back(examplar2);
    //find the iteration best
    this->ibest = screenlist2[0];
    for (int i = 1; i < (int)screenlist2.size(); i++) {
        if (ants[this->ibest]->fit > ants[screenlist2[i]]->fit) {
            this->ibest = screenlist2[i];
        }
    }
    //update
    if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
    for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 0; i < ants[ibest]->routeNum; i++) {
		if (ants[ibest]->nodeNum[i] <= 3) {
			int cur = ants[ibest]->route[i][0];
			int nex = ants[ibest]->route[i][1];
			pher[cur][nex] += (1.0 / ants[ibest]->fit);
			pher[nex][cur] = pher[cur][nex];
		}
		else {
			for (int j = 1; j < ants[ibest]->nodeNum[i]; j++) {
				int cur = ants[ibest]->route[i][j];
				int nex = ants[ibest]->route[i][j - 1];
				pher[cur][nex] += (1.0 / ants[ibest]->fit);
				pher[nex][cur] = pher[cur][nex];
			}
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
    this->refined += screenlist1.size();
    this->repaired += screenlist2.size();
}

void CACO::buildSolutions2() {
	vector<int> allpoints;
	for (int i = 0; i < cdnumber; i++) {
		allpoints.push_back(i);
	}
	double* prob = new double[cdnumber];
	for (int i = 0; i < antno; i++) {
		vector<int> alltemp = allpoints;
		int curone = (int)(udis(gen) * alltemp.size());
		int counter = 0;
		while (alltemp[curone] == 0)
		{
			curone = (int)(udis(gen) * alltemp.size());
		}
		ants[i]->circle[counter] = alltemp[curone];
		alltemp.erase(alltemp.begin() + curone);
		while (!alltemp.empty())
		{
			counter++;
			int lastone = ants[i]->circle[counter - 1];
			memset(prob, 0, sizeof(double) * cdnumber);
			for (int j = 0; j < (int)alltemp.size(); j++) {
				double heu1 = 1.0 / instance->getDistance(lastone, alltemp[j]);
				prob[j] = pher[lastone][alltemp[j]] * pow(heu1, 2.0);
			}
			for (int j = 1; j < (int)alltemp.size(); j++) {
				prob[j] += prob[j - 1];
			}
			double threshold = udis(gen) * prob[alltemp.size() - 1];
			int theind = 0;
			for (int j = 0; j < (int)alltemp.size(); j++) {
				if (prob[j] >= threshold || j == (int)alltemp.size() - 1) {
					theind = j;
					break;
				}
			}
			ants[i]->circle[counter] = alltemp[theind];
			alltemp.erase(alltemp.begin() + theind);
		}
	}
	delete[] prob;

	int* temp = new int[cdnumber];
	for (int i = 0; i < antno; i++) {
		memcpy(temp, ants[i]->circle, sizeof(int) * cdnumber);
		int zeropos = 0;
		for (int j = 0; j < cdnumber; j++) {
			if (ants[i]->circle[j] == 0) {
				zeropos = j;
				break;
			}
		}
		for (int j = 0; j < cdnumber; j++) {
			if (zeropos == cdnumber) {
				zeropos = 0;
			}
			ants[i]->circle[j] = temp[zeropos];
			zeropos++;
		}
	}
	delete[] temp;


	for (int i = 0; i < antno; i++) {
		prinsSplitAnt(ants[i], instance);
	}
}

void CACO::buildSolutionsFromCandi2() {
	double* prob = new double[cdnumber];
	int* havebeenchosen = new int[cdnumber];

	for (int i = 0; i < antno; i++) {
		int alllen = cdnumber;
		int curone = (int)(udis(gen) * alllen);
		int counter = 0;
		ants[i]->circle[counter] = curone;
		alllen--;
		memset(havebeenchosen, 0, sizeof(int) * cdnumber);
		havebeenchosen[ants[i]->circle[counter]] = 1;
		while (alllen > 0)
		{
			counter++;
			int lastone = ants[i]->circle[counter - 1];
			vector<int> tobechosen(30);
			int length = 0;
			for (auto e : candidatelist[lastone]) {
				if (havebeenchosen[e] == 0) {
					tobechosen[length] = e;
					length++;
				}
			}
			if (length == 0) {
				tobechosen.resize(alllen);
				int cursor = 0;
				for (int j = 0; j < cdnumber; j++) {
					if (havebeenchosen[j] == 0) {
						tobechosen[cursor] = j;
						cursor++;
					}
				}
				length = cursor;
			}
			memset(prob, 0, sizeof(double) * cdnumber);
			for (int j = 0; j < length; j++) {
				double heu1 = 1.0 / instance->getDistance(lastone, tobechosen[j]);
				prob[j] = pher[lastone][tobechosen[j]] * pow(heu1, 2.0);
			}
			for (int j = 1; j < length; j++) {
				prob[j] += prob[j - 1];
			}
			double threshold = udis(gen) * prob[length - 1];
			int theind = 0;
			for (int j = 0; j < length; j++) {
				if (prob[j] >= threshold || j == length - 1) {
					theind = j;
					break;
				}
			}
			ants[i]->circle[counter] = tobechosen[theind];
			havebeenchosen[ants[i]->circle[counter]] = 1;
			alllen--;
		}
	}
	delete[] prob;
	delete[] havebeenchosen;

	int* temp = new int[cdnumber];
	for (int i = 0; i < antno; i++) {
		memcpy(temp, ants[i]->circle, sizeof(int) * cdnumber);
		int zeropos = 0;
		for (int j = 0; j < cdnumber; j++) {
			if (ants[i]->circle[j] == 0) {
				zeropos = j;
				break;
			}
		}
		for (int j = 0; j < cdnumber; j++) {
			if (zeropos == cdnumber) {
				zeropos = 0;
			}
			ants[i]->circle[j] = temp[zeropos];
			zeropos++;
		}
	}
	delete[] temp;

	for (int i = 0; i < antno; i++) {
		prinsSplitAnt(ants[i], instance);
	}
}

void CACO::evaluateAll2() {
	for (int i = 0; i < antno; i++) {
		opt2ToAnt(ants[i], instance);
		ants[i]->fit = fixOneSolution(ants[i]);
	}
	ibest = 0;
	for (int i = 1; i < antno; i++) {
		if (ants[i]->fit < ants[ibest]->fit) {
			ibest = i;
		}
	}
	if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->circle, ants[ibest]->circle, sizeof(int) * cdnumber);
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 1; i < cdnumber - 1; i++) {
		int cur = ants[ibest]->circle[i];
		int nex = ants[ibest]->circle[i + 1];
		pher[cur][nex] += (1.0 / ants[ibest]->fit);
		pher[nex][cur] = pher[cur][nex];
	}
	int cur = ants[ibest]->circle[cdnumber - 1];
	int nex = 0;
	pher[cur][nex] += (1.0 / ants[ibest]->fit);
	pher[nex][cur] = pher[cur][nex];
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
	this->refined += antno;
	this->repaired += antno;
}

void CACO::evaluateSome2() {
	int examplar1 = 0;
	double minFitBeforeLS = DBL_MAX;
	for (int i = 0; i < antno; i++) {
		ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
		if (minFitBeforeLS > ants[i]->fit) {
			examplar1 = i;
			minFitBeforeLS = ants[i]->fit;
		}
	}
	//screen the promising solutions to do local search
	unordered_set<double> uniq1;
	uniq1.insert(ants[examplar1]->fit);
	opt2ToAnt(ants[examplar1], instance);
	ants[examplar1]->fit = round(ants[examplar1]->fit * 1000000.0) / 1000000.0;
	if (refiningImprovement.size() < 30) {
		this->maxRefine = BIGVALUE;
	}
	else {
		this->maxRefine = getLargest(refiningImprovement);
	}
	if (this->maxRefine < (minFitBeforeLS - ants[examplar1]->fit)) {
		this->maxRefine = (minFitBeforeLS - ants[examplar1]->fit) * 1.2;
	}
	vector<int> screenlist1;
	screenlist1.reserve(antno);
	for (int i = 0; i < antno; i++) {
        if (i != examplar1 && ants[i]->fit - this->maxRefine <= ants[examplar1]->fit && uniq1.find(ants[i]->fit) == uniq1.end()) {
            screenlist1.push_back(i);
            uniq1.insert(ants[i]->fit);
        }
    }
	//do local search one the promising solutions
	double iterMaxRef = minFitBeforeLS - ants[examplar1]->fit;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        double orifit = ants[ind]->fit;
        opt2ToAnt(ants[ind], instance);
        if (iterMaxRef < orifit - ants[ind]->fit) {
            iterMaxRef = orifit - ants[ind]->fit;
        }
        ants[ind]->fit = round(ants[ind]->fit * 1000000.0) / 1000000.0;
    }
    refiningImprovement.push_back(iterMaxRef);
    if (refiningImprovement.size() > 30) refiningImprovement.erase(refiningImprovement.begin());
    screenlist1.push_back(examplar1);
	//screen the promising solutions to fix
	int examplar2 = 0;
    double minFitBeforeFix = DBL_MAX;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        if (ants[ind]->fit < minFitBeforeFix) {
            minFitBeforeFix = ants[ind]->fit;
            examplar2 = ind;
        }
    }
	unordered_set<double> uniq2;
    uniq2.insert(ants[examplar2]->fit);
    ants[examplar2]->fit = fixOneSolution(ants[examplar2]);
    if (this->minRepair >= ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = (ants[examplar2]->fit - minFitBeforeFix) * 0.8;
    }
	vector<int> screenlist2;
    screenlist2.reserve(antno);
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        if (ind != examplar2 && ants[ind]->fit + this->minRepair <= ants[examplar2]->fit && uniq2.find(ants[ind]->fit) == uniq2.end()) {
            screenlist2.push_back(ind);
            uniq2.insert(ants[ind]->fit);
        }
    }
	//fix the selected solutions
	for (int i = 0; i < (int)screenlist2.size(); i++) {
        int ind = screenlist2[i];
        double orifit = ants[ind]->fit;
        ants[ind]->fit = fixOneSolution(ants[ind]);
        if (ants[ind]->fit - orifit < this->minRepair || this->minRepair == 0) {
            this->minRepair = ants[ind]->fit - orifit;
        }
    }
    if (this->minRepair > ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = ants[examplar2]->fit - minFitBeforeFix;
    }
    screenlist2.push_back(examplar2);
	//find the iteration best
	this->ibest = screenlist2[0];
    for (int i = 1; i < (int)screenlist2.size(); i++) {
        if (ants[this->ibest]->fit > ants[screenlist2[i]]->fit) {
            this->ibest = screenlist2[i];
        }
    }
	//update
	if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->circle, ants[ibest]->circle, sizeof(int) * cdnumber);
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 1; i < cdnumber - 1; i++) {
		int cur = ants[ibest]->circle[i];
		int nex = ants[ibest]->circle[i + 1];
		pher[cur][nex] += (1.0 / ants[ibest]->fit);
		pher[nex][cur] = pher[cur][nex];
	}
	int cur = ants[ibest]->circle[cdnumber - 1];
	int nex = 0;
	pher[cur][nex] += (1.0 / ants[ibest]->fit);
	pher[nex][cur] = pher[cur][nex];
	// for (int i = 0; i < ants[ibest]->routeNum; i++) {
	// 	int cur = 0;
	// 	int nex = ants[ibest]->route[i][1];
	// 	pher[cur][nex] += (1.0 / ants[ibest]->fit);
	// 	pher[nex][cur] = pher[cur][nex];
	// 	if (ants[ibest]->nodeNum[i] > 3) {
	// 		int cur = 0;
	// 		int nex = ants[ibest]->route[i][ants[ibest]->nodeNum[i] - 2];
	// 		pher[cur][nex] += (1.0 / ants[ibest]->fit);
	// 		pher[nex][cur] = pher[cur][nex];
	// 	}
	// }
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
	this->refined += screenlist1.size();
	this->repaired += screenlist2.size();
}

void CACO::evaluateSomeForOnlyLocalSearch0() {
	int examplar1 = 0;
    double minFitBeforeLS = DBL_MAX;
    for (int i = 0; i < antno; i++) {
        ants[i]->fit = 0;
        for (int j = 0; j < ants[i]->routeNum; j++) {
            for (int k = 0; k < ants[i]->nodeNum[j] - 1; k++) {
                ants[i]->fit += instance->getDistance(ants[i]->route[j][k], ants[i]->route[j][k + 1]);
            }
        }
        ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
        if (minFitBeforeLS > ants[i]->fit) {
            examplar1 = i;
            minFitBeforeLS = ants[i]->fit;
        }
    }
    //screen the promising solutions to do local search
    unordered_set<double> uniq1;
    uniq1.insert(ants[examplar1]->fit);
    opt2ToAnt(ants[examplar1], instance);
    opt2starNoStation2(ants[examplar1], instance);
    orToAnt(ants[examplar1], instance);
    ants[examplar1]->fit = round(ants[examplar1]->fit * 1000000.0) / 1000000.0;
    if (refiningImprovement.size() < 50) {
        this->maxRefine = BIGVALUE;
    }
    else {
        this->maxRefine = getLargest(refiningImprovement);
    }
    if (this->maxRefine < (minFitBeforeLS - ants[examplar1]->fit)) {
        this->maxRefine = (minFitBeforeLS - ants[examplar1]->fit) * 1.4;
    }
    vector<int> screenlist1;
    screenlist1.reserve(antno);
    for (int i = 0; i < antno; i++) {
        if (i != examplar1 && ants[i]->fit - this->maxRefine <= ants[examplar1]->fit && uniq1.find(ants[i]->fit) == uniq1.end()) {
            screenlist1.push_back(i);
            uniq1.insert(ants[i]->fit);
        }
    }
    //do local search on the promising solutions
    double iterMaxRef = minFitBeforeLS - ants[examplar1]->fit;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        double orifit = ants[ind]->fit;
        opt2ToAnt(ants[ind], instance);
        opt2starNoStation2(ants[ind], instance);
        orToAnt(ants[ind], instance);
        if (iterMaxRef < orifit - ants[ind]->fit) {
            iterMaxRef = orifit - ants[ind]->fit;
        }
        ants[ind]->fit = round(ants[ind]->fit * 1000000.0) / 1000000.0;
    }
    refiningImprovement.push_back(iterMaxRef);
    if (refiningImprovement.size() > 50) refiningImprovement.erase(refiningImprovement.begin());
    screenlist1.push_back(examplar1);

	for (int i = 0; i < (int)screenlist1.size(); i++) {
		int ind = screenlist1[i];
		ants[ind]->fit = fixOneSolution(ants[ind]);
	}
	this->ibest = screenlist1[0];
	for (int i = 1; i < (int)screenlist1.size(); i++) {
		if (ants[this->ibest]->fit > ants[screenlist1[i]]->fit) {
			this->ibest = screenlist1[i];
		}
	}
	if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
    for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 0; i < ants[ibest]->routeNum; i++) {
		if (ants[ibest]->nodeNum[i] <= 3) {
			int cur = ants[ibest]->route[i][0];
			int nex = ants[ibest]->route[i][1];
			pher[cur][nex] += (1.0 / ants[ibest]->fit);
			pher[nex][cur] = pher[cur][nex];
		}
		else {
			for (int j = 1; j < ants[ibest]->nodeNum[i]; j++) {
				int cur = ants[ibest]->route[i][j];
				int nex = ants[ibest]->route[i][j - 1];
				pher[cur][nex] += (1.0 / ants[ibest]->fit);
				pher[nex][cur] = pher[cur][nex];
			}
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
    this->refined += screenlist1.size();
    this->repaired += screenlist1.size();
}

void CACO::evaluateSomeForOnlyLocalSearch1() {
	int examplar1 = 0;
	double minFitBeforeLS = DBL_MAX;
	for (int i = 0; i < antno; i++) {
		ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
		if (minFitBeforeLS > ants[i]->fit) {
			examplar1 = i;
			minFitBeforeLS = ants[i]->fit;
		}
	}
	//screen the promising solutions to do local search
	unordered_set<double> uniq1;
	uniq1.insert(ants[examplar1]->fit);
	opt2ToAnt(ants[examplar1], instance);
	ants[examplar1]->fit = round(ants[examplar1]->fit * 1000000.0) / 1000000.0;
	if (refiningImprovement.size() < 50) {
		this->maxRefine = BIGVALUE;
	}
	else {
		this->maxRefine = getLargest(refiningImprovement);
	}
	if (this->maxRefine < (minFitBeforeLS - ants[examplar1]->fit)) {
		this->maxRefine = (minFitBeforeLS - ants[examplar1]->fit) * 1.4;
	}
	vector<int> screenlist1;
	screenlist1.reserve(antno);
	for (int i = 0; i < antno; i++) {
        if (i != examplar1 && ants[i]->fit - this->maxRefine <= ants[examplar1]->fit && uniq1.find(ants[i]->fit) == uniq1.end()) {
            screenlist1.push_back(i);
            uniq1.insert(ants[i]->fit);
        }
    }
	//do local search one the promising solutions
	double iterMaxRef = minFitBeforeLS - ants[examplar1]->fit;
    for (int i = 0; i < (int)screenlist1.size(); i++) {
        int ind = screenlist1[i];
        double orifit = ants[ind]->fit;
        opt2ToAnt(ants[ind], instance);
        if (iterMaxRef < orifit - ants[ind]->fit) {
            iterMaxRef = orifit - ants[ind]->fit;
        }
        ants[ind]->fit = round(ants[ind]->fit * 1000000.0) / 1000000.0;
    }
    refiningImprovement.push_back(iterMaxRef);
    if (refiningImprovement.size() > 50) refiningImprovement.erase(refiningImprovement.begin());
    screenlist1.push_back(examplar1);

	for (int i = 0; i < (int)screenlist1.size(); i++) {
		int ind = screenlist1[i];
		ants[ind]->fit = fixOneSolution(ants[ind]);
	}

	this->ibest = screenlist1[0];
	for (int i = 0; i < (int)screenlist1.size(); i++) {
		if (ants[this->ibest]->fit > ants[screenlist1[i]]->fit) {
			this->ibest = screenlist1[i];
		}
	}
	if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->circle, ants[ibest]->circle, sizeof(int) * cdnumber);
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 1; i < cdnumber - 1; i++) {
		int cur = ants[ibest]->circle[i];
		int nex = ants[ibest]->circle[i + 1];
		pher[cur][nex] += (1.0 / ants[ibest]->fit);
		pher[nex][cur] = pher[cur][nex];
	}
	int cur = ants[ibest]->circle[cdnumber - 1];
	int nex = 0;
	pher[cur][nex] += (1.0 / ants[ibest]->fit);
	pher[nex][cur] = pher[cur][nex];
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
	this->refined += screenlist1.size();
	this->repaired += screenlist1.size();
}

void CACO::evaluateSomeForOnlyFixing0() {
	for (int i = 0; i < antno; i++) {
		ants[i]->fit = 0;
		for (int j = 0; j < ants[i]->routeNum; j++) {
            for (int k = 0; k < ants[i]->nodeNum[j] - 1; k++) {
                ants[i]->fit += instance->getDistance(ants[i]->route[j][k], ants[i]->route[j][k + 1]);
            }
        }
        ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
	}
	for (int i = 0; i < antno; i++) {
		opt2ToAnt(ants[i], instance);
		opt2starNoStation2(ants[i], instance);
		orToAnt(ants[i], instance);
		ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
	}
	int examplar2 = 0;
    double minFitBeforeFix = DBL_MAX;
    for (int i = 0; i < antno; i++) {
        int ind = i;
        if (ants[ind]->fit < minFitBeforeFix) {
            minFitBeforeFix = ants[ind]->fit;
            examplar2 = ind;
        }
    }
    unordered_set<double> uniq2;
    uniq2.insert(ants[examplar2]->fit);
    ants[examplar2]->fit = fixOneSolution(ants[examplar2]);
    if (this->minRepair >= ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = (ants[examplar2]->fit - minFitBeforeFix) * afr;
    }
	vector<int> screenlist2;
    screenlist2.reserve(antno);
	for (int i = 0; i < antno; i++) {
		int ind = i;
        if (ind != examplar2 && ants[ind]->fit + this->minRepair <= ants[examplar2]->fit && uniq2.find(ants[ind]->fit) == uniq2.end()) {
            screenlist2.push_back(ind);
            uniq2.insert(ants[ind]->fit);
        }
	}
	for (int i = 0; i < (int)screenlist2.size(); i++) {
        int ind = screenlist2[i];
        double orifit = ants[ind]->fit;
        ants[ind]->fit = fixOneSolution(ants[ind]);
        if (ants[ind]->fit - orifit < this->minRepair || this->minRepair == 0) {
            this->minRepair = ants[ind]->fit - orifit;
        }
    }
    if (this->minRepair > ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = ants[examplar2]->fit - minFitBeforeFix;
    }
    screenlist2.push_back(examplar2);
	this->ibest = screenlist2[0];
    for (int i = 1; i < (int)screenlist2.size(); i++) {
        if (ants[this->ibest]->fit > ants[screenlist2[i]]->fit) {
            this->ibest = screenlist2[i];
        }
    }
    //update
    if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
    for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 0; i < ants[ibest]->routeNum; i++) {
		if (ants[ibest]->nodeNum[i] <= 3) {
			int cur = ants[ibest]->route[i][0];
			int nex = ants[ibest]->route[i][1];
			pher[cur][nex] += (1.0 / ants[ibest]->fit);
			pher[nex][cur] = pher[cur][nex];
		}
		else {
			for (int j = 1; j < ants[ibest]->nodeNum[i]; j++) {
				int cur = ants[ibest]->route[i][j];
				int nex = ants[ibest]->route[i][j - 1];
				pher[cur][nex] += (1.0 / ants[ibest]->fit);
				pher[nex][cur] = pher[cur][nex];
			}
		}
	}
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
    this->refined += antno;
    this->repaired += screenlist2.size();
}

void CACO::evaluateSomeForOnlyFixing1() {
	for (int i = 0; i < antno; i++) {
		opt2ToAnt(ants[i], instance);
		ants[i]->fit = round(ants[i]->fit * 1000000.0) / 1000000.0;
	}

	int examplar2 = 0;
    double minFitBeforeFix = DBL_MAX;
    for (int i = 0; i < antno; i++) {
        int ind = i;
        if (ants[ind]->fit < minFitBeforeFix) {
            minFitBeforeFix = ants[ind]->fit;
            examplar2 = ind;
        }
    }
	unordered_set<double> uniq2;
    uniq2.insert(ants[examplar2]->fit);
    ants[examplar2]->fit = fixOneSolution(ants[examplar2]);
    if (this->minRepair >= ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = (ants[examplar2]->fit - minFitBeforeFix) * afr;
    }
	vector<int> screenlist2;
    screenlist2.reserve(antno);
    for (int i = 0; i < antno; i++) {
        int ind = i;
        if (ind != examplar2 && ants[ind]->fit + this->minRepair <= ants[examplar2]->fit && uniq2.find(ants[ind]->fit) == uniq2.end()) {
            screenlist2.push_back(ind);
            uniq2.insert(ants[ind]->fit);
        }
    }
	//fix the selected solutions
	for (int i = 0; i < (int)screenlist2.size(); i++) {
        int ind = screenlist2[i];
        double orifit = ants[ind]->fit;
        ants[ind]->fit = fixOneSolution(ants[ind]);
        if (ants[ind]->fit - orifit < this->minRepair || this->minRepair == 0) {
            this->minRepair = ants[ind]->fit - orifit;
        }
    }
    if (this->minRepair > ants[examplar2]->fit - minFitBeforeFix) {
        this->minRepair = ants[examplar2]->fit - minFitBeforeFix;
    }
    screenlist2.push_back(examplar2);
	this->ibest = screenlist2[0];
    for (int i = 1; i < (int)screenlist2.size(); i++) {
        if (ants[this->ibest]->fit > ants[screenlist2[i]]->fit) {
            this->ibest = screenlist2[i];
        }
    }
	//update
	if (ants[ibest]->fit < bestSolution->fit) {
        bestSolution->fit = ants[ibest]->fit;
		bestSolution->routeNum = ants[ibest]->routeNum;
		memcpy(bestSolution->circle, ants[ibest]->circle, sizeof(int) * cdnumber);
		memcpy(bestSolution->nodeNum, ants[ibest]->nodeNum, sizeof(int) * bestSolution->routeNum);
		for (int i = 0; i < bestSolution->routeNum; i++) {
			memcpy(bestSolution->route[i], ants[ibest]->route[i], sizeof(int) * bestSolution->nodeNum[i]);
		}
		this->maxt = 1.0 / (0.02 * bestSolution->fit);
		double px = exp(log(0.05) / cdnumber);
		this->mint = 1.0 * (1.0 - px) / (px * (((double)antno + 1) / 2.0));
		this->mint = this->maxt * this->mint;
    }
	for (int i = 0; i < cdnumber; i++) {
		for (int j = 0; j < cdnumber; j++) {
			pher[i][j] *= 0.98;
		}
	}
	for (int i = 1; i < cdnumber - 1; i++) {
		int cur = ants[ibest]->circle[i];
		int nex = ants[ibest]->circle[i + 1];
		pher[cur][nex] += (1.0 / ants[ibest]->fit);
		pher[nex][cur] = pher[cur][nex];
	}
	int cur = ants[ibest]->circle[cdnumber - 1];
	int nex = 0;
	pher[cur][nex] += (1.0 / ants[ibest]->fit);
	pher[nex][cur] = pher[cur][nex];
	for (int i = 0; i < cdnumber; i++) {
		for (int j = i + 1; j < cdnumber; j++) {
			if (pher[i][j] > maxt) pher[i][j] = maxt;
			if (pher[i][j] < mint) pher[i][j] = mint;
			pher[j][i] = pher[i][j];
		}
	}
	this->refined += antno;
	this->repaired += screenlist2.size();
}

double CACO::fixOneSolution(Ant* anant) {
    double afit = 0;
    for (int i = 0; i < anant->routeNum; i++) {
        double xx = insertStationBySimpleEnumerationArray(anant->route[i], anant->nodeNum[i], instance);
        if (xx == -1) {
            double yy = insertStationByRemoveArray(anant->route[i], anant->nodeNum[i], instance);
            if (yy == -1) {
                afit += INFEASIBLE;
            }
            else {
                afit += yy;
            }
        }
        else {
            afit += xx;
        }
    }
    afit = round(afit * 1000000.0) / 1000000.0;
	anant->fit = afit;
    return afit;
}

bool CACO::create_directories_if_not_exists(const string &directoryPath) {
    if (!fs::exists(directoryPath)) {
        try {
            fs::create_directories(directoryPath);
//            std::cout << "Directory created successfully: " << directoryPath << std::endl;
            return true;
        } catch (const std::exception& e) {
//            std::cerr << "Error creating directory: " << e.what() << std::endl;
            return false;
        }
    } else {
//        std::cout << "Directory already exists: " << directoryPath << std::endl;
        return true;
    }
}
