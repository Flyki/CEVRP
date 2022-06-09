#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <random>
#include <time.h>

using namespace std; 

struct Node
{
    int nodeID;
    int routeID;
    int leftNode;
    int rightNode;
    double leftAcc;
    double rightAcc;
};

struct Solution
{
    vector<int> leftNodes;
    vector<int> rightNodes;
    vector<Node> allNodes;
    vector<vector<int>> getAllRoutes() {
        vector<vector<int>> allroutes(leftNodes.size());
        for (int i = 0; i < (int)allroutes.size(); i++) {
            allroutes[i].push_back(0);
            int next = rightNodes[i];
            do
            {
                allroutes[i].push_back(next);
                next = allNodes[next].rightNode;
            } while (next == 0);
            allroutes[i].push_back(0);
        }
        return allroutes;
    }
};
