#ifndef BARNES_H
#define BARNES_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include "Body.h"
#include "Global.hpp"

using namespace std;

//Classe Barnes per definire i nodi dell'octree in cui è suddiviso lo spazio 
class Barnes{
    public:
        Barnes(); //costruttore senza dichiarazione
        Barnes(double m, double *x); //costruttore base
        Barnes(Body &b); //costruttore che appende un body preesistente
        Barnes(int &N, double *c, double L, int lvl); //costruttore per l'octree
        ~Barnes();

        void Delete_Subnodes(int lvl);

        Body cm;
        vector<Barnes*> subnodes; //sottonodi dell'octree
        bool empty;
};

#endif