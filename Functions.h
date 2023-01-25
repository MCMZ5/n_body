#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "Body.h"
#include "Barnes.h"

using namespace std;

int count(string file);
vector<Body> load(string file, int dimension);
void save(string file, vector<Body> b, int dimension);
double Distance(double *a, double *b);
double Distance(Body b1, Body b2);
void Vectorial_Distance(Body b1, Body b2, double *d);
Barnes* Generate_Octree(int &N, double c[], double L);
void Delete_Octree(Barnes*);
void acc(double m, double d, double *x, double *a);
void Compute_Gravity(int i, Barnes* current, double *g, double L);

#endif