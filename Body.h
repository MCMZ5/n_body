#ifndef BODY_H
#define BODY_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

//Classe Body per stivare i dati di ognuno degli N corpi
class Body{
    public:
    Body();
    Body(double, double *x, double *v);
    ~Body();

    double GetMass();
    void GetCoordinates(double *x);
    void GetVelocity(double *v);

    void SetMass(double);
    void SetCoordinates(double *x);
    void SetVelocity(double *v);

    bool IsInCube(double *c, double l);

    private:
    double mass;
    double vel[3];
    double coord[3];
};

#endif