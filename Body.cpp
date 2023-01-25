#include "Body.h"

using namespace std;

Body::Body(){}

Body::Body(double m, double v[], double x[]){
    mass = m;
    coord[0]=x[0];
    coord[1]=x[1];
    coord[2]=x[2];
    vel[0]=v[0];
    vel[1]=v[1];
    vel[2]=v[2];
}

Body::~Body(){

}

double Body::GetMass(){
    return mass;
}
void Body::GetCoordinates(double x[]){
    x[0]=coord[0];
    x[1]=coord[1];
    x[2]=coord[2];
}
void Body::GetVelocity(double v[]){
    v[0]=vel[0];
    v[1]=vel[1];
    v[2]=vel[2];
}
void Body::SetMass(double m){
    mass = m;
}
void Body::SetCoordinates(double x[]){
    coord[0]=x[0];
    coord[1]=x[1];
    coord[2]=x[2];
}
void Body::SetVelocity(double v[]){
    vel[0]=v[0];
    vel[1]=v[1];
    vel[2]=v[2];
}

bool Body::IsInCube(double c[], double l){
    bool in = false;
    double x_min = min(c[0]+l/2,c[0]-l/2);
    double x_max = max(c[0]+l/2,c[0]-l/2);
    double y_min = min(c[1]+l/2,c[1]-l/2);
    double y_max = max(c[1]+l/2,c[1]-l/2);
    double z_min = min(c[2]+l/2,c[2]-l/2);
    double z_max = max(c[2]+l/2,c[2]-l/2);
    if(coord[0]>=x_min && coord[0]<x_max &&
       coord[1]>=y_min && coord[1]<y_max &&
       coord[2]>=z_min && coord[2]<z_max){
        in = true;
    }
    return in;
}

