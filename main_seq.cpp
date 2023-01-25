#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "Body.h"
#include "Functions.h"
#include "Barnes.h"
#include "Global.hpp"

using namespace std;

//INTERFACCIA: VARIABILI DA MODIFICARE
#define L          30   //Dimensione dello spazio dell'octree
#define T          0.1  //Tempo totale di simulazione

//VARIABILI GLOBALI
vector<Body> bodies;

//MAIN
int main(){

//Dichiarazione variabili necessarie
  double   t     = 0.;
  double   dt    = 0.05; //salto temporale
  string input_file = "nbody_start.txt";
  string output_file = "nbody_end.txt";
  Barnes* octree;
  double center[3] = {0,0,0};

//--------------------------------------------------------------------------+
//|                     IMPORTO GLI N CORPI DEL SISTEMA                     |
//+-------------------------------------------------------------------------+ 

  int N = count(input_file);
  bodies = load(input_file,N);

//--------------------------------------------------------------------------+
//|                         EVOLUZIONE DEL SISTEMA                          |
//+-------------------------------------------------------------------------+ 

  while(t<=T){
    //Using position Verlet and Barnes-Hut algorithms 
    //Opening Angle for Barnes-Hut = 1

    cout << t << " Passo 1" << endl;
    for(int i=0;i<N;i++){ //parallelizzare
      double x[3];
      double v[3];
      bodies[i].GetCoordinates(x);
      bodies[i].GetVelocity(v);
      x[0] = x[0] + 0.5*dt*v[0];
      x[1] = x[1] + 0.5*dt*v[1];
      x[2] = x[2] + 0.5*dt*v[2];
      bodies[i].SetCoordinates(x);
      bodies[i].SetVelocity(v);
    }

    cout << t << " Passo 2" << endl;
    octree = Generate_Octree(N, center, L); 
    double g[N][3];
    for(int i=0;i<N;i++){  //parallelizzare
      g[i][0] = 0; //inizializzo gravità
      g[i][1] = 0;
      g[i][2] = 0;
      Compute_Gravity(i, octree, g[i], L);
    }
    Delete_Octree(octree); 

    cout << t << " Passo 3" << endl;
    for(int i=0;i<N;i++){  //parallelizzare
      double x[3];
      double v[3];
      bodies[i].GetCoordinates(x);
      bodies[i].GetVelocity(v);
      v[0] = v[0] + dt*g[i][0];
      v[1] = v[1] + dt*g[i][1];
      v[2] = v[2] + dt*g[i][2];
      x[0] = x[0] + 0.5*dt*v[0];
      x[1] = x[1] + 0.5*dt*v[1];
      x[2] = x[2] + 0.5*dt*v[2];
      bodies[i].SetCoordinates(x);
      bodies[i].SetVelocity(v);
    }
    t += dt;
  }

//--------------------------------------------------------------------------+
//|                           ESPORTO I RISULTATI                           |
//+-------------------------------------------------------------------------+ 

  save(output_file,bodies,N);
}