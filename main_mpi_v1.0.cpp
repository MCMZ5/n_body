//compilare con mpic++

#include "mpi.h"
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <math.h>

using namespace std;

//--------------------------------------------------------------------------+
//|                   INTERFACCIA: VARIABILI DA MODIFICARE                  |
//+-------------------------------------------------------------------------+ 
#define L          100   //Dimensione dello spazio dell'octree
#define T          2  //Tempo totale di simulazione

//--------------------------------------------------------------------------+
//|                              STRUCT BODY                                |
//+-------------------------------------------------------------------------+ 
struct Body{
	double mass; 
  double coord[3]; 
  double vel[3]; 

  Body(){};
  Body(double m, double v[3], double x[3]){
    mass = m;
    for(int i=0;i<3;i++){
      coord[i] = x[i];
      vel[i] = v[i];
    }
  };
  Body(double m, double x[3]){
    mass = m;
    for(int i=0;i<3;i++){
      coord[i] = x[i];
      vel[i] = 0;
    }
  };
  bool IsInCube(double c[], double l){
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
  };
  bool IsInCubes(double c[], double l, int n){
    int a = 0;
    for(int i=0;i<n;i++){
      if(IsInCube(&c[3*i], l)==true) a++;
    }
    if(a==0){
      return false;
    }
    else{
      return true;
    }
  };
  void Print(){
    cout << "m = " << mass << "\n";
    cout << "x = " << coord[0] << "," << coord[1] << "," << coord[2] << "\n";
    cout << "v = " << vel[0] << "," << vel[1] << "," << vel[2] << "\n";
  }
};

//--------------------------------------------------------------------------+
//|                             STRUCT VECT3D                               |
//+-------------------------------------------------------------------------+ 
struct Vect3D{
  double v[3];
};

//--------------------------------------------------------------------------+
//|                             CLASSE BARNES                               |
//+-------------------------------------------------------------------------+ 
class Barnes{
    public:
        Barnes(){
          empty = false;
        }
        Barnes(double m, double *x){
          cm.mass=m;
          cm.coord[0]=x[0];
          cm.coord[1]=x[1];
          cm.coord[2]=x[2];
          cm.vel[0]=0.;
          cm.vel[1]=0.;
          cm.vel[2]=0.;
          empty = false;
        }
        Barnes(Body &b){
          cm = b;
          empty = false;
        }
        Barnes(int &N, double c[], double Len, Body bodies[]){
          //creo il nodo base
          double m = 0;
          double mm;
          double p[3] = {0,0,0};  //per il centro di massa
          int nb = 0;
          int lf;
          for(int i=0;i<N;i++){
            if(bodies[i].IsInCube(c,Len)==true){
              nb++;
              if(nb==1){
                lf=i;
              }
              mm = bodies[i].mass;
              m += mm;
              p[0] += bodies[i].coord[0]*mm;
              p[1] += bodies[i].coord[1]*mm;
              p[2] += bodies[i].coord[2]*mm;
            }
          }
          if(nb==0){
            empty = true;
          }
          if(nb>0){  //lo considero solo se contiene corpi
            p[0]=p[0]/m;
            p[1]=p[1]/m;
            p[2]=p[2]/m;
            if (nb==1){
              cm = bodies[lf];
              empty = false;
            }
            else{
              cm.mass = m;
              cm.coord[0]=p[0];
              cm.coord[1]=p[1];
              cm.coord[2]=p[2];
              empty = false;
            }
          }
          if(nb>1){
            //divido in ottanti e chiamo il costruttore
            double l = Len/2;
            //salvo il centro madre
            double cc[3] = {c[0],c[1],c[2]};
            //primo ottante
            c[0] = cc[0] - l/2;
            c[1] = cc[1] - l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd1 = new Barnes(N,c,l,bodies);
            if(nd1->empty==true){
              delete nd1;
            }
            else{
              subnodes.push_back(nd1);
            }
            //secondo ottante
            c[2] = cc[2] + l/2;
            Barnes* nd2 = new Barnes(N,c,l,bodies);
            if(nd2->empty==true){
              delete nd2;
            }
            else{
              subnodes.push_back(nd2);
            }
            //terzo ottante
            c[1] = cc[1] + l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd3 = new Barnes(N,c,l,bodies);
            if(nd3->empty==true){
              delete nd3;
            }
            else{
              subnodes.push_back(nd3);
            }
            //quarto ottante
            c[2] = cc[2] + l/2;
            Barnes* nd4 = new Barnes(N,c,l,bodies);
            if(nd4->empty==true){
              delete nd4;
            }
            else{
              subnodes.push_back(nd4);
            }
            //quinto ottante
            c[0] = cc[0] + l/2;
            c[1] = cc[1] - l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd5 = new Barnes(N,c,l,bodies);
            if(nd5->empty==true){
              delete nd5;
            }
            else{
              subnodes.push_back(nd5);
            }
            //sesto ottante
            c[2] = cc[2] + l/2;
            Barnes* nd6 = new Barnes(N,c,l,bodies);
            if(nd6->empty==true){
              delete nd6;
            }
            else{
              subnodes.push_back(nd6);
            }
            //settimo ottante
            c[1] = cc[1] + l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd7 = new Barnes(N,c,l,bodies);
            if(nd7->empty==true){
              delete nd7;
            }
            else{
              subnodes.push_back(nd7);
            }
            //ottavo ottante
            c[2] = cc[2] + l/2;
            Barnes* nd8 = new Barnes(N,c,l,bodies);
            if(nd8->empty==true){
              delete nd8;
            }
            else{
              subnodes.push_back(nd8);
            }
          }
        }
        ~Barnes(){

        }
        void Delete_Subnodes(){
          if(subnodes.size()>0){
            for(int i=0;i<subnodes.size();i++){
              subnodes[i]->Delete_Subnodes();
              delete subnodes[i];
            }
          }
        }

        Body cm;
        vector<Barnes*> subnodes; //sottonodi dell'octree
        bool empty;
};

//--------------------------------------------------------------------------+
//|                        DICHIARAZIONI FUNZIONI                           |
//+-------------------------------------------------------------------------+ 
//funzione per sapere quanti corpi sono descritti dal file
int count(string file);
//funzione per leggere il file con un array di corpi
void load(string file, Body bodies[], int dimension);
//funzione per salvare su file i valori di un array di corpi
void save(string file, Body bodies[], int dimension);
//funzione per costruire la parte di octree locale di ogni processo
Barnes* Generate_Octree(int &N, double c[], double Len, Body bodies[], int myrank);
//funzione per gestire l'allocazione dinamica dell'octree locale
void Delete_Octree(Barnes* root);
//funzioni Distance per il calcolo di distanze tra Body
double Distance(Body b1, Body b2);
void Vectorial_Distance(Body b1, Body b2, double *d);
//funzione per l'equazione di Newton
void acc(double m, double d, double *x, double *a);
//funzioni per il calcolo della gravità
void Compute_Gravity(Body body, Barnes* current, double *g, double Len);
void Compute_Gravity_cm(Body body, int i, Body nodes_cm[], int Ncm, int Np, int rank, int rank0, vector<int> v[], double *g, double Len);

//--------------------------------------------------------------------------+
//|                                 MAIN                                    |
//+-------------------------------------------------------------------------+ 
int main(int argc, char **argv){

//Inizializzo MPI
  MPI_Init(&argc, &argv);
  double t0,t1,t2,t3,t4,t5;
  t0 = MPI_Wtime();
  int myrank[2];
  int Np[2];
  MPI_Comm COMM_LVL1;
  MPI_Comm COMM_LVL2;
  MPI_Comm_rank(MPI_COMM_WORLD, myrank);
  MPI_Comm_size(MPI_COMM_WORLD, Np);
  int color;
  int color1;
  if(myrank[0]<8){
    color = myrank[0]%8;
    color1 = 1;
  }
  if(myrank[0]>=8){
    color = (myrank[0]-8)/((Np[0]-8)/8);
    color1 = 0;
  }
  MPI_Comm_split(MPI_COMM_WORLD, color1, myrank[0], &COMM_LVL1);
  MPI_Comm_split(MPI_COMM_WORLD, color, myrank[0], &COMM_LVL2);
  MPI_Comm_rank(COMM_LVL2, &myrank[1]);
  MPI_Comm_size(COMM_LVL2, &Np[1]);
  

//Dichiarazione variabili necessarie
  double   t     = 0.;
  double   dt    = 0.0001; //salto temporale
  string input_file = "nbody_start_16000.txt";
  string output_file = "nbody_end.txt";
  int N;
  MPI_Status status;
  double cc[3] = {0,0,0}; //centro del sistema
  double ccc[3] = {0,0,0};
  double l;
  int Npc;
  if(Np[0]<=8){
    Npc=Np[0];
  }
  else{
    Npc=Np[1];
  }

//definizione tipo per mandare le struct Body e vect3D
  MPI_Datatype Bodytype;
  MPI_Type_contiguous(7,MPI_DOUBLE,&Bodytype);
  MPI_Type_commit(&Bodytype);
  MPI_Aint extent,lb;
  MPI_Type_get_extent(Bodytype,&lb, &extent);

  MPI_Datatype type3D;
  MPI_Type_contiguous(3,MPI_DOUBLE,&type3D);
  MPI_Type_commit(&type3D);
  MPI_Aint extent1,lb1;
  MPI_Type_get_extent(type3D,&lb1, &extent1);

//--------------------------------------------------------------------------+
//|                     IMPORTO I CORPI DEL SISTEMA                         |
//+-------------------------------------------------------------------------+
  if (myrank[0]==0) {
    N = count(input_file);
    cout << "Dati contati: " << N << "\n";
  }
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD); //serve per allocare vettori abbastanza grandi
  
  Body * bodies = NULL;
  bodies = new Body[N];
  Body * bodies_end = NULL;
  bodies_end = new Body[N];
  Body * octree = NULL;
  octree = new Body[8*N];
  Barnes* local_octree[8/Npc]; 
  Body cm1[8/Npc];
  double g[N][3];
  Body *cm = NULL;
  int Ncm;
  if(Np[0]<=8){
    Ncm=8;
  }
  if(Np[0]>8){
    Ncm=64;
  }
  cm = new Body[Ncm];
  double *ctr = NULL;
  ctr = new double[3*(8/Npc)];
  Body send[N];
  Body newbodies[N];
  int n_recv;
  int n[8];
  int n_end[Np[0]];
  int shift[Np[0]];
  int a;
  if (myrank[0]==0) {
    load(input_file,bodies,N);
  }
  t1 = MPI_Wtime();

//--------------------------------------------------------------------------+
//|                       DIVISIONE DEL PROBLEMA                            |
//+-------------------------------------------------------------------------+ 

int lvl=0;
if(Np[0]>8){
  lvl=1;
}
l=L;
for(int k=0;k<=lvl;k++){
  int Npp;
  if(Np[k]<=8){
    Npp=Np[k];
  }
  else{
    Npp=8;
  }
  if (myrank[k]==0){
    //divido in ottanti 
    l = l/2;
    //posizioni dei centri
    if(k>0){
      ccc[0]=ctr[0];
      ccc[1]=ctr[1];
      ccc[2]=ctr[2];
    }
    double c[8][3] = {{ccc[0] - l/2,ccc[1] - l/2,ccc[2] - l/2},
                      {ccc[0] - l/2,ccc[1] - l/2,ccc[2] + l/2},
                      {ccc[0] - l/2,ccc[1] + l/2,ccc[2] - l/2},
                      {ccc[0] - l/2,ccc[1] + l/2,ccc[2] + l/2},
                      {ccc[0] + l/2,ccc[1] - l/2,ccc[2] - l/2},
                      {ccc[0] + l/2,ccc[1] - l/2,ccc[2] + l/2},
                      {ccc[0] + l/2,ccc[1] + l/2,ccc[2] - l/2},
                      {ccc[0] + l/2,ccc[1] + l/2,ccc[2] + l/2}};
    //spedisco ai processi fratelli i centri di loro competenza
    if(k==0){
      MPI_Scatter(&c, 3*(8/Npp), MPI_DOUBLE, ctr, 3*(8/Npp), MPI_DOUBLE, 0, COMM_LVL1);
    }
    if(k==1){
      MPI_Scatter(&c, 3*(8/Npp), MPI_DOUBLE, ctr, 3*(8/Npp), MPI_DOUBLE, 0, COMM_LVL2);
    }
    vector<int> o[8];
    for(int j=0;j<8;j++){
      n[j]=0;
    }
    for(int i=0;i<N;i++){
      for(int j=0;j<8;j++){
        if(bodies[i].IsInCube(c[j],l)==true){
          o[j].push_back(i);
          n[j]++;
        }
      }
    }
    int aa = 0;
    for(int j=0;j<8;j++){
      for(int y=0;y<o[j].size();y++){
        octree[aa+y]=bodies[o[j][y]];
      }
      aa+=o[j].size();
    }
  }
  //ricevo nei processi fratelli i centri di loro competenza
  if(myrank[k]>0 && myrank[k]<8){
    l = L/2;
    double c; //serve solo per il compilatore
    if(k==0){
      MPI_Scatter(&c, 3*(8/Npp), MPI_DOUBLE, ctr, 3*(8/Npp), MPI_DOUBLE, 0, COMM_LVL1);
    }
    if(k==1){
      MPI_Scatter(&c, 3*(8/Npp), MPI_DOUBLE, ctr, 3*(8/Npp), MPI_DOUBLE, 0, COMM_LVL2);
    }
  }
  //condivido la visione globale
  if(k==0 && myrank[0]<8){
    MPI_Bcast(&n, 8, MPI_INT, 0, COMM_LVL1); 
  }
  if(k==1){
    MPI_Bcast(&n, 8, MPI_INT, 0, COMM_LVL2); 
  }
  for(int j=0;j<Npp;j++){
    n[j]=n[(j*(8/Npp))];
    for(int i=1;i<(8/Npp);i++){
      n[j]+=n[(j*(8/Npp))+i];
    }
  }
  N=n[myrank[k]];
  shift[0]=0;
  for(int i=1;i<Npp;i++){
    shift[i]=shift[i-1]+n[i-1];
  } 
  if(k==0 && myrank[0]<8){
    MPI_Scatterv(octree,n,shift,Bodytype,bodies,N,Bodytype,0,COMM_LVL1);
  }
  if(k==1){
    MPI_Scatterv(octree,n,shift,Bodytype,bodies,N,Bodytype,0,COMM_LVL2);
  }
}
t2 = MPI_Wtime();

//--------------------------------------------------------------------------+
//|                           OCTREE E GRAVITA'                             |
//+-------------------------------------------------------------------------+ 

  for(int i=0;i<(8/Npc);i++){
    local_octree[i] = Generate_Octree(N, &ctr[3*i], l, bodies, myrank[0]); 
    cm1[i] = local_octree[i]->cm;
  }
  MPI_Allgather(&cm1, 8/Npc, Bodytype, cm, 8/Npc, Bodytype, MPI_COMM_WORLD);

  vector<int> to[Np[0]];
  int rrank;
  if(Np[0]<=8){
    rrank=myrank[0];
  }
  if(Np[0]>8){
    rrank=myrank[1];
  }

  for(int i=0;i<N;i++){
    //inizializzo gravità
    g[i][0] = 0;
    g[i][1] = 0;
    g[i][2] = 0;
    for(int j=0;j<(8/Npc);j++){
      Compute_Gravity(bodies[i], local_octree[j], g[i], l);
    }
    Compute_Gravity_cm(bodies[i], i, cm, Ncm,  Np[1], rrank, myrank[0], to, g[i], l);
  }

  //Computazione gravità negli altri processi
  int n_rec[Np[0]];
  //invio dimensioni
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    MPI_Gather(&sn,1,MPI_INT,&n_rec,1,MPI_INT,i,MPI_COMM_WORLD);
  }
  //operazioni parallele
  int tot=0;
  for(int j=0;j<Np[0];j++){
    tot+=n_rec[j];
    if(j==0){
      shift[j]=0;
    }
    else{
      shift[j]=shift[j-1]+n_rec[j-1];
    }
  }
  Body in[tot];
  Vect3D gy[tot];
  //invio corpi
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    Body out[sn];  //cosa invio
    for(int j=0;j<sn;j++){
      out[j]=bodies[to[i][j]];
    }
    MPI_Gatherv(&out,sn,Bodytype,&in,n_rec,shift,Bodytype,i,MPI_COMM_WORLD);
  }
  //Computo gravità
  for(int j=0;j<tot;j++){
    //inizializzo gravità
    gy[j].v[0] = 0;
    gy[j].v[1] = 0;
    gy[j].v[2] = 0;
    for(int k=0;k<(8/Npc);k++){
      Compute_Gravity(in[j], local_octree[k], gy[j].v, l);
    }
  }
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    Vect3D gg[sn];
    MPI_Scatterv(&gy,n_rec,shift,type3D,&gg,sn,type3D,i,MPI_COMM_WORLD);
    for(int j=0;j<sn;j++){
      g[to[i][j]][0]+=gg[j].v[0];
      g[to[i][j]][1]+=gg[j].v[1];
      g[to[i][j]][2]+=gg[j].v[2];
    }
  }

  for(int i=0;i<(8/Npc);i++){
    Delete_Octree(local_octree[i]); 
  } 
t3 = MPI_Wtime();
//--------------------------------------------------------------------------+
//|                            CICLO SUL TEMPO                              |
//+-------------------------------------------------------------------------+

  while(t<=T){
    if(myrank[0]==0){
      cout << "ciclo al tempo " << t << endl;
    }

//--------------------------------------------------------------------------+
//|                      VELOCITY VERLET PRIMA PARTE                        |
//+-------------------------------------------------------------------------+ 
    for(int i=0;i<N;i++){
      bodies[i].vel[0] = bodies[i].vel[0] + 0.5*dt*g[i][0];
      bodies[i].vel[1] = bodies[i].vel[1] + 0.5*dt*g[i][1];
      bodies[i].vel[2] = bodies[i].vel[2] + 0.5*dt*g[i][2];
      bodies[i].coord[0] = bodies[i].coord[0] + dt*bodies[i].vel[0];
      bodies[i].coord[1] = bodies[i].coord[1] + dt*bodies[i].vel[1];
      bodies[i].coord[2] = bodies[i].coord[2] + dt*bodies[i].vel[2];
    } 

//--------------------------------------------------------------------------+
//|                              SCAMBIO CORPI                              |
//+-------------------------------------------------------------------------+ 
    int n_send;
    vector<int> s;
    vector<int> h; 
    //inizializzo newbodies e n_recv sullo stato attuale
    for(int i=0;i<N;i++){
      newbodies[i]=bodies[i];
    }
    n_recv = N;
    N = 0; //lo ricalcolo nello scambio 
    for(int j=0;j<Np[0];j++){
      for(int i=0;i<n_recv;i++){
        if(newbodies[i].IsInCube(cc,L)==false){
           //bodies fuori dal sistema, non li considero oltre
        }
        else if(newbodies[i].IsInCubes(ctr, l, 8/Npc)==false){
          s.push_back(i); //salvo l'indice dei bodies fuori dal sottosistema da inviare
        }
        else{
          h.push_back(i); //salvo l'indice dei bodies da tenere
        }
      }
      for(int i=0;i<s.size();i++){
        send[i]= newbodies[s[i]];
      }
      for(int i=0;i<h.size();i++){
        bodies[N]= newbodies[h[i]];   
        N++;                          
      }
      if(j<(Np[0]-1)){
        //spedizioni di scambio a cerchio
        n_send=s.size();
        int id_s = myrank[0]+1;
        int id_r = myrank[0]-1;
        if(myrank[0]==(Np[0]-1)){
          id_s = 0;
        }
        if(myrank[0]==0){
          id_r = Np[0]-1;
        }
        if(myrank[0]%2==0){
          MPI_Send(&n_send, 1, MPI_INT, id_s, j+1, MPI_COMM_WORLD);
          MPI_Recv(&n_recv, 1, MPI_INT, id_r, j+1, MPI_COMM_WORLD, &status);

          MPI_Send(&send, n_send, Bodytype, id_s, j+8, MPI_COMM_WORLD);
          MPI_Recv(&newbodies, n_recv, Bodytype, id_r, j+8, MPI_COMM_WORLD, &status);
        }
        else{
          MPI_Recv(&n_recv, 1, MPI_INT, id_r, j+1, MPI_COMM_WORLD, &status);
          MPI_Send(&n_send, 1, MPI_INT, id_s, j+1, MPI_COMM_WORLD);

          MPI_Recv(&newbodies, n_recv, Bodytype, id_r, j+8, MPI_COMM_WORLD, &status);
          MPI_Send(&send, n_send, Bodytype, id_s, j+8, MPI_COMM_WORLD);
        }
      }
      s.clear();
      h.clear();
    }

//--------------------------------------------------------------------------+
//|                           OCTREE E GRAVITA'                             |
//+-------------------------------------------------------------------------+ 

  //aggiorno l'octree e la gravità
  for(int i=0;i<(8/Npc);i++){
    local_octree[i] = Generate_Octree(N, &ctr[3*i], l, bodies, myrank[0]); 
    cm1[i] = local_octree[i]->cm;
  }
  MPI_Allgather(&cm1, 8/Npc, Bodytype, cm, 8/Npc, Bodytype, MPI_COMM_WORLD);

  vector<int> to[Np[0]];
  int rrank;
  if(Np[0]<=8){
    rrank=myrank[0];
  }
  if(Np[0]>8){
    rrank=myrank[1];
  }

  for(int i=0;i<N;i++){
    //inizializzo gravità
    g[i][0] = 0;
    g[i][1] = 0;
    g[i][2] = 0;
    for(int j=0;j<(8/Npc);j++){
      Compute_Gravity(bodies[i], local_octree[j], g[i], l);
    }
    Compute_Gravity_cm(bodies[i], i, cm, Ncm,  Np[1], rrank, myrank[0], to, g[i], l);
  }

  //Computazione gravità negli altri processi
  int n_rec[Np[0]];
  //invio dimensioni
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    MPI_Gather(&sn,1,MPI_INT,&n_rec,1,MPI_INT,i,MPI_COMM_WORLD);
  }
  //operazioni parallele
  int tot=0;
  for(int j=0;j<Np[0];j++){
    tot+=n_rec[j];
    if(j==0){
      shift[j]=0;
    }
    else{
      shift[j]=shift[j-1]+n_rec[j-1];
    }
  }
  Body in[tot];
  Vect3D gy[tot];
  //invio corpi
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    Body out[sn];  //cosa invio
    for(int j=0;j<sn;j++){
      out[j]=bodies[to[i][j]];
    }
    MPI_Gatherv(&out,sn,Bodytype,&in,n_rec,shift,Bodytype,i,MPI_COMM_WORLD);
  }
  //Computo gravità
  for(int j=0;j<tot;j++){
    //inizializzo gravità
    gy[j].v[0] = 0;
    gy[j].v[1] = 0;
    gy[j].v[2] = 0;
    for(int k=0;k<(8/Npc);k++){
      Compute_Gravity(in[j], local_octree[k], gy[j].v, l);
    }
  }
  //restituisco valori della gravità
  for(int i=0;i<Np[0];i++){
    int sn = to[i].size(); //Numero di invii
    Vect3D gg[sn];
    MPI_Scatterv(&gy,n_rec,shift,type3D,&gg,sn,type3D,i,MPI_COMM_WORLD);
    for(int j=0;j<sn;j++){
      g[to[i][j]][0]+=gg[j].v[0];
      g[to[i][j]][1]+=gg[j].v[1];
      g[to[i][j]][2]+=gg[j].v[2];
    }
  }

  for(int i=0;i<(8/Npc);i++){
    Delete_Octree(local_octree[i]); 
  } 

//--------------------------------------------------------------------------+
//|                      VELOCITY VERLET SECONDA PARTE                      |
//+-------------------------------------------------------------------------+
    //velocity Verlet parte 2
    for(int i=0;i<N;i++){
      bodies[i].vel[0] = bodies[i].vel[0] + 0.5*dt*g[i][0];
      bodies[i].vel[1] = bodies[i].vel[1] + 0.5*dt*g[i][1];
      bodies[i].vel[2] = bodies[i].vel[2] + 0.5*dt*g[i][2]; 
    }
    t += dt;
  }
  t4 = MPI_Wtime();

//--------------------------------------------------------------------------+
//|                        ESPORTO I RISULTATI                              |
//+-------------------------------------------------------------------------+
  MPI_Gather(&N,1,MPI_INT,n_end,1,MPI_INT,0,MPI_COMM_WORLD); //raccolgo n dai vari processi
  shift[0]=0;
  for(int i=1;i<Np[0];i++){
    shift[i]=shift[i-1]+n_end[i-1];
  }
  MPI_Gatherv(bodies,N,Bodytype,bodies_end,n_end,shift,Bodytype,0,MPI_COMM_WORLD);
  if (myrank[0]==0) {
    N=0;
    for(int i=0;i<Np[0];i++){
        N += n_end[i];
    }
    save(output_file,bodies_end,N);
    cout << "Dati salvati: " << N << "\n";
  } 
  delete octree;
  delete cm;
  delete ctr;
  delete bodies_end;
  delete bodies;
  t5 = MPI_Wtime();
  MPI_Finalize();
  if(myrank[0]==0){
    ofstream time_data;
	  time_data.open("Performance_"+to_string(Np[0])+".txt");
    time_data << "Number of bodies: " << N << "\n";    
    time_data << "Total time (MPI) is " << t5-t0 << "\n";
    time_data << "Import time (MPI) is " << t1-t0 << "\n";
    time_data << "Division time (MPI) is " << t2-t1 << "\n";
    time_data << "First Octree time (MPI) is " << t3-t2 << "\n";
    time_data << "Export time (MPI) is " << t5-t4 << "\n";
  }
  return 0;
}

//--------------------------------------------------------------------------+
//|                              FUNZIONI                                   |
//+-------------------------------------------------------------------------+
int count(string file){
	int c;
	c = 0;
	ifstream txt;
	txt.open(file);
	string s;
	while(txt.eof() == false){
		getline(txt,s);
		if(s != "" && s.at(0) != '#'){
			c++;
		}
	}
	txt.close();
return c;
}

void load(string file, Body bodies[], int dimension){
	ifstream txt;
	txt.open(file);
	string s;
	string p[3];			
	string v[3];	
	string m;		
	int i;
	i = 0;  //numero riga
    int k;
	k = 0;	//numero colonna
	while(txt.eof() == false){
		getline(txt,s);
		if(s != "" && s.at(0) != '#'){
			if(i < dimension){
				for(int j=0; j < s.length(); j++){
					if(s[j] != ' ' && k==0){
						m=m+s[j];
					}
					else if(s[j] != ' ' && k>0 && k<4){
						p[k-1]=p[k-1]+s[j];
					}
          			else if(s[j] != ' ' && k>3 && k<7){
						v[k-4]=v[k-4]+s[j];
					}
					else{
            			k++;
					}
				}
				double x[3] = {stod(p[0]),stod(p[1]),stod(p[2])};
				double vl[3] = {stod(v[0]),stod(v[1]),stod(v[2])};
				Body bd(stof(m),vl,x);
				bodies[i]=bd;
				m = "";
       			p[0] = "";
				p[1] = "";
        		p[2] = "";
				v[0] = "";
				v[1] = "";
        		v[2] = "";
				i++;
			}
			else{
				cout << "Errore: file oltre il limite del vettore" << endl;
				cout << "Non tutti i dati sono stati raccolti" << endl;
				break;		
			}
			k = 0;
		}
	}
	txt.close();
}

void save(string file, Body bodies[], int dimension){
	ofstream txt;
	txt.open(file);
	for(int i=0;i<dimension;i++){
		txt << setiosflags(ios::scientific);
		txt << bodies[i].mass << " ";
		txt << bodies[i].coord[0] << " ";
		txt << bodies[i].coord[1] << " ";
		txt << bodies[i].coord[2] << " ";
		txt << bodies[i].vel[0] << " ";
		txt << bodies[i].vel[1] << " ";
		txt << bodies[i].vel[2] << endl;
	}
	txt.close();
}

Barnes* Generate_Octree(int &N, double c[], double Len, Body bodies[], int myrank){
    //creo il nodo root
    double m = 0;
    double mm;
    double p[3] = {0,0,0};  //per il centro di massa
    for(int i=0;i<N;i++){
        mm = bodies[i].mass;
        m += mm;
        p[0] += bodies[i].coord[0]*mm;
        p[1] += bodies[i].coord[1]*mm;
        p[2] += bodies[i].coord[2]*mm;
    }
    p[0]=p[0]/m;
    p[1]=p[1]/m;
    p[2]=p[2]/m;
    Barnes* root = new Barnes(m,p);
    if(N>1){
      //divido in ottanti e chiamo il costruttore
      double l = Len/2;
      //salvo il centro madre
      double cc[3] = {c[0],c[1],c[2]};
      //primo ottante
      c[0] = cc[0] - l/2;
      c[1] = cc[1] - l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd1 = new Barnes(N,c,l,bodies);
      if(nd1->empty==true){
          delete nd1;
      }
      else{
          root->subnodes.push_back(nd1);
      }
      //secondo ottante
      c[2] = cc[2] + l/2;
      Barnes* nd2 = new Barnes(N,c,l,bodies);
      if(nd2->empty==true){
          delete nd2;
      }
      else{
          root->subnodes.push_back(nd2);
      }
      //terzo ottante
      c[1] = cc[1] + l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd3 = new Barnes(N,c,l,bodies);
      if(nd3->empty==true){
          delete nd3;
      }
      else{
          root->subnodes.push_back(nd3);
      }
      //quarto ottante
      c[2] = cc[2] + l/2;
      Barnes* nd4 = new Barnes(N,c,l,bodies);
      if(nd4->empty==true){
          delete nd4;
      }
      else{
          root->subnodes.push_back(nd4);
      }
      //quinto ottante
      c[0] = cc[0] + l/2;
      c[1] = cc[1] - l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd5 = new Barnes(N,c,l,bodies);
      if(nd5->empty==true){
          delete nd5;
      }
      else{
          root->subnodes.push_back(nd5);
      }
      //sesto ottante
      c[2] = cc[2] + l/2;
      Barnes* nd6 = new Barnes(N,c,l,bodies);
      if(nd6->empty==true){
          delete nd6;
      }
      else{
          root->subnodes.push_back(nd6);
      }
      //settimo ottante
      c[1] = cc[1] + l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd7 = new Barnes(N,c,l,bodies);
      if(nd7->empty==true){
          delete nd7;
      }
      else{
          root->subnodes.push_back(nd7);
      }
      //ottavo ottante
      c[2] = cc[2] + l/2;
      Barnes* nd8 = new Barnes(N,c,l,bodies);
      if(nd8->empty==true){
          delete nd8;
      }
      else{
          root->subnodes.push_back(nd8);
      }
    //risistemo la variabile centro passata by reference
      c[0] = cc[0];
      c[1] = cc[1];
      c[2] = cc[2];
    }
	return root;
}

void Delete_Octree(Barnes* root){ //gestione allocazione dinamica
	root->Delete_Subnodes();
	delete root;
}

double Distance(Body b1, Body b2){

	double a[3] = {b1.coord[0],b1.coord[1],b1.coord[2]};
	double b[3] = {b2.coord[0],b2.coord[1],b2.coord[2]};

    double d = sqrt(abs(((a[0]-b[0])*(a[0]-b[0]))+
                    ((a[1]-b[1])*(a[1]-b[1]))+
                    ((a[2]-b[2])*(a[2]-b[2]))));

    return d;
}

void Vectorial_Distance(Body b1, Body b2, double *d){

	double a[3] = {b1.coord[0],b1.coord[1],b1.coord[2]};
	double b[3] = {b2.coord[0],b2.coord[1],b2.coord[2]};

  d[0] = b[0]-a[0];
	d[1] = b[1]-a[1];
	d[2] = b[2]-a[2]; 

}

void acc(double m, double d, double *x, double *a){
  double mr = m/(d*d*d);
  a[0] += x[0]*mr;
  a[1] += x[1]*mr;
  a[2] += x[2]*mr;
}

void Compute_Gravity(Body body, Barnes* current, double *g, double Len){
    double dx[3]; //per la distanza vettoriale tra corpi
    double m;
    double l;
	  double d;
    for(int k=0;k<current->subnodes.size();k++){
        current=current->subnodes[k];
        d = Distance(body,current->cm);
        l = Len/(2);
        if(d<l && current->subnodes.size()>0){
			    Compute_Gravity(body,current,g,l);
        }
        else if(d<=pow(10,-16)){ //escludo i casi di distanze tra un corpo e se stesso con precisione double
        
        }
        else{
            m = current->cm.mass;
            Vectorial_Distance(body,current->cm,dx);
            acc(m,d,dx,g);  //aggiungo gravità 
          }
      }
}

void Compute_Gravity_cm(Body body, int i, Body nodes_cm[], int Ncm, int Np, int rank, int rank0, vector<int> v[], double *g, double Len){
  double dx[3]; //per la distanza vettoriale tra corpi
  double m;
  double l;
	double d;
  for(int k=0;k<Ncm;k++){ 
        if(k>=((rank*(8/Np))+(rank0/8)) && k<(((rank+1)*(8/Np)))+(rank0/8)){
          //escludo i centri di massa dell'octree locale
        }
        else{
          m = nodes_cm[k].mass;
          if(m>0){ //escludo i nodi vuoti
            d = Distance(body,nodes_cm[k]);
            l = Len/(2);
            if(d<l){
              v[(k/(8/Np))].push_back(i);
            }
            else if(d<=pow(10,-16)){ //escludo i casi di distanze tra un corpo e se stesso con precisione double

            }
            else{
              Vectorial_Distance(body,nodes_cm[k],dx);
              acc(m,d,dx,g);
            }
          }
        }
  }
}