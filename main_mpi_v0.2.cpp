//compilare con mpiCC

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

//INTERFACCIA: VARIABILI DA MODIFICARE
#define L          40   //Dimensione dello spazio dell'octree
#define T          0.04  //Tempo totale di simulazione

//STRUCT BODY
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
  void Print(){
    cout << "m = " << mass << "\n";
    cout << "x = " << coord[0] << "," << coord[1] << "," << coord[2] << "\n";
    cout << "v = " << vel[0] << "," << vel[1] << "," << vel[2] << "\n";
  }
};

//CLASSE BARNES
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
            //if(m == 16.632){
            //  cout << m << "," << p[0] << "," << p[1] << "," << p[2] << "\n";
            //}
            if (nb==1){
              cm = bodies[lf];
              empty = false;
            }
            else{
              cm.mass = m;
              cm.coord[0]=p[0];
              cm.coord[0]=p[1];
              cm.coord[0]=p[2];
              empty = false;
            }
            //cout << "Generating node of level " << lvl << endl;
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
            Barnes* nd1 = new Barnes(N,c,l,bodies);//static Barnes nd1(N,c,l,lvl+1);
            if(nd1->empty==true){
              delete nd1;
            }
            else{
              subnodes.push_back(nd1);
            }
            //secondo ottante
            c[2] = cc[2] + l/2;
            Barnes* nd2 = new Barnes(N,c,l,bodies);//static Barnes nd2(N,c,l,lvl+1);
            if(nd2->empty==true){
              delete nd2;
            }
            else{
              subnodes.push_back(nd2);
            }
            //terzo ottante
            c[1] = cc[1] + l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd3 = new Barnes(N,c,l,bodies);//static Barnes nd3(N,c,l,lvl+1);
            if(nd3->empty==true){
              delete nd3;
            }
            else{
              subnodes.push_back(nd3);
            }
            //quarto ottante
            c[2] = cc[2] + l/2;
            Barnes* nd4 = new Barnes(N,c,l,bodies);//static Barnes nd4(N,c,l,lvl+1);
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
            Barnes* nd5 = new Barnes(N,c,l,bodies);//static Barnes nd5(N,c,l,lvl+1);
            if(nd5->empty==true){
              delete nd5;
            }
            else{
              subnodes.push_back(nd5);
            }
            //sesto ottante
            c[2] = cc[2] + l/2;
            Barnes* nd6 = new Barnes(N,c,l,bodies);//static Barnes nd6(N,c,l,lvl+1);
            if(nd6->empty==true){
              delete nd6;
            }
            else{
              subnodes.push_back(nd6);
            }
            //settimo ottante
            c[1] = cc[1] + l/2;
            c[2] = cc[2] - l/2;
            Barnes* nd7 = new Barnes(N,c,l,bodies);//static Barnes nd7(N,c,l,lvl+1);
            if(nd7->empty==true){
              delete nd7;
            }
            else{
              subnodes.push_back(nd7);
            }
            //ottavo ottante
            c[2] = cc[2] + l/2;
            Barnes* nd8 = new Barnes(N,c,l,bodies);//static Barnes nd8(N,c,l,lvl+1);
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

//DICHIARAZIONI FUNZIONI
//funzione per sapere quanti corpi sono descritti dal file
int count(string file);
//funzione per leggere il file con un array di corpi
void load(string file, Body bodies[], int dimension);
//funzione per salvare su file i valori di un array di corpi
void save(string file, Body bodies[], int dimension);
//funzione per costruire la parte di octree locale di ogni processo
Barnes* Generate_Octree(int &N, double c[], double Len, Body bodies[]);
//funzione per gestire l'allocazione dinamica dell'octree locale
void Delete_Octree(Barnes* root);
//funzioni Distance per il calcolo di distanze tra Body
double Distance(Body b1, Body b2);
void Vectorial_Distance(Body b1, Body b2, double *d);
//funzione per l'equazione di Newton
void acc(double m, double d, double *x, double *a);
//funzione per il calcolo della gravità
void Compute_Gravity(Body body, Barnes* current, double *g, double Len);
void Compute_Gravity(Body body, Barnes* current, double *g, double Len, Body nodes_cm[], int N_nodes, int self);

//MAIN
int main(int argc, char **argv){

//Inizializzo MPI
  MPI_Init(&argc, &argv);
  double t0,t1;
  t0 = MPI_Wtime();
  int myrank;
  int Np;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &Np);

//Dichiarazione variabili necessarie
  double   t     = 0.;
  double   dt    = 0.05; //salto temporale
  string input_file = "nbody_start.txt";
  string output_file = "nbody_end.txt";
  int N;
  MPI_Status status;
  Barnes* local_octree;
  double cc[3] = {0,0,0}; //centro del sistema
  double l;
  int exam = 7; //for debug
  Body cm1;

//definizione tipo per mandare le struct points
  MPI_Datatype Bodytype;
  MPI_Type_contiguous(7,MPI_DOUBLE,&Bodytype);
  MPI_Type_commit(&Bodytype);
  MPI_Aint extent,lb;
  MPI_Type_get_extent(Bodytype,&lb, &extent);

//Importo i corpi del sistema
  if (myrank==0) {
    N = count(input_file);
    //cout << "Dati contati: " << N << "\n";
  }
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  Body * bodies = NULL;
  bodies = new Body[N];
  Body * bodies_end = NULL;
  bodies_end = new Body[N];
  Body * octree = NULL;
  octree = new Body[8*N];
  double g[N][3];
  Body *cm = NULL;
  cm = new Body[8];
  double *ctr = NULL;
  ctr = new double[3];
  int n[Np];
  for(int i=0;i<Np;i++){
    n[i]=0;
  }
  int a;
  if (myrank==0) {
    load(input_file,bodies,N);
    //cout << "Dati caricati" << "\n";
  }

//inizializzo l'octree dividendo i corpi
  if (myrank==0) {
    for(int i=0;i<N;i++){ //ignoro i corpi fuori dal sistema
      if(bodies[i].IsInCube(cc,L)==false){ 
        bodies[i]=bodies[N-1];
        N--;
        //cout << "Corpo espulso" << "\n";
      }
    }
    if(Np>=2){
      //divido in ottanti 
      l = L/2;
      double c[8][3] = {{cc[0] - l/2,cc[1] - l/2,cc[2] - l/2},
                        {cc[0] - l/2,cc[1] - l/2,cc[2] + l/2},
                        {cc[0] - l/2,cc[1] + l/2,cc[2] - l/2},
                        {cc[0] - l/2,cc[1] + l/2,cc[2] + l/2},
                        {cc[0] + l/2,cc[1] - l/2,cc[2] - l/2},
                        {cc[0] + l/2,cc[1] - l/2,cc[2] + l/2},
                        {cc[0] + l/2,cc[1] + l/2,cc[2] - l/2},
                        {cc[0] + l/2,cc[1] + l/2,cc[2] + l/2}}; //posizioni dei centri
      //spedisco ai processi fratelli la visione globale e i centri di loro competenza
      MPI_Scatter(&c, 3, MPI_DOUBLE, ctr, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      for(int i=0;i<N;i++){
        for(int j=0;j<8;j++){
          if(bodies[i].IsInCube(c[j],l)==true){
            octree[j*N+n[j]]=bodies[i];
            n[j]++;
          }
        }
      }
      //cout << "Corpi divisi" << "\n";
    }
    else{
      n[0]=N;
      ctr[0]=cc[0];
      ctr[1]=cc[1];
      ctr[2]=cc[2];
      l = L;
    }
  }
  //spedisco ai processi fratelli la visione globale e i centri di loro competenza
  if(myrank>0 && myrank<8){
    double c; //serve solo per il compilatore
    MPI_Scatter(&c, 3, MPI_DOUBLE, ctr, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  MPI_Bcast(&n, 8, MPI_INT, 0, MPI_COMM_WORLD);  //può servire per il load balancing

  if(myrank==0){
    if(Np>=8){
      for(int i=1;i<8;i++){
        if(n[i]>0){
          MPI_Send(&octree[(i)*N], n[i], Bodytype, i, 0, MPI_COMM_WORLD); //meglio invio ad albero?
        }
      }
          //cout << "Corpi inviati" << "\n";
      double m = 0;
      double p[3] = {0,0,0};  //per il centro di massa
      double vv[3] = {0,0,0}; 
      for(int i=0;i<n[myrank];i++){
        bodies[i]=octree[i];
        m += bodies[i].mass;
        p[0] += bodies[i].coord[0]*bodies[i].mass;
        p[1] += bodies[i].coord[1]*bodies[i].mass;
        p[2] += bodies[i].coord[2]*bodies[i].mass;
      }
      for(int i=0;i<3;i++){
        p[i] = p[i]/m;
      }
      cm1 = Body(m,vv,p);
      //cm1.Print();
      //cout << "Centro di massa calcolato" << "\n";
      MPI_Allgather(&cm1, 1, Bodytype, cm, 1, Bodytype, MPI_COMM_WORLD);
      //cout << "Centri di massa condivisi" << "\n";
    }
    //se ho meno di 8 processi cosa faccio?
  }
  
  if(myrank>0 && myrank<8){
    l = L/2;
    MPI_Recv(bodies, n[myrank], Bodytype, 0, 0, MPI_COMM_WORLD, &status);
    double m = 0;
    double p[3] = {0,0,0};  //per il centro di massa
    double vv[3] = {0,0,0}; 
    for(int i=0;i<n[myrank];i++){
      m += bodies[i].mass;
      p[0] += bodies[i].coord[0]*bodies[i].mass;
      p[1] += bodies[i].coord[1]*bodies[i].mass;
      p[2] += bodies[i].coord[2]*bodies[i].mass;
    }
    for(int i=0;i<3;i++){
        p[i] = p[i]/m;
    }
    cm1 = Body(m,vv,p);
    //cm1.Print();
    MPI_Allgather(&cm1, 1, Bodytype, cm, 1, Bodytype, MPI_COMM_WORLD);
  }  

  if(myrank>=0){

    local_octree = Generate_Octree(n[myrank], ctr, l, bodies);  

    for(int i=0;i<n[myrank];i++){
      g[i][0] = 0; //inizializzo gravità
      g[i][1] = 0;
      g[i][2] = 0;
      //if(Np==1) Compute_Gravity(bodies[i], local_octree, g[i], L);
      //if(Np==8) Compute_Gravity(bodies[i], local_octree, g[i], L, cm, 8, myrank);
    }

    Delete_Octree(local_octree); 

  }

  //ciclo sul tempo
  while(t<=T){
    //Using velocity Verlet algorithm
    //velocity Verlet parte 1
    for(int i=0;i<n[myrank];i++){
      bodies[i].vel[0] = bodies[i].vel[0] + 0.5*dt*g[i][0];
      bodies[i].vel[1] = bodies[i].vel[1] + 0.5*dt*g[i][1];
      bodies[i].vel[2] = bodies[i].vel[2] + 0.5*dt*g[i][2];
      bodies[i].coord[0] = bodies[i].coord[0] + dt*bodies[i].vel[0];
      bodies[i].coord[1] = bodies[i].coord[1] + dt*bodies[i].vel[1];
      bodies[i].coord[2] = bodies[i].coord[2] + dt*bodies[i].vel[2];
    } 

    //scambio bodies tra i processi ed elimino quelli fuoriusciti dal sistema
    if(Np==1){
      vector<int> h;
      for(int i=0;i<n[myrank];i++){
        if(bodies[i].IsInCube(ctr,L)==true){
          h.push_back(i); //salvo l'indice dei bodies da tenere
        }
      }
      n[myrank]=0;
      for(int i=0;i<h.size();i++){
        bodies[n[myrank]]= bodies[h[i]];
        n[myrank]++;
      }
    }
    if(Np==8){
      int n_send;
      int n_recv;
      vector<int> s;
      vector<int> h; 
      Body send[N];
      Body newbodies[N];
      //inizializzo newbodies e n_recv sullo stato attuale
      n_recv = n[myrank];
      for(int i=0;i<n[myrank];i++){
        newbodies[i]=bodies[i];
      }
      n[myrank] = 0; //lo ricalcolo nello scambio
      for(int j=0;j<7;j++){
        for(int i=0;i<n_recv;i++){
          if(newbodies[i].IsInCube(cc,L)==false){
             //bodies fuori dal sistema, non li considero oltre
          }
          else if(newbodies[i].IsInCube(ctr,l)==false){
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
          bodies[n[myrank]]= newbodies[h[i]];
          n[myrank]++;
        }
        
      //spedizioni di scambio a cerchio
        n_send=s.size();
        int id_s = myrank+1;
        int id_r = myrank-1;
        if(myrank==7){
          id_s = 0;
        }
        if(myrank==0){
          id_r = 7;
        }
        if(myrank%2==0){
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
        s.clear();
        h.clear();
      }
    }

    //aggiorno l'octree e la gravità
    local_octree = Generate_Octree(n[myrank], ctr, l, bodies);  

    cm1 = local_octree->cm;
    MPI_Allgather(&cm1, 1, Bodytype, cm, 1, Bodytype, MPI_COMM_WORLD);

    for(int i=0;i<n[myrank];i++){
      g[i][0] = 0; //inizializzo gravità
      g[i][1] = 0;
      g[i][2] = 0;
      //if(Np==1) Compute_Gravity(bodies[i], local_octree, g[i], L);
      //if(Np==8) Compute_Gravity(bodies[i], local_octree, g[i], L, cm, 8, myrank);
    }

    Delete_Octree(local_octree);  

    //velocity Verlet parte 1
    for(int i=0;i<n[myrank];i++){
      bodies[i].vel[0] = bodies[i].vel[0] + 0.5*dt*g[i][0];
      bodies[i].vel[1] = bodies[i].vel[1] + 0.5*dt*g[i][1];
      bodies[i].vel[2] = bodies[i].vel[2] + 0.5*dt*g[i][2]; 
    }    
    t += dt;
  }

//Esporto i risultati
  int shift[8];
  shift[0]=0;
  for(int i=1;i<8;i++){
    shift[i]=shift[i-1]+n[i-1];
  } 
  if (myrank==0) {
    N=0;
    for(int i=0;i<Np;i++){
        N += n[i];
    }
    if(Np>=2){
      MPI_Gatherv(bodies,n[myrank],Bodytype,bodies_end,n,shift,Bodytype,0,MPI_COMM_WORLD);
      save(output_file,bodies_end,N);
    }
    if(Np==1) save(output_file,bodies,N);
    //cout << "Dati salvati" << "\n";
  }
  if (myrank>0 && myrank<8){
    MPI_Gatherv(bodies,n[myrank],Bodytype,bodies_end,n,shift,Bodytype,0,MPI_COMM_WORLD);
  }
  delete octree;
  delete cm;
  delete ctr;
  delete bodies_end;
  t1 = MPI_Wtime();
  MPI_Finalize();
  if(myrank==0){
    cout << "Total time (MPI) " << myrank << " is " << t1-t0 << "\n";
  }
  return 0;
}

//FUNZIONI
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

Barnes* Generate_Octree(int &N, double c[], double Len, Body bodies[]){
    //creo il nodo root
    //cout << "Generating root" << endl;
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
    Barnes* root = new Barnes(m,p);//static Barnes root(m,p);
    if(N>1){
      //divido in ottanti e chiamo il costruttore
      double l = Len/2;
      //salvo il centro madre
      double cc[3] = {c[0],c[1],c[2]};
      //primo ottante
      c[0] = cc[0] - l/2;
      c[1] = cc[1] - l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd1 = new Barnes(N,c,l,bodies);//static Barnes nd1(N,c,l,1);
      if(nd1->empty==true){
          delete nd1;
      }
      else{
          root->subnodes.push_back(nd1);
      }
      //secondo ottante
      c[2] = cc[2] + l/2;
      Barnes* nd2 = new Barnes(N,c,l,bodies);//static Barnes nd2(N,c,l,1);
      if(nd2->empty==true){
          delete nd2;
      }
      else{
          root->subnodes.push_back(nd2);
      }
      //terzo ottante
      c[1] = cc[1] + l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd3 = new Barnes(N,c,l,bodies);//static Barnes nd3(N,c,l,1);
      if(nd3->empty==true){
          delete nd3;
      }
      else{
          root->subnodes.push_back(nd3);
      }
      //quarto ottante
      c[2] = cc[2] + l/2;
      Barnes* nd4 = new Barnes(N,c,l,bodies);//static Barnes nd4(N,c,l,1);
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
      Barnes* nd5 = new Barnes(N,c,l,bodies);//static Barnes nd5(N,c,l,1);
      if(nd5->empty==true){
          delete nd5;
      }
      else{
          root->subnodes.push_back(nd5);
      }
      //sesto ottante
      c[2] = cc[2] + l/2;
      Barnes* nd6 = new Barnes(N,c,l,bodies);//static Barnes nd6(N,c,l,1);
      if(nd6->empty==true){
          delete nd6;
      }
      else{
          root->subnodes.push_back(nd6);
      }
      //settimo ottante
      c[1] = cc[1] + l/2;
      c[2] = cc[2] - l/2;
      Barnes* nd7 = new Barnes(N,c,l,bodies);//static Barnes nd7(N,c,l,1);
      if(nd7->empty==true){
          delete nd7;
      }
      else{
          root->subnodes.push_back(nd7);
      }
      //ottavo ottante
      c[2] = cc[2] + l/2;
      Barnes* nd8 = new Barnes(N,c,l,bodies);//static Barnes nd8(N,c,l,1);
      if(nd8->empty==true){
          delete nd8;
      }
      else{
          root->subnodes.push_back(nd8);
      }
    }
	return root;
}

void Delete_Octree(Barnes* root){//gestione allocazione dinamica
	root->Delete_Subnodes();
	delete root;
}

double Distance(Body b1, Body b2){

	double a[3] = {b1.coord[0],b1.coord[1],b1.coord[2]};
	double b[3] = {b2.coord[0],b2.coord[1],b2.coord[2]};

    double d = sqrt(abs(((a[0]-b[0])*(a[0]-b[0]))+
                    ((a[1]-b[1])*(a[1]-b[1]))+
                    ((a[2]-b[2])*(a[2]-b[2]))));
    
    if(isnan(d)==1){
        cout << "Nan rilevato" << "\n";
        b1.Print();
        b2.Print();
    }

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
  //cout << d << "\n";
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
        //cout << d << "\n";
        l = Len/(2);
        if(d<=l && current->subnodes.size()>0){
          //cout << "Going down" << "\n";
          //cout << g[0] << "\n";
			    Compute_Gravity(body,current,g,l);
          //cout << g[0] << "\n";
        }
        else if(d<=pow(10,-16)){ //escludo i casi di distanze tra un corpo e se stesso con precisione double
          //cout << d << "\n";
        }
        else{
            m = current->cm.mass;
            Vectorial_Distance(body,current->cm,dx);
            acc(m,d,dx,g);  //aggiungo gravità 
            //cout << g[0] << "\n";
          }
      }
}

void Compute_Gravity(Body body, Barnes* current, double *g, double Len, Body nodes_cm[], int N_nodes, int self){
  double dx[3]; //per la distanza vettoriale tra corpi
  double m;
  double l;
	double d;
  for(int k=0;k<current->subnodes.size();k++){
    current=current->subnodes[k];
    d = Distance(body,current->cm);
    //cout << d << "\n";
    l = Len/(2);
    if(d<=l && current->subnodes.size()>0){
          //cout << "Going down" << "\n";
          //cout << g[0] << "\n";
			Compute_Gravity(body,current,g,l);
          //cout << g[0] << "\n";
    }
    else if(d<=pow(10,-16)){ //escludo i casi di distanze tra un corpo e se stesso con precisione double
          //cout << d << "\n";
    }
    else{
      m = current->cm.mass;
      Vectorial_Distance(body,current->cm,dx);
      acc(m,d,dx,g);  //aggiungo gravità 
            //cout << g[0] << "\n";
    }
  }
  for(int i=0;i<N_nodes;i++){
    if(i!=self){
      m = nodes_cm[i].mass;
      Vectorial_Distance(body,nodes_cm[i],dx);
      acc(m,d,dx,g);
    }
  }
}