#include "Functions.h"
#include "Barnes.h"

//funzione per sapere quanti corpi sono descritti dal file
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

//funzione per leggere il file con un array di corpi
vector<Body> load(string file, int dimension){
	vector<Body> vector;
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
				//vector[i].SetMass(stof(m));
				double x[3] = {stod(p[0]),stod(p[1]),stod(p[2])};
				double vl[3] = {stod(v[0]),stod(v[1]),stod(v[2])};
				//vector[i].SetCoordinates(x);
				//vector[i].SetVelocity(vl);
				Body bd(stof(m),vl,x);
				vector.push_back(bd);
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
	return vector;
}

//funzione per salvare su file i valori di un array di corpi
void save(string file, vector<Body> vector, int dimension){
	ofstream txt;
	txt.open(file);
	double x[3];
	double v[3];
	for(int i=0;i<dimension;i++){
		vector[i].GetCoordinates(x);
		vector[i].GetVelocity(v);
		txt << setiosflags(ios::scientific);
		txt << vector[i].GetMass() << " ";
		txt << x[0] << " ";
		txt << x[1] << " ";
		txt << x[2] << " ";
		txt << v[0] << " ";
		txt << v[1] << " ";
		txt << v[2] << endl;
	}
	txt.close();
}

//Funzione Distance per il calcolo di distanze 3D
double Distance(double *a, double *b){

    double d = sqrt(((a[0]-b[0])*(a[0]-b[0]))+
                    ((a[1]-b[1])*(a[1]-b[1]))+
                    ((a[2]-b[2])*(a[2]-b[2])));

    return d;
}

double Distance(Body b1, Body b2){

	double a[3];
	double b[3];
	b1.GetCoordinates(a);
	b2.GetCoordinates(b);

    double d = sqrt(((a[0]-b[0])*(a[0]-b[0]))+
                    ((a[1]-b[1])*(a[1]-b[1]))+
                    ((a[2]-b[2])*(a[2]-b[2])));

    return d;
}

void Vectorial_Distance(Body b1, Body b2, double *d){

	double a[3];
	double b[3];
	b1.GetCoordinates(a);
	b2.GetCoordinates(b);

    d[0] = b[0]-a[0];
	d[1] = b[1]-a[1];
	d[2] = b[2]-a[2];
}

//Per la creazione dell'Octree
Barnes* Generate_Octree(int &N, double c[], double L){
    if(bodies.size()==1){
        cout << "Error: single body problem!" << endl;
		exit(0);
    }
    if(bodies.size()==0){
        cout << "Error: no bodies!" << endl;
		exit(0);
    }
    //creo il nodo root
    double m = 0;
    double mm;
    double x[3]; //per le coordinate dei corpi
    double p[3] = {0,0,0};  //per il centro di massa
    int nb = 0;
    for(int i=0;i<N;i++){
        if(bodies[i].IsInCube(c,L)==true){
            nb++;
            bodies[i].GetCoordinates(x);
            mm = bodies[i].GetMass();
            m += mm;
            p[0] += x[0]*mm;
            p[1] += x[1]*mm;
            p[2] += x[2]*mm;
        }
        else{
            bodies.erase(bodies.begin()+i);
            N--;
            cout << "A body left the system" << endl;
			if(bodies.size()<3){
        		cout << "Error: not enough bodies!" << endl;
				exit(0);
    		}
        }
    }
    p[0]=p[0]/m;
    p[1]=p[1]/m;
    p[2]=p[2]/m;
    Barnes* root = new Barnes(m,p);//static Barnes root(m,p);
    //cout << "Generating root" << endl;
    //divido in ottanti e chiamo il costruttore
    double l = L/2;
    //salvo il centro madre
    double cc[3] = {c[0],c[1],c[2]};
    //primo ottante
    c[0] = cc[0] - l/2;
    c[1] = cc[1] - l/2;
    c[2] = cc[2] - l/2;
    Barnes* nd1 = new Barnes(N,c,l,1);//static Barnes nd1(N,c,l,1);
	if(nd1->empty==true){
        delete nd1;
    }
    else{
        root->subnodes.push_back(nd1);
    }
    //secondo ottante
    c[2] = cc[2] + l/2;
    Barnes* nd2 = new Barnes(N,c,l,1);//static Barnes nd2(N,c,l,1);
	if(nd2->empty==true){
        delete nd2;
    }
    else{
        root->subnodes.push_back(nd2);
    }
    //terzo ottante
    c[1] = cc[1] + l/2;
    c[2] = cc[2] - l/2;
    Barnes* nd3 = new Barnes(N,c,l,1);//static Barnes nd3(N,c,l,1);
	if(nd3->empty==true){
        delete nd3;
    }
    else{
        root->subnodes.push_back(nd3);
    }
    //quarto ottante
    c[2] = cc[2] + l/2;
    Barnes* nd4 = new Barnes(N,c,l,1);//static Barnes nd4(N,c,l,1);
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
    Barnes* nd5 = new Barnes(N,c,l,1);//static Barnes nd5(N,c,l,1);
	if(nd5->empty==true){
        delete nd5;
    }
    else{
        root->subnodes.push_back(nd5);
    }
    //sesto ottante
    c[2] = cc[2] + l/2;
    Barnes* nd6 = new Barnes(N,c,l,1);//static Barnes nd6(N,c,l,1);
	if(nd6->empty==true){
        delete nd6;
    }
    else{
        root->subnodes.push_back(nd6);
    }
    //settimo ottante
    c[1] = cc[1] + l/2;
    c[2] = cc[2] - l/2;
    Barnes* nd7 = new Barnes(N,c,l,1);//static Barnes nd7(N,c,l,1);
	if(nd7->empty==true){
        delete nd7;
    }
    else{
        root->subnodes.push_back(nd7);
    }
    //ottavo ottante
    c[2] = cc[2] + l/2;
    Barnes* nd8 = new Barnes(N,c,l,1);//static Barnes nd8(N,c,l,1);
	if(nd8->empty==true){
        delete nd8;
    }
    else{
        root->subnodes.push_back(nd8);
    }

	return root;
}

void Delete_Octree(Barnes* root){//gestione allocazione dinamica
	root->Delete_Subnodes(1);
	//cout << "Delete root" << endl;
	delete root;
}

void acc(double m, double d, double *x, double *a){
  double mr = m/(d*d*d);
  a[0] += x[0]*mr;
  a[1] += x[1]*mr;
  a[2] += x[2]*mr;
}

void Compute_Gravity(int i, Barnes* current, double *g, double L){
	double dx[3]; //per la distanza vettoriale tra corpi
    double m;
    double l;
	double d;
    for(int k=0;k<current->subnodes.size();k++){
        current=current->subnodes[k];
        d = Distance(bodies[i],current->cm);
        l = L/(2);
        if(d<=l && current->subnodes.size()>0){
			Compute_Gravity(i,current,g,l);
        }
        else{
            m = current->cm.GetMass();
            Vectorial_Distance(bodies[i],current->cm,dx);
            acc(m,d,dx,&g[i]);  //aggiungo gravità 
          }
      }
}