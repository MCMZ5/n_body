#include "Barnes.h"

using namespace std;

Barnes::Barnes(){
    empty = false;
}

Barnes::Barnes(double m, double *x){
    cm.SetMass(m);
    cm.SetCoordinates(x);
    empty = false;
}

Barnes::Barnes(Body &b){
    cm = b;
    empty = false;
}

Barnes::Barnes(int &N, double c[], double L, int lvl){
    //creo il nodo base
    double m = 0;
    double mm;
    double x[3]; //per le coordinate dei corpi
    double p[3] = {0,0,0};  //per il centro di massa
    int nb = 0;
    int lf;
    for(int i=0;i<N;i++){
        if(bodies[i].IsInCube(c,L)==true){
            nb++;
            if(nb==1){
                lf=i;
            }
            bodies[i].GetCoordinates(x);
            mm = bodies[i].GetMass();
            m += mm;
            p[0] += x[0]*mm;
            p[1] += x[1]*mm;
            p[2] += x[2]*mm;
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
            cm.SetMass(m);
            cm.SetCoordinates(p);
            empty = false;
        }
        //cout << "Generating node of level " << lvl << endl;
    }
    if(nb>1){
        //divido in ottanti e chiamo il costruttore
        double l = L/2;
        //salvo il centro madre
        double cc[3] = {c[0],c[1],c[2]};
        //primo ottante
        c[0] = cc[0] - l/2;
        c[1] = cc[1] - l/2;
        c[2] = cc[2] - l/2;
        Barnes* nd1 = new Barnes(N,c,l,lvl+1);//static Barnes nd1(N,c,l,lvl+1);
        if(nd1->empty==true){
            delete nd1;
        }
        else{
            subnodes.push_back(nd1);
        }
        //secondo ottante
        c[2] = cc[2] + l/2;
        Barnes* nd2 = new Barnes(N,c,l,lvl+1);//static Barnes nd2(N,c,l,lvl+1);
        if(nd2->empty==true){
            delete nd2;
        }
        else{
            subnodes.push_back(nd2);
        }
        //terzo ottante
        c[1] = cc[1] + l/2;
        c[2] = cc[2] - l/2;
        Barnes* nd3 = new Barnes(N,c,l,lvl+1);//static Barnes nd3(N,c,l,lvl+1);
        if(nd3->empty==true){
            delete nd3;
        }
        else{
            subnodes.push_back(nd3);
        }
        //quarto ottante
        c[2] = cc[2] + l/2;
        Barnes* nd4 = new Barnes(N,c,l,lvl+1);//static Barnes nd4(N,c,l,lvl+1);
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
        Barnes* nd5 = new Barnes(N,c,l,lvl+1);//static Barnes nd5(N,c,l,lvl+1);
        if(nd5->empty==true){
            delete nd5;
        }
        else{
            subnodes.push_back(nd5);
        }
        //sesto ottante
        c[2] = cc[2] + l/2;
        Barnes* nd6 = new Barnes(N,c,l,lvl+1);//static Barnes nd6(N,c,l,lvl+1);
        if(nd6->empty==true){
            delete nd6;
        }
        else{
            subnodes.push_back(nd6);
        }
        //settimo ottante
        c[1] = cc[1] + l/2;
        c[2] = cc[2] - l/2;
        Barnes* nd7 = new Barnes(N,c,l,lvl+1);//static Barnes nd7(N,c,l,lvl+1);
        if(nd7->empty==true){
            delete nd7;
        }
        else{
            subnodes.push_back(nd7);
        }
        //ottavo ottante
        c[2] = cc[2] + l/2;
        Barnes* nd8 = new Barnes(N,c,l,lvl+1);//static Barnes nd8(N,c,l,lvl+1);
        if(nd8->empty==true){
            delete nd8;
        }
        else{
            subnodes.push_back(nd8);
        }
    }
}

Barnes::~Barnes(){

}

void Barnes::Delete_Subnodes(int lvl){
    if(subnodes.size()>0){
        for(int i=0;i<subnodes.size();i++){
            subnodes[i]->Delete_Subnodes(lvl+1);
            //cout << "Deleting node of level " << lvl << endl;
            delete subnodes[i];
        }
    }
}