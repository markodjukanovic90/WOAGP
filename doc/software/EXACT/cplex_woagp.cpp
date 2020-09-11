/* the CPLEX modes for a specific 2-LCIPS problem */
/* before compilation makes sure that CPLEX software path is correct */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "Timer.h"
#include <vector>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <vector>
#include <set>
#include <limits>
#include <unordered_map>
#include <map>
#include <ilcp/cp.h>


// the following "include" is necessary for the correct working/compilation of CPLEX. You should adapt this path to your installation of CPLEX
#include "/home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/include/ilcplex/ilocplex.h"
typedef pair<float, float> Point_2;
/** pomocne strukture **/

vector<float> Surface; // for a vertex indexed by i, we return the surface of V(i)
vector<float> Cost; // at index i --> w_i
vector<vector<float>>Intersection; // Intersection[i][j] = Volume of V(i) intersected by V(j)
vector<set<Point_2>> S; // i-th position of S is a set of all points that are visible from vertex i
set<Point_2> D_P;
vector<Point_2> Vertices; // vector of vertices
#define INFEASIBLE 1000000
int cardinalityD; // |D(P)|
float wTotal = 0.0; // total weight
int w_type = 0; // ti ptezine koju pozivamo; 0: tezina proporcionalna velicini vidljivosti svakog vrha; 1: --; 2: --

using namespace std; 

ILOSTLBEGIN

// global variables concerning the random number generator
time_t t;

// vector for keeping the names of the input files
string inputFile;
int sigma = 0;
// time limit for CPLEX (can be changed at run time via the -t comand line parameter)
double time_limit = 900.0;
string outdir="/home/djukanovic/Desktop/projects/rflcsp/MIP-RFLCS-2d/CP-instances/"; // the directory where .out files will be saved
int m = 0; // number of sets
int n = 0; // cardnality of X, i.e., |X|
int alg = 0; // koji algoritam pozivamo 0: CPLEX (WSC); 1: Greedy (WSC)

// end pomocne strukture

inline int stoi(string &s) {

       return atoi(s.c_str());
}

inline double stof(string &s) {

       return atof(s.c_str());
}

float dist(int i, int j)
{
     return sqrt( pow( Vertices[i].first - Vertices[j].first, 2) + pow( Vertices[ i ].second - Vertices[ j ].second, 2) );

}

void read_from_file(std:: string path)
{
     // TODO: read from file
     std::ifstream infile(path);
     std::string line;
     std::getline(infile, line);
     int n = stoi(line); // number of vertices

    int temp = 0;

    while(temp<n){
        std::getline(infile, line);

        std::string delim1 = ":";
         auto start = 0U;
         auto end = line.find(delim1);

         //izdvajanje koordinata tjemena (vrha) poligona
         std:: string  vrh =  line.substr(start, end - start);

         std::string delim2 = ",";

         auto start1 = 0U;
         auto end1 = vrh.find(delim2);

         std:: string  vrh1 =  vrh.substr(start1+1, end1 - start1-1);
         std:: string  vrh2 =  vrh.substr(end1+1, vrh.length()-end1-2);


         Point_2 vertex;
         vertex.first = stof(vrh1);
         vertex.second = stof(vrh2);
         D_P.insert(vertex);
         Vertices.push_back(vertex);

        //izdvjanje skupa tacaka koje se vide iz tog vrha

        std:: string ostatak = line.substr(end+1, line.length() - end);

        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;
        set<Point_2> vertexSet;

        while ((pos = ostatak.find(delimiter)) != std::string::npos) {
            token = ostatak.substr(0, pos);
            std::cout << token << std::endl;

            //obrada svakog tokena da bi se dobile koordinate tacaka


         std::string delim3 = ",";

         auto start3 = 0U;
         auto end3 = vrh.find(delim3);

         std:: string  x =  token.substr(start3+1, end3 - start3-1);
         std:: string  y =  token.substr(end3+1, token.length()-end3-2);

         Point_2 tacka;
         tacka.first = stof(x);
         tacka.second = stof(y); 
         D_P.insert(tacka);

        vertexSet.insert(tacka);
        ostatak.erase(0, pos + delimiter.length());
    }
    S.push_back(vertexSet);
    temp++;
   }
    //izdvajanje povrsina
    temp = 0;
    while(temp<n){
        std::getline(infile, line);

        std::string delim1 = ":";
         auto start = 0U;
         auto end = line.find(delim1);

        std:: string ostatak = line.substr(end+1, line.length() - end);

        float povrsina = stof(ostatak);
         Surface.push_back(povrsina);
        temp++;

    }


   //izdvjanje presjeka

   temp = 1;
   vector<float> tempVect;
    while(std::getline(infile, line))
    {

         std::string delim1 = ":";
         auto start = 0U;
         auto end = line.find(delim1);

          std:: string  vrhovi =  line.substr(start, end - start);
          std:: string ostatak = line.substr(end+1, line.length() - end);



        float povrsinaPresjeka = stof(ostatak);


        if(temp<n){
            tempVect.push_back(povrsinaPresjeka);
            temp++;
        }

        if(temp==n){
            Intersection.push_back(tempVect);
            tempVect.clear();
            temp=1;
            n--;
        }
    }
    //cout << S.size() << " n: " << n << " Surface " << Surface.size() << endl;
}


void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-f")==0) inputFile = argv[++iarg];
        else if (strcmp(argv[iarg],"-t")==0) time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-alg")==0) alg = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-w_type")==0) w_type = atof(argv[++iarg]);
        iarg++;
    } 
     cout << "w_type " << w_type << endl;
}

// write during the CPLEX call
ILOMIPINFOCALLBACK6(loggingCallback,
                    Timer&, timer,
                    vector<double>&, times,
                    vector<double>&, results,
                    vector<double>&, gaps,
                    int, iter,
                    IloNum,         lastIncumbent) {

    IloNum nv = getIncumbentObjValue();
    double newTime = timer.elapsed_time(Timer::VIRTUAL);
    double newGap = 100.0*getMIPRelativeGap();
    if (nv < lastIncumbent) {
        cout << "CPLEX sol: " << nv << "\ttime " << newTime <<  "\tgap " << newGap << endl;
        results[iter] = nv;
        times[iter] = newTime;
        gaps[iter] = newGap;
    }
    lastIncumbent = nv;
}

void run_cplex(vector<set<Point_2>>& S, set<Point_2>& D_P, int n){

  IloEnv env;  cout << "run_cplex " << endl;
  env.setOut(env.getNullStream());
  try
  {
   IloModel model(env);
   IloObjective obj = IloMinimize(env);
   // defining the set of binary variables Z
   vector<IloNumVar> Z; cout << "S.size() ----> " << S.size() << n << endl;
   for(int i = 0; i < n; ++i){
       //cout << "cost i " << Cost[i] << endl;
       IloNumVar myIntVar(env, 0, 1, ILOINT);
       Z.push_back(myIntVar); // x_i  \in {0, 1}
       obj.setLinearCoef(Z[i], Cost[i]); // sum_i c_i x_i
   }
   // constraints 
   cout << "#Vars: " << Z.size() << endl;
   // constraints 
   int index = 0;
   for(Point_2 p : D_P) 
   {
       IloExpr expr_i(env); bool add = false;
       for(int j = 0; j < n; ++j)   
       {
           if( S[j].count( p ) ){ // item i se nalazi u skupu j
               expr_i += Z[j];
               add = true;
           }
       }
       if(add)
          model.add(expr_i >= 1);      
       index++ ;
   }  

   //std::cout << "#Constraints: "<< constraints_num << std::endl;
   //solve the modelIloLinearNumExpr objective = cplex.linearNumExpr();
   model.add(obj);
   IloCplex cplex(model);
   //int time_limit = 900;
   // pass the time limit to CPLEX
   cplex.setParam(IloCplex::TiLim, time_limit);
   // the following two parameters should always be set in the way as shown
   cplex.setParam(IloCplex::NodeFileInd, 2);
   cplex.setParam(IloCplex::Threads, 1);
   IloNum lastObjVal = std::numeric_limits<double>::min();
   // tell CPLEX to make use of the function 'loggingCallback' for writing out information to the screen
   //cplex.use(loggingCallback(env, timer, times, results, gaps, iter, lastObjVal));
   cout << "CPLEX pokrenut" << endl;
   // cplex.exportModel("ws.lp");
   cplex.solve();
  
   if (cplex.getStatus() == IloAlgorithm::Optimal or cplex.getStatus() == IloAlgorithm::Feasible)
   {
       if(cplex.getStatus() == IloAlgorithm::Optimal)
          cout << "CPLEX finds optimal" << endl;
       else
          cout << "CPLEX finds feasible solution" << endl;

       double lastVal = double(cplex.getObjValue());
       cout << "\n CPLEX sol: " << lastVal << endl;
       // print the objective point
       cout << "Sets in the solution: {" <<endl;
       bool first = true;  
       for(int i = 0; i < Z.size(); ++i){
           IloNum xval = cplex.getValue(Z[i]);
           // the reason for 'xval > 0.9' instead of 'xval == 1.0' will be explained in class
           if (xval > 0.9) {
               cout << "S_" << ( i + 1 ) << ", ";
               //myfileOut << (*it).first;
           }

       }
       cout << "}" << endl;
       //myfileOut << "value: " << lastVal << endl;
   }
   else cout << "Nema rjesenja" << endl;
}
catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
}
env.end();
}



int run_cp(vector<set<Point_2>>& S, set<Point_2>& D_P, int n){

  IloEnv env; cout << "Run CP" << endl;
  env.setOut(env.getNullStream());
  try
  {
   // define the LCIS problem model...
   IloModel model(env);
   IloObjective obj = IloMinimize(env);
   // defining the set of binary variables Z
   vector<IloNumVar> Z; // cout << "S.size() " << S.size() << endl;
   for(int i = 0; i < n; ++i){
       //cout << "w_i: " << Cost[i] << endl;
       IloNumVar myIntVar(env, 0, 1, ILOINT);
       Z.push_back(myIntVar); // x_i  \in {0, 1}
       obj.setLinearCoef(Z[i], Cost[i]); // sum_i c_i x_i
   }
   cout << "#Vars: " << Z.size() << " n=" << n << " " << S.size() << endl;
   // constraints 
   int index = 0;
   for(Point_2 p : D_P) 
   {
       IloExpr expr_i(env); bool add = false;
       for(int j = 0; j < n; ++j)   
       {
           if( S[j].count( p )){ // item i se nalazi u skupu j
               //cout << "j: " << j << endl;
               expr_i += Z[j];
               add = true;
           }
       }
       if(add)
          model.add(expr_i >= 1);      
        index++;
   }  
   model.add(obj);
   IloCP cp(model);
   cp.setParameter(IloCP::FailLimit, 3000000);
   // the following two parameters should always be set in the way as shown

   if (cp.solve())
   {
       double lastVal = double(cp.getObjValue());
       cout << "Cp sol: " << lastVal << endl;
       // print the objective point
       cout << "Sets in the solution: {" <<endl;
       bool first = true;  
       for(int i = 0; i < Z.size(); ++i){
           IloNum xval = cp.getValue(Z[i]);
           // the reason for 'xval > 0.9' instead of 'xval == 1.0' will be explained in class
           if (xval > 0.9) {
               cout << "S_" << ( i + 1 ) << ", ";
               //myfileOut << (*it).first;
           }

       }
       cout << "}" << endl;
   }
 }catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
}
env.end();
}


void fill_cost()
{
// ----------------------- dodjela tezina (proporcionalno broju pokrivenih tacaka iz D(P) ------------------------------------
    cout << "fill cost " << w_type << endl;
    switch(w_type){
 
          case 0: { // tezina proporcionalna sa velicinom skupa S[i]
                  for(auto& X: S){
                      cardinalityD += X.size();
                  }
                  for (size_t i = 0; i < n; i++){
                      float w_i = n * n * ( ((float) S[i].size()) / cardinalityD );
                      wTotal += w_i;cout << "w_i " << w_i << endl;
                      Cost.push_back(w_i); 
                  }; break;}
          case 1: {     
                  float dist0  = (dist(0, n-1) + dist(0,1) ) / 2; 
                  Cost.push_back(dist0); wTotal += dist0;
                   for(int i = 1; i < n; ++i){ 
                      float w_i = ( (0.0 + dist(i-1, i) + dist(i, i+1) ) / 2 );
                      Cost.push_back( w_i );
                      wTotal += w_i;
                      //cout << "w_i= " << Cost[i] << endl;
                   }
                     break;}
          default: {
                    for(size_t i = 0; i < n; i++) // non-weighted version of the problem
                       Cost.push_back(1);
                   }
   }


}
/**********
Main function
**********/


template<typename T>
void print_vector(vector<T>& vec)
{
     for(const T& c: vec)
         cout << c <<" ";
     cout << endl;  
}

int main(int argc, char **argv ) {
        
        cout << " Reading parameters " << endl;
        read_parameters(argc,argv);
        //setting the output format for doubles to 2 decimals after the comma
        std::cout << std::setprecision(2) << std::fixed;
        Timer timer;
        // generisanje stream-a
        double cur_time = timer.elapsed_time(Timer::VIRTUAL);
        // opening the corresponding input file and reading the problem data
        read_parameters(argc, argv);
        read_from_file(inputFile); n = S.size();// fill Cost, Intersection and Surface
        fill_cost(); cout << "cost--->" << Cost.size() << " n: " << n << endl;
        if(alg == 0) 
        { 
            cout << "Run CPLEX..." << endl;
            run_cplex(S, D_P, n); //, myfileOut);
            double end_time = timer.elapsed_time(Timer::VIRTUAL);
            cout << "CPLEX time:" << (end_time - cur_time ) <<"\n";
            // greedy procedure:
        }else{
        
            Timer timer_cp;
            // generisanje stream-a
            double cur_time_cp = timer_cp.elapsed_time(Timer::VIRTUAL);
            run_cp(S, D_P, n);; 
            double end_time_cp = timer_cp.elapsed_time(Timer::VIRTUAL);
            cout <<"time: " << (end_time_cp - cur_time_cp) << endl;
            //myfileOut << "time: " << (end_time - cur_time ) <<"\n";*/
            //myfileOut.close();
        }
}



