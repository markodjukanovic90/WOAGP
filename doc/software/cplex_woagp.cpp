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

// the following "include" is necessary for the correct working/compilation of CPLEX. You should adapt this path to your installation of CPLEX
#include "/home/djukanovic/Desktop/projects/LCAPS_software/cplex-12.5/include/ilcplex/ilocplex.h"

using namespace std; 

ILOSTLBEGIN

// global variables concerning the random number generator
time_t t;

// vector for keeping the names of the input files
vector<string> inputFile;
int sigma = 0;
// time limit for CPLEX (can be changed at run time via the -t comand line parameter)
double time_limit = 900.0;
string outdir="/home/djukanovic/Desktop/projects/rflcsp/MIP-RFLCS-2d/CP-instances/"; // the directory where .out files will be saved
int m = 0; // number of sets
int n = 0; // cardnality of X, i.e., |X|
int alg = 0; // koji algoritam pozivamo 0: CPLEX (WSC); 1: Greedy (WSC)

inline int stoi(string &s) {

       return atoi(s.c_str());
}

inline double stof(string &s) {

       return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-i")==0) inputFile.push_back(argv[++iarg]);
        else if (strcmp(argv[iarg],"-t")==0) time_limit = atof(argv[++iarg]);
        else if (strcmp(argv[iarg],"-alg")==0) alg = atof(argv[++iarg]);
        iarg++;
    }
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

int f(vector<set<int>>& C)
{
       std::set<int> d;
       for(set<int> di : C){
     
           for (int dii: di){
               //cout << "dii" << dii << endl;
               d.insert(dii);
           }
        
       }
       //cout <<"D" << d.size() << endl; 
       //for(int diii: d)
       //  cout << diii << "\t" <<endl;

       return (int)d.size();
}

int f_minus(vector<set<int>>& C, set<int>& ss ) // f(C u {s}) - f(S)
{
    int count = 0; // ako i \in ss nije niti u jednom od skupova u S, onda count++;
    if(C.empty())
       return ss.size();

    for(int i: ss) 
    {
        bool presented = false;
        for(set<int>& di : C){
            if(di.count(i) > 0){
               presented = true;
               break; 
            }
        }
        if(presented) 
           count++;
    }
    return ss.size() - count; 
}

float greedy_criterion(vector<int>& Cost, vector<set<int>>& S, int i, vector<set<int>>& C) // take s_i from S
{     
    set<int> s = S[i];
    int f_m = f_minus(C, s);  //cout << "f_m: " << f_m << endl;
    if(f_m == 0)
       return 100000;

    float val = ((float)Cost[i]) / (f_minus(C, s));     //cout << "val: " << val << endl; 
    return val;
}

bool findA(vector<int>& s, int a)
{

     for(int si : s)
         if(si == a)
            return true;

      return false;  
}

int min_greedy(vector<int>& Cost, vector<set<int>>& S, vector<set<int>>& C, vector<int>& indeks)
{         //cout << "min_greedy" << endl;
          float g_m = 10000000; int dodaj;
          for(int i = 0; i < S.size(); ++i)
          {  // cout << "i " << i << endl;
              float g_mi = greedy_criterion(Cost, S, i, C); //cout << "gmi: " << g_mi << endl;
              if(g_mi < g_m and g_mi != 100000 and !findA(indeks, i)) 
              { 
                 dodaj = i;
                 g_m = g_mi;   
              }
          }//cout << "dodati....." << dodaj << endl;
          return dodaj;
}


float greedy_procedure(vector<set<int>>& S, vector<int>& Cost)
{
     vector<set<int>> C; vector<int> indeks;
     int f_S = f(S); //cout << "f_S: " << f_S << endl;
     int f_C = 0;
     while(f_C != f_S) 
     {
           int index_set = min_greedy( Cost, S, C, indeks);  
           
           if(!findA(indeks, index_set)){ //jos nije dodan
               C.push_back(S[index_set]);cout << "dodaj ----> " << index_set << endl;
              indeks.push_back(index_set);
           } 
           f_C = f(C); 
           cout << "f_C=" << f_C << " |C|=" << C.size() << endl;
     }
     float greedy_val = 0;
     for(int ix: indeks){
         greedy_val += Cost[ ix ];
     }
     return greedy_val; 
     //return C.size();
}


void run_cplex(vector<int>& C, vector<set<int>>& S){

  IloEnv env;
  env.setOut(env.getNullStream());
  try
  {
   // define the LCIS problem model...
   IloModel model(env);
   IloObjective obj = IloMinimize(env);
   // defining the set of binary variables Z
   vector<IloNumVar> Z; cout << "S.size() " << S.size() << endl;
   for(int i = 0; i < S.size(); ++i){

       IloNumVar myIntVar(env, 0, 1, ILOINT);
       Z.push_back(myIntVar); // x_i  \in {0, 1}
       obj.setLinearCoef(Z[i], C[i]); // sum_i c_i x_i
   }

   cout << "#Vars: " << Z.size() << " n=" << n << " " << S.size() << endl;
   // constraints 
   for(int i = 0; i < n; ++i) 
   {
       IloExpr expr_i(env); bool add = false;
       for(int j = 0; j < S.size(); ++j)   
       {
           if( S[j].count(i+1) ){ // item i se nalazi u skupu j
               expr_i += Z[j];
               add = true;
           }
       }
       if(add)
          model.add(expr_i >= 1);      
       
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
       cout << "\n CPLEX sol: " << lastVal << endl;
       //myfileOut << "value: " << lastVal << endl;
   }
   else cout << "Nema rjesenja" << endl;
}
catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
}
env.end();
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

	vector< set<int> > S;
	vector <int> C; // randomly select cost of each set S_i
         
	ifstream input;
	input.open(inputFile[0].c_str());
        //input >> m >> n; cout << "m=" << m << " n=" << n << endl;
        //std::string str1;
        //std::getline(input, str1);  
        std::string str;
        string content; int x;
        int line = 0;
        while(std::getline(input, content)) {
              //cout << content << ' ';
              if(line >= 1) // prva linija je vec procitana
              {   
                 std::stringstream   linestream(content);
                 set<int> S_i;
                 while(linestream >> x){
                    //cout << x << endl;
                    S_i.insert(x);
                 }
                 if(line % 2 == 0){ // svaka druga linija ocitava jedan skup:
                     S.push_back(S_i);
                     S_i.clear();
                     int k = 1 + rand() % 10;
                     C.push_back(k);    // random cijena od [0,99]
                     //cout << k << endl;
                 }
              }
              else{
                  std::stringstream   linestream(content);
                  linestream >> m >> n; cout <<"n---> " << n << endl;
              }
              line++;
        }
	
        input.close();
        cout << "End of input " << inputFile[0].c_str() << ".out" << endl;
        cout << "|S|=" << S.size() <<endl; 

        for(int i = 0; i < S.size(); ++i){
            cout << "read S_" << i << endl;
            for(set<int>::iterator it = S[i].begin(); it != S[i].end(); ++it){
                cout << (*it) << endl;
            }
         }
         // write the results in the output...
         std::stringstream filenameOut;
         filenameOut << outdir << inputFile[0].c_str() << ".out";
         if(alg == 0) 
         { 
            cout << "Run CPLEX..." << endl;
            run_cplex(C, S); //, myfileOut);
            double end_time = timer.elapsed_time(Timer::VIRTUAL);
            cout << "CPLEX time:" << (end_time - cur_time ) <<"\n";
            // greedy procedure:
         }else{
        
            Timer timer_greedy;
            // generisanje stream-a
            double cur_time_greedy = timer_greedy.elapsed_time(Timer::VIRTUAL);
            float greedy_sol = greedy_procedure( S, C); 
            cout <<"GREEDY sol: " << greedy_sol << endl;
            double end_time_greedy = timer_greedy.elapsed_time(Timer::VIRTUAL);
            cout <<"GREEDY time: " << (end_time_greedy - cur_time_greedy) << endl;
            //myfileOut << "time: " << (end_time - cur_time ) <<"\n";*/
            //myfileOut.close();
         }
}



