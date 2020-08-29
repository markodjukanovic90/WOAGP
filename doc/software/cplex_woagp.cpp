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
        cout << "value " << nv << "\ttime " << newTime <<  "\tgap " << newGap << endl;
        results[iter] = nv;
        times[iter] = newTime;
        gaps[iter] = newGap;
    }
    lastIncumbent = nv;
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
           //else
           //    expr_i += 0;//Z[j];
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
   cplex.exportModel("ws.lp");
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
       cout << "\nvalue: " << lastVal << endl;
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
        
         cout << "Run CPLEX..." << endl;
         run_cplex(C, S); //, myfileOut);
         double end_time = timer.elapsed_time(Timer::VIRTUAL);
         cout << "time: " << (end_time - cur_time ) <<"\n";
         //myfileOut << "time: " << (end_time - cur_time ) <<"\n";*/
         //myfileOut.close();
}



