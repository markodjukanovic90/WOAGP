#include <list>
// #include <minisat/core/Solver.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <sstream>
#include <chrono>
#include <bits/stdc++.h>
#include <utility>
#include <string>
#include <math.h>       /* sqrt */
#include <unordered_map> 
#include <unordered_set> 
#include<string>

// the following "include" is necessary for the correct working/compilation of CPLEX. You should adapt this path to your installation of CPLEX
#include "/home/marko/Desktop/CPLEX_Studio127/cplex/include/ilcplex/ilocplex.h"

using namespace std;
using namespace std::chrono;
ILOSTLBEGIN

/** pomocne strukture za CMSA metod **/
vector<set<int>> Sets;

#define INFEASIBLE 1000000
 
int t_lim = 0; // time limit 
int t_solver = 100;
std::string path ;
std::string out_path = "/home/marko/Desktop/kod/";
int alg = 0; // 0: CPLEX, 1: CMSA
int n = 0;  
int m = 0;  

/** another try to spped up operations over CoveredPoints structure **/
vector<bool>covered;
  

/** end of the structure **/
 
inline int stoi(string &s) {

     return atoi(s.c_str());
}

inline double stof(string &s) {

     return atof(s.c_str());
}

//CPLEX model: returns a solution

vector<int> run_cplex(vector<int>&  V, int t_solver = 0)
{

  vector<int> solution; // save solution here 
  IloEnv env;  env.setOut(env.getNullStream());
  cout << "run_cplex " << n << " " << m << endl;
  auto startTime = high_resolution_clock::now(); // time measurement
  
  try
  {
   IloModel model(env);
   IloObjective obj = IloMinimize(env); 
   
   // defining the set of binary variables Z
   vector<IloNumVar> X;  
   for(int i = 0; i < m; ++i){
       //cout << "cost i " << Cost[i] << endl;
       IloNumVar myIntVar(env, 0, 1, ILOINT);
       X.push_back(myIntVar); // x_i  \in {0, 1}
       obj.setLinearCoef(X[i], 1.0); // sum_i c_i x_i
   }
   // constraints 
   cout << "#Vars: " << X.size() << endl;
   
   // constraints 
   int index = 0;
   
   for(int i = 0; i < n; ++i) 
   {
       IloExpr expr_i(env); 
       bool add = false;
       for(int j = 0; j < m; ++j)   
       {
           if( Sets[j].count( i+1) ){ // item i se nalazi u skupu j
               expr_i += X[j];
               add = true;
           }
       }
       if(add)
       {   
          model.add(expr_i >= 1);      
          index++ ;
           
       } 
   }  
 
   // add constraints to solve SUBINSTNACE in CMSA
   if(V.size() > 0)
   {
       cout << "Add subproblem constraints " << endl;
       for(int i = 0; i < m; ++i)
       {
            if(count(V.begin(), V.end(), i ) == 0)
               model.add(X[i] == 0);
       }  
    
   }  
   cout <<"Number of constr " << index << endl;
   //solve the modelIloLinearNumExpr objective = cplex.linearNumExpr();
   model.add(obj);
   // add model to CPLEX
   IloCplex cplex(model);
   
   //int time_limit = 900;
   // pass the time limit to CPLEX
   if(t_solver == 0)
      cplex.setParam(IloCplex::TiLim, t_lim);
   else
      cplex.setParam(IloCplex::TiLim, t_solver);
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
         
       for(int i = 0; i < X.size(); ++i){
           IloNum xval = cplex.getValue(X[i]);
   
           if (xval > 0.9) {
               cout << "S_" << ( i + 1 ) << ", ";
               solution.push_back( i );
               //myfileOut << (*it).first;
           }
       }
       cout << "}" << endl;
       
       
       if(alg == 0) // if a CPLEX is called 
       {// write the results in a file 
       	string outname = "" ;
       	outname = out_path + "scp_" + std::to_string(n) + "_" + std::to_string(m) + ".txt.out";
       	std::ofstream outfile; // stream
       
       	outfile.open(outname);
       	auto stopTime = high_resolution_clock::now();
       	auto duration = duration_cast<std::chrono::duration<float>>(stopTime - startTime); 
       	outfile << "obj: " << lastVal << " time: " << duration.count() ;
       	
       
      }   
      return solution; //lastVal;
   }
   else{
       cout << "Nema rjesenja" << endl;
       return solution;//-1;
   }  
}
catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
}
env.end();
}

// end CPLEX model

vector<string> split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);
        if (pos == string::npos) pos = str.length();
        string token = str.substr(prev, pos-prev);
        if (!token.empty()) tokens.push_back(token);
        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return tokens;
}
 
void read_parameters(int argc, char **argv) {

     int iarg=1;
     cout << "Reding params..." << endl;
     while (iarg < argc) {
        if (strcmp(argv[iarg],"-f") == 0)          path = (argv[++iarg]);
        else if(strcmp(argv[iarg],"-t") == 0)      t_lim = atoi(argv[++iarg]);
        else if(strcmp(argv[iarg],"-t_solve") == 0)      t_solver = atoi(argv[++iarg]);
        else if(strcmp(argv[iarg], "-alg") == 0)  alg = atoi(argv[++iarg]);
        else ++iarg;
     }
}

void read_from_file(std:: string path)
{
     // TODO: read from file
          cout << "Reding from file..." << endl;
     const string delim = " ";
     std::ifstream infile(path);
     std::string line;
     std::getline(infile, line);
     
     vector<string> nm = split(line, delim);
     n = stoi(nm[0]); // cardinality of X
     m = stoi(nm[1]); // m -- number of sets 
     cout << n << " " << m << endl;
     int temp = 0;
     
     while(temp < m){
        std::getline(infile, line);
        vector<string> S_i = split(line, delim);
        set<int> Si;
        for(string& str : S_i )
            Si.insert(stoi(str));
        
        Sets.push_back(Si);  
        temp++;
    }
}

int f(set<int> &  Cover, int index , vector<set<int>>&  Sets)
{   
     std::set<int> result = Cover;
     result.insert(Sets[index].begin(), Sets[index].end());
     return result.size() - Cover.size();

}   


bool sortFunction(pair<int, int>&  p1, pair<int, int>&  p2)
{
     return(p1.second > p2.second);
   
} 
// returns index with best greedy value:
int greedy_iteration(set<int> &  Cover, float d_rate = 0, float prob = 1)
{
     
     int max = 0;
     int index = -1;
     
     vector<pair<int, int>> extensions; // need to be sort out 
     
     for(int i = 0; i < m; ++i)
     {
        
         if(prob == 1)
         {
            int fi = f(Cover, i, Sets);
        
            if(fi  > max)
            {
               max = fi;
               index = i;
            } 
         }else { // introduce randomness 
            
             int fi = f(Cover, i, Sets);
             if( fi > 0 ) // relevant extension
                extensions.push_back(std::make_pair( i, fi ));
         }    
       }     
        // sort extensions:
       std::sort(extensions.begin(), extensions.end(), sortFunction);
       // random index in extensions:
       
       srand( time(NULL) ); // reset SEED val
       float prGen = rand();
       if(prGen <= prob) // choose best among extensions
          index = extensions[0].first; 
       else{
       
          int range = int( extensions.size() * d_rate) - 0  ;
          srand( time(NULL) ); // reset SEED val
          index = rand() % range + 0;
          //cout  << extensions[ 0 ].second << "----" << extensions[ 1 ].second  << endl;
          //cout << "rand " << index << " " << extensions.size() << endl;
          index =   extensions[ index ].first; 
       }      
       return index;
} 


vector<int> greedy_procedure(float drate, float prob)
{   
    //TODO best cover among others 
    set<int> Cover;
    vector<int> partial;   
    while( Cover.size() < n )
    {   
        
        int index = greedy_iteration( Cover, drate, prob);
        partial.push_back( index ); // cout << index << endl;
        // update Cover
       Cover.insert(Sets[index].begin(), Sets[index].end()); 
    }
    //cout << "Objective: " << Cover.size() << endl;
    return partial; //Cover.size();
} 
// RANDOM GREEDY PROCEDURE:

vector<int> greedy_procedure_random(double d_rate, double prob)
{   
    //TODO best cover among others 
    set<int> Cover;
    vector<int> partial; // cout << "n " << n << endl;
    while( Cover.size() < n )
    {   
        
        int index = greedy_iteration(Cover, d_rate, prob );
        partial.push_back( index ); // cout << index << endl;
        // update Cover
       Cover.insert(Sets[index].begin(), Sets[index].end()); 
    }
    //cout << "Objective: " << Cover.size() << endl;
    return partial; //Cover.size();
} 

void CMSA(int agemax, int na, float drate, float prob, int t_solver = 0)
{
    cout << " RUN CMSA " << endl;
    float sbest = 10000000;
    set<int> Vprime;
    vector<int> age(m, 0);
    
    auto start = high_resolution_clock::now();
    auto stop = high_resolution_clock::now();
    
    auto duration = duration_cast<std::chrono::duration<float>>(stop - start); 
    while(duration.count() < t_lim)
    {
    
          for(int i=1; i <= na; ++i)
          {
              vector<int> S(greedy_procedure(drate, prob));
              for(int v : S)
              {
                  if( Vprime.count(v) == 0 )
                  {      
                     age[v] = 0;
                     Vprime.insert( v );
                  }
              }         
          }
          // Apply Exact solver:
          vector<int> partial(Vprime.begin(), Vprime.end());
          
          vector<int> SprimeOpt(run_cplex(partial, t_solver)); // solve subinstance via CPLEX
          
          if(SprimeOpt.size() < sbest)
             sbest = SprimeOpt.size();
          // ADAPT mechanism depends on: Vprime, SprimeOpt, agemax
          for(int v: Vprime)
              age[v]++;
          
          for(int v: SprimeOpt)
              age[v] = 0;
          // discard sol. components > agemax
          int index = 0;
          for(int v: age)
          { 
              if( v > agemax ) // remove index from Vprime
              {    auto it = Vprime.find ( index);
                   Vprime.erase (it, Vprime.end());
              }
              index++;
          }
          auto stop = high_resolution_clock::now();
          duration = duration_cast<std::chrono::duration<float>>(stop - start); 
    }
    
    // write the results in a file 
    string outname = "" ;
    outname = out_path + "scp_" + std::to_string(n) + "_" + std::to_string(m) + ".txt.out";
    std::ofstream outfile; // stream
    outfile.open(outname);
    outfile << "obj: " << sbest << " time: " << duration.count() ;

}
void write_test(string output, string tekst)
{
  std::ofstream outfile;

  outfile.open(output, std::ios_base::app); // append instead of overwrite
  outfile << tekst;
}


int main( int argc, char **argv ) 
{

    read_parameters(argc, argv);
    read_from_file(path);  
    // instace data print
    
    for(set<int> & set : Sets )
    {
       for(std::set<int>::iterator it=set.begin(); it!=set.end(); ++it)
           std::cout << *it << " ";
        std::cout << std::endl;
    }
    vector<int>  Vx;
    
    cout << "alg: " << alg << endl;
    if(alg == 0)
       run_cplex(Vx);
    else
       CMSA(2, 10, 0.3, 0.5, t_solver);
   /* 
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::duration<float>>(stop - start); 
   cout << "Time Execution: " << duration.count() << "seconds" << endl;
   
 */   
    return EXIT_SUCCESS;
}
