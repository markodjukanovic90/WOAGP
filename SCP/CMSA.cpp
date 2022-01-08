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

// the following "include" is necessary for the correct working/compilation of CPLEX. You should adapt this path to your installation of CPLEX
#include "/home/marko/Desktop/CPLEX_Studio127/cplex/include/ilcplex/ilocplex.h"

using namespace std;
using namespace std::chrono;
ILOSTLBEGIN

/** pomocne strukture za CMSA metod **/
vector<set<int>> Sets;
int alg = 0; // if alg = 1: we execute a Greedy+LS procedure

#define INFEASIBLE 1000000
 
int t_lim = 0; // time limit intersection-based;
std::string path ;
std::string output = "";

int n = 0;  
int m = 0;  

/** another try to spped up operations over CoveredPoints structure **/
vector<bool>covered;
 /**
struct hashFunc{
    size_t operator()(const Point_2 &k) const{
    size_t h1 = std::hash<double>()(k.first);
    size_t h2 = std::hash<double>()(k.second);
    return (h1 ^ (h2 << 1));
    }
};

struct equalsFunc{
  bool operator()( const Point_2& lhs, const Point_2& rhs ) const{
    return (lhs.first == rhs.first) && (lhs.second == rhs.second);
  }
};

unordered_map<Point_2, bool, hashFunc, equalsFunc> pointDPMapping;

**/

/** end of the structure **/
 
inline int stoi(string &s) {

     return atoi(s.c_str());
}

inline double stof(string &s) {

     return atof(s.c_str());
}

 

// begin CPLEX model

float run_cplex(vector<set<int>>&  Sets, int m, int n,  vector<int>&  V)
{

  IloEnv env;  // cout << "run_cplex " << endl;
  env.setOut(env.getNullStream());
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
       IloExpr expr_i(env); bool add = false;
       for(int j = 0; j < m; ++j)   
       {
           if( Sets[j].count( i ) ){ // item i se nalazi u skupu j
               expr_i += X[j];
               add = true;
           }
       }
       if(add)
          model.add(expr_i >= 1);      
       
       index++ ;
   }  
   // add constraints to solve SUBINSTNACE in CMSA
   if(V.size() > 0)
   {
       for(int i = 0; i < m; ++i)
       {
            if(count(V.begin(), V.end(), i ) == 0)
               model.add(X[i] == 0);
       }  
    
   }  
   
   //solve the modelIloLinearNumExpr objective = cplex.linearNumExpr();
   model.add(obj);
   // add model to CPLEX
   IloCplex cplex(model);
   
   //int time_limit = 900;
   // pass the time limit to CPLEX
   cplex.setParam(IloCplex::TiLim, t_lim);
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
               //myfileOut << (*it).first;
           }

       }
       cout << "}" << endl;
       //myfileOut << "value: " << lastVal << endl;
       return lastVal;
   }
   else{
       cout << "Nema rjesenja" << endl;
       return -1;
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
     if (strcmp(argv[iarg],"-f") == 0) path = (argv[++iarg]);
     else if(strcmp(argv[iarg],"-t") == 0) t_lim = atoi(argv[++iarg]);
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
    vector<int> partial; cout << "n " << n << endl;
    while( Cover.size() < n )
    {   
        
        int index = greedy_iteration( Cover, drate, prob);
        partial.push_back( index ); // cout << index << endl;
        // update Cover
       Cover.insert(Sets[index].begin(), Sets[index].end()); 
    }
    cout << "Objective: " << Cover.size() << endl;
    return partial; //Cover.size();
} 
// RANDOM GREEDY PROCEDURE:

vector<int> greedy_procedure_random(double d_rate, double prob)
{   
    //TODO best cover among others 
    set<int> Cover;
    vector<int> partial; cout << "n " << n << endl;
    while( Cover.size() < n )
    {   
        
        int index = greedy_iteration(Cover, d_rate, prob );
        partial.push_back( index ); // cout << index << endl;
        // update Cover
       Cover.insert(Sets[index].begin(), Sets[index].end()); 
    }
    cout << "Objective: " << Cover.size() << endl;
    return partial; //Cover.size();
} 

void CMSA(int na, float drate, float prob)
{
    //TODO:
    float s = 10000000;
    vector<int> partial(greedy_procedure(drate, prob));
 
    for(int i: partial)
       cout << "  " << i << endl;
    run_cplex( Sets, m, n, partial); // solve subinstance via CPLEX
    //

}
void write_test(string tekst)
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

    
    auto start = high_resolution_clock::now();
    
    CMSA(10, 1.0, 0.9);
    
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::duration<float>>(stop - start); 
   cout << "Time Execution: " << duration.count() << "seconds" << endl;
    
    // if(output.compare("") != 0){
      //   string name_polygon = split(split(path, "/")[4], "_")[0];
       //  write_test("Solution");
    // }
    return EXIT_SUCCESS;
}
