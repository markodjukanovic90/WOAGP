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
int iter = 0; 
int line_number = 0;
int sigma = 0;
int alg = 0; // type of the algorithm (0: MIP; 1: DP)
// time limit for CPLEX (can be changed at run time via the -t comand line parameter)
double time_limit = 900.0;
string outdir="/home/djukanovic/Desktop/projects/rflcsp/MIP-RFLCS-2d/CP-instances/"; // the directory where .out files will be saved

/** preprocessing structures **/
int n_len; int n_of_sequences; int alphabet_size = 0; // lenght of the largest input string
std::vector<std::vector<std::vector<int16_t>>> occurance_positions; ///< structure for reading the number of occurrences of letter <a> of some string i<=m, starting from pl to the end of string
std::vector<std::vector<int16_t>> M;
std::vector<std::vector<int16_t>> M_revs;
int primal_sol = 20;
int primal_sol_time = 0; // time of BS execution

/* two input strings */
int m = 2; // number of strings in input
vector<int> s1;
vector<int> s2; 



// structures for keys of Hash-mapwhich will store binary variables of the 2-LCPS model

struct Point2d // point for structure in 2D-LCPS
{

   public:
         int a,b, lett;
         Point2d(): a(0), b(0), lett(-1) {}; // default constructor...
         Point2d(int _a,int _b, int _lett): a(_a), b(_b), lett(_lett) {};

         bool operator == (const Point2d & p2) const
         {
              return (a == p2.a and b == p2.b);
         }

         friend std::ostream&  operator << (std::ostream & os, const Point2d&  p)
         {
                os<< "(" << p.a << ", " << p.b << ") " << endl;
	          return os;
         }  
};

class Hash2d  //hashing the point
{

      // the hash-key function works for strings up to the lenght of 1000 
      // works only if the length of strings <= 5000
      public:
	    std::size_t operator()(const Point2d & pl ) const
	    {
		return 5000 * pl.a + + pl.b;  
	    }
};


inline int stoi(string &s) {

       return atoi(s.c_str());
}

inline double stof(string &s) {

       return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

    int iarg = 1;

    while (iarg < argc) {
        if (strcmp(argv[iarg],"-alg")==0) alg = atoi(argv[++iarg]);
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

    cout << "ispis" << endl;
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


/** following the notation given in our concept **/
bool ordering_conflict_1(const Point2d& p1, const Point2d& p2)
{
     if(p1.a > p2.a and  p1.b > p2.b and p1.lett < p2.lett) //vars in conflict
	return true;
     else
        return false;
}

bool ordering_conflict_2(const Point2d& p1, const Point2d& p2)
{   
     if(p1.a < p2.a and  p1.b < p2.b and p1.lett > p2.lett)  //vars in conflict
	return true;
     else
        return false;
}

bool intersection_conflict(const Point2d& p1, const Point2d& p2)
{
     if((p1.a >= p2.a and  p1.b <= p2.b) or (p1.a <= p2.a and  p1.b >= p2.b)) // vars in conflict
	 return true;
     else
         return false;
}


// redukcija broja varijabli & preprocessing structures:

void M_fill() // creating M_ij score matrix (of dimension |s_i| x |s_j| M[x, y] --> LCS(s_0[x, |s_0|], s_1[y, |s_1|])
{
   M = vector<vector<int16_t>>((int)s1.size() + 1, vector<int16_t>((int)s2.size() + 1, 0));
   for (int x =  (int) s1.size(); x >= 0; x--)
   {
 	for (int  y = (int) s2.size() ; y >= 0; y--){
             if(x == (int) s1.size())
 	        M[x][y] = 0;
 	     else
 	     if(y == (int) s2.size())
 	        M[x][y] = 0;
 	     else
             if( s1[ x ]  == s2[ y ]) // x + 1 and y + 1
 		 M[x][y]  = M[x + 1][y + 1] + 1;
             else
 	         M[x][y] = std::max(M[x][y + 1], M[x + 1][y]);
        }
   }
}
// parallel based structure

void M_rev() // M_revs[i][j]: corresponds to the lenght of LCS({s1[1, i], s2[1, j]})
{    cout << s1.size() << " " << s2.size() << endl;
     M_revs = vector<vector<int16_t>>((int)s1.size() + 1, vector<int16_t>((int)s2.size() + 1, 0));
     for (int i = 0; i <= s1.size(); ++i)  
     {  
          for (int j = 0; j <= s2.size(); ++j)  
          {  
               if (i == 0 || j == 0)  
                   M_revs[i][j] = 0;  
               else if (s1[i - 1] == s2[j - 1])  
                   M_revs[i][j] = 1 + M_revs[i - 1][j - 1];  
               else
                   M_revs[i][j] = std::max(M_revs[i - 1][j], M_revs[i][j - 1]);  
          }  
     }  
}

vector<int16_t> occurances_string_letter(std::vector<int>& str, int a)
{
     // occurances of letter <a> in string str: 0 - means does not appear: vector[pl] denoted the number of occurrences of letter a starting from s[pl, end];
     vector<int16_t> occurances_letter(n_len, 0);
     // predecessor of S_i for  letter <letter> :
     int number = 0;  
     for (int i = n_len - 1;  i >= 0; i--)// traverse through X
     {
          if(i < str.size() and str[i] == a)
             number++;
           occurances_letter[i] = number;  // occ[ a ][ 0 ]: means all #of occ. of a in X
     } 
     return occurances_letter;
}

vector<vector<int16_t>> occurances_all_letters(vector<int>& str)
{
 
    vector<vector<int16_t>> occurances_all(alphabet_size, vector<int16_t>(n_len, 0));  
    for (int a = 0; a < alphabet_size; a++) // iterate for all letters from alphabet
         occurances_all[a] = occurances_string_letter(str, a);

    return occurances_all;
}

void structure_occurances()
{
    occurance_positions = vector<vector<vector<int16_t> > >(2, vector<vector<int16_t>> (alphabet_size, vector<int16_t> (n_len, 0)));
    occurance_positions[ 0 ] =  occurances_all_letters(s1);
    occurance_positions[ 1 ] =  occurances_all_letters(s2);
}

// occurance based UB for RFLCS

int ub1_left(Point2d& c1) // UB1(S[(1,1), (left, right)])
{   
    int left = 0;
    for(int a = 0; a <  alphabet_size; a++) {
        if(a != c1.lett){ // repetition-free constraint:
           if ((occurance_positions[ 0 ][ a ][ 0 ] - occurance_positions[ 0 ][ a ][ c1.a ] > 0 )
               and ( occurance_positions[ 1 ][ a ][ 0 ] - occurance_positions[ 1 ][ a ][ c1.b ] > 0 )){ // presented in both left part of strings
               left++;
           }
        }
    }
    return left;
}

int ub1_right(Point2d& c1) // UB1(S[(1,1), (left, right)])
{
    
    int right = 0;
    for(int a = 0; a <  alphabet_size; a++) {
        if(a != c1.lett ) { // repetition-free constraint:
           if (occurance_positions[ 0 ][ a ][ c1.a ] > 0 and  occurance_positions[ 1 ][ a ][ c1.b ] > 0) // letter <a> presented in both strings right parts of s_1 and s_2
               right++;
       }
    }
    return right;
}

//UB_2 bound for RFLCS:
int ub2_left(Point2d& c1)
{
    int left = M[ 0  ][ 0 ] - M[ c1.a  ][ c1.b ]; // approximation
    return left;
}

int ub2_left_rev(Point2d& c1) // using revers scoring matrix M_rev
{
   int left =  M_revs[ c1.a ][ c1.b ]; // approximation
   return left;
}

int ub2_right(Point2d& c1)
{
    int right =  M[ c1.a + 1  ][ c1.b + 1 ]; // exact value;
    return right;
}
/* reduction 2 can be seen as a kind of tightenning the reduction1 rule: */
bool reduction1(Point2d& c1){
 
     int ub_min_right = std::min(ub1_right(c1), ub2_right(c1));
     int ub_min_left  = std::min(ub1_left(c1), ub2_left(c1));

     return (ub_min_left + ub_min_right + 1 < primal_sol); // if true: remove node corresp. to c1
}
// using parallel scoring matrix M_rev
bool reduction2(Point2d& c1){
 
     int ub_min_right = std::min(ub1_right(c1), ub2_right(c1));
     ub_min_right = std::min(ub_min_right, s1.size() -  );
     int ub_min_left  = std::min(ub1_left(c1), ub2_left_rev(c1));
     // cout << ub_min_right<< " " << ub_min_left << endl;
     return (ub_min_left + ub_min_right + 1 < primal_sol); // if true: remove node corresp. to c1
}


bool reduction(Point2d& c1){
      
     int l_c1 = std::min(c1.a, c1.b);
     int r_c2 = std::min( s1.size() - c1.a - 1, s2.size() - c1.b - 1);
     return ( l_c1 + r_c2 + 1 < primal_sol);  	
}

/* end of preprocessing structures  */




void run_cplex_lcisp(vector<int>& s1,vector<int>& s2, Timer& timer){ //, ofstream& myfileOut) {

  IloEnv env;
  env.setOut(env.getNullStream());
  try
  {
   // define the LCIS problem model...
   IloModel model(env);
   int vars_num = 0;
   IloObjective obj = IloMaximize(env);
   // defining the set of binary variables Z
   unordered_map<Point2d, IloNumVar, Hash2d> Z;
   for(int i = 0; i < s1.size(); ++i){
       for(int j = 0; j < s2.size(); ++j){
          if(s1[i] == s2[j] ){ // 2d-matching...
             
             bool reducing_ij = false;
             Point2d p(i, j, s1[i]);
             //if (reduction2(p)) // with the parallel structure
             //    reducing_ij = true;   
             if( !reducing_ij )
             {   IloNumVar myIntVar(env, 0, 1, ILOINT);
                 Z[p] = myIntVar;
                 obj.setLinearCoef(Z[p], 1.0); // define objective function
                 ++vars_num;
             }
          }
       }
   }
   cout << "#Vars: " << vars_num << endl;
   map<int, IloExpr> cs;  
   // repetition-free constraints:
   cout << "Add repetition-free--constraint: " << endl;
   for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator itx = Z.begin(); itx != Z.end(); ++itx){
      
       if(cs.find(s1[(*itx).first.a]) == cs.end() ){
          IloExpr xpr(env);
          cs.insert({s1[(*itx).first.a], xpr });
       }
       cs[s1[(*itx).first.a]] += (*itx).second;
   }
   https://bitbucket.org/dashboard/overview
   int constraints_num = 0;
   for(int i = 0; i < alphabet_size; i++){
       if(cs.find(i) != cs.end())
          model.add(cs[i] <= 1); 
       constraints_num++; 
   }

   std::cout << "Set of varables Z is generated. |Z| = " << Z.size() << std::endl;
   // define a set of constraints:
   for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator it1 = Z.begin(); it1 != Z.end(); ++it1){
       for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator it2 = Z.begin(); it2 != Z.end(); ++it2){
          if( !( (*it1).first == (*it2).first ) and 
               ( intersection_conflict( (*it1).first, (*it2).first) 
                 or ordering_conflict_1((*it1).first, (*it2).first)    
                 or ordering_conflict_2((*it1).first, (*it2).first) 
              ))
          {   // conflict constraints
               model.add( (*it1).second + (*it2).second <= 1);
               constraints_num++;
          }
       }
   }

   std::cout << "#Constraints: "<< constraints_num << std::endl;
   // run CPLEX set up parameters
   vector<double> times;
   vector<double> results;
   vector<double> gaps;
   //solve the model IloLinearNumExpr objective = cplex.linearNumExpr();
   model.add(obj);
   IloCplex cplex(model);
   IloExpr expression(env);
   for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator itx = Z.begin(); itx != Z.end(); ++itx){
       expression += (*itx).second;
   }
   // additional constriants: conditional constrints:
   for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator it = Z.begin(); it != Z.end(); ++it){
       int a = (*it).first.a; 
       int b = (*it).first.b; 
       Point2d pd = (*it).first;
       int lett =  (*it).first.lett; // ovo pogledati...
       int right2 = std::min( (int) s1.size() -  a - 1, (int) s2.size() - b - 1 );
       int right = std::min( alphabet_size - s1[a] - 1, right2 ); 
       int left =  std::min( lett - 1, std::min(a, b) );  // cout << left << " " << right << " " << right + left + 1.0 << endl;
       //strengthened:
       int ub_min_right = std::min(ub1_right( pd ), ub2_right( pd ));
       int ub_min_left  = std::min(ub1_left( pd ), ub2_left( pd ));
       int right_s = min(right, ub_min_right);
       int left_s  = min(left, ub_min_left);
       //IloConstraint ir(env, 0.0, expression, right + left + 1.0);
       model.add( IloIfThen(env, (*it).second > 0, expression <= right_s + left_s + 1.0 ));
   }
    cout << "run..." << endl;
   //int time_limit = 900;
   // pass the time limit to CPLEX
   cplex.setParam(IloCplex::TiLim, time_limit);
   // the following two parameters should always be set in the way as shown
   cplex.setParam(IloCplex::NodeFileInd, 2);
   cplex.setParam(IloCplex::Threads, 1);
   IloNum lastObjVal = std::numeric_limits<double>::max();
   ofstream fout("mylog.log");
   // tell CPLEX to make use of the function 'loggingCallback' for writing out information to the screen
   //cplex.use(loggingCallback(env, timer, times, results, gaps, iter, lastObjVal));
   cout <<"Solve..." << endl;
   // cplex set up:
   cplex.setParam(IloCplex::EpGap,0.000000001); 
   cplex.setParam(IloCplex::MIPDisplay, 2);    // display every MIPInterval nodes
   cplex.setParam(IloCplex::WorkMem,700);    // set working memory to 1.5GB 
   cplex.setParam(IloCplex::MIPInterval, 3);
   //cplex.setOut(fout);
   cplex.solve();
   
   if (cplex.getStatus() == IloAlgorithm::Optimal or cplex.getStatus() == IloAlgorithm::Feasible)
   {
       if(cplex.getStatus() == IloAlgorithm::Optimal)
          cout << "CPLEX finds optimal" << endl;
       else
          cout << "CPLEX finds feasible solution" << endl;

       double lastVal = double(cplex.getObjValue());
       // print the objective point
       cout << "nodes/vertices in the solution: {" <<endl;
       bool first = true;
       for(unordered_map<Point2d, IloNumVar, Hash2d >::iterator it = Z.begin(); it != (Z.end()); ++it){
           IloNum xval = cplex.getValue((*it).second);
           // the reason for 'xval > 0.9' instead of 'xval == 1.0' will be explained in class
           if (xval > 0.9) {
               cout << (*it).first;
               //myfileOut << (*it).first;
           }
       }
       cout << "}" << endl;
       cout << "value: " << lastVal << endl;
       //myfileOut << "value: " << lastVal << endl;
   }
}
catch(IloException& e) {
        cerr  << " ERROR: " << e << endl;
}
env.end();
}

// DP method in O(mn) complexity
int DP_LCIS(const vector<int>& arr1, int n, const vector<int>& arr2,  int m,int* table, int& ind) 
{ 
    // table[j] is going to store length of LCIS 
    // ending with arr2[j]. We initialize it as 0, 
    for (int j=0; j<m; j++) 
        table[j] = 0; 
  
    // Traverse all elements of arr1[] 
    for (int i=0; i<n; i++) 
    { 
        // Initialize current length of LCIS 
        int current = 0; 
  
        // For each element of arr1[], traverse all 
        // elements of arr2[]. 
        for (int j=0; j<m; j++) 
        { 
            // If both the array have same elements. 
            // Note that we don't break the loop here. 
            if (arr1[i] == arr2[j]) 
                if (current + 1 > table[j]) 
                    table[j] = current + 1; 
  
            /* Now seek for previous smaller common 
               element for current element of arr1 */
            if (arr1[i] > arr2[j]) 
                if (table[j] > current) 
                    current = table[j]; 
        } 
    } 
    // The maximum value in table[] is out result 
    int result = 0; 
    for (int i=0; i<m; i++) 

        if (table[i] > result) {
           result = table[i]; 
  	   ind = i;
        }

    return result; 
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
        iter = line_number = 1;
        ifstream input(inputFile[0].c_str());
        //vector<int> str1; vector<int> str2; // input strings
        int m = 2;
        string line; int n1, n2;
        // opening the corresponding input file and reading the problem data
        while (std::getline(input, line)) {
               std::cout << line << "\n";
               std::string buf;
               std::stringstream ss(line.c_str());
               if(line_number == 1){
                  int c = 0;
                  while(ss >> buf){
                        if(c == 0)
                           m = std::atoi(buf.c_str());
                        else
                           alphabet_size = std::atoi(buf.c_str());
                         c++;
                  }
               }
               else 
               if(line_number == 3) //first string 
               {   
                   while(ss >> buf)
                         s1.push_back(atoi(buf.c_str()));
               }
               else
               if(line_number == 5){ // second string
                  while(ss >> buf)
                        s2.push_back(atoi(buf.c_str()));
               }
               line_number++;
         }
         input.close();
         n_len = s1.size();
         cout << "End of input " << inputFile[0].c_str() << ".out" << endl;
         /*cout << "Input data: " << "sigma: " << sigma << endl;
         print_vector(str1);
         cout << " str2: "<< endl;
         print_vector(str2);
         */
         Timer timer;
         double cur_time = timer.elapsed_time(Timer::VIRTUAL);
         // write the results in the output...
         std::stringstream filenameOut;
         filenameOut << outdir << inputFile[0].c_str() << ".out";
         std::string filenameO = filenameOut.str();
         ofstream myfileOut (filenameO);

         if(alg == 0) 
         {
            structure_occurances();
            M_fill();
            M_rev();
            
            cout << "Run Cplex..." << endl;
            run_cplex_lcisp(s1, s2, timer); //, myfileOut);

         }else{

            cout << "Run DP recursion to solve the LCIS problem: " << endl;
  	    int table[ (int)s2.size() ];
  	    int ind = 0;
            int res = DP_LCIS(s1, (int)s1.size(), s2 , (int)s2.size(), table, ind); 
            cout << "value: " << res << endl;
         }
         double end_time = timer.elapsed_time(Timer::VIRTUAL);
         cout << "time: " << (end_time - cur_time) <<"\n";
         //myfileOut << "time: " << (end_time - cur_time ) <<"\n";
         //myfileOut.close();
}

