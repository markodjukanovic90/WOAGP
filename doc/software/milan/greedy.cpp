#include <list>
// #include <minisat/core/Solver.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <bits/stdc++.h>
#include <utility>

using namespace std;
using namespace std::chrono;

typedef pair<float, float> Point_2;


/** pomocne strukture za gridi metod **/

vector<float> Surface; // for a vertex indexed by i, we return the surface of V(i)
vector<float> Cost; // at index i --> w_i
vector<vector<float>>Intersection; // Intersection[i][j] = Volume of V(i) intersected by V(j)
vector<set<Point_2>> S; // i-th position of S is a set of all points that are visible from vertex i
#define INFEASIBLE 1000000
int cardinalityD; // |D(P)|
int wTotal = 0; // total weight 
int t_lim = 0; // time limit
int greedy = 0; // tip gridija => 0: Lovasz; 1: price-per-unit; 2: intersection-based; 3: Greedy by Dragan
std::string path;
int n = 0; // |V(P)|
int w_type = 0; // ti ptezine koju pozivamo; 0: tezina proporcionalna velicini vidljivosti svakog vrha; 1: --; 2: --
/** kraj pomocnih str. za gridi metode **/


inline int stoi(string &s) {

     return atoi(s.c_str());
}

inline double stof(string &s) {

     return atof(s.c_str());
}

void read_parameters(int argc, char **argv) {

     int iarg=1;
     if (strcmp(argv[iarg],"-f") == 0) path = (argv[++iarg]);
     else if(strcmp(argv[iarg],"-t") == 0) t_lim = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-greedy") == 0) greedy = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-w_type") == 0) w_type = atoi(argv[++iarg]);
     else ++iarg;
}

void read_from_file(std:: string path)
{
     // TODO: read from file
}

// unija skupova iz @S sa indeksima iz @index
/**
param:
@index: skup indeksa cije vrijednosti treba sojiti uniom
@S: Skuovi tacaka vidljivosti za odgovarajuce strazare
**/
set<Point_2> union_by_index(vector<int> &index)
{
  set<Point_2> result;
  for (size_t i = 0; i < index.size(); i++)
    for(Point_2 p : S[i])
       result.insert(p);
  
  return result;
}

// kardinalnost skupa tacaka koje su pokrivene sa tjemenima ciji su indeksi dati sa @indeks
/**
param:
@index: skup indeksa cije vrijednosti treba spojiti uniom
@S: Skuovi tacaka vidljivosti za odgovarajuce strazare
**/
int cardinality_by_index(vector<int> &index)
{
  return union_by_index(index).size();
}

int f(vector<set<Point_2>> &C)
{
  std::set<Point_2> d;
       for(set<Point_2> di : C){
     
           for (Point_2 dii: di){
               //cout << "dii" << dii << endl;
               d.insert(dii);
           }
       }
       return (int)d.size();
}

int f_minus(vector<set<Point_2>>& C, set<Point_2>& ss) // f(C u {s}) - f(S)
{
    int count = 0; // ako i \in ss nije niti u jednom od skupova u S, onda count++;
    if(C.empty())
       return ss.size();

    for(Point_2 i: ss) 
    {
        bool presented = false;
        for(set<Point_2>& di : C){
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

float greedy_criterion_1(int i)
{
    return Cost[i] / Surface[i];

}

/**   Gridi funkcija (6)  prema indeksiranju u radu
Param:
@i: indeks cvora za koji se racuna gridi vrijednost 
@indeks: trenutni parcijalni skup  --> samo indeksi skupova uzeti radi optimizacije
@S: skup svih tjemena
 **/

float greedy_criterion_2(vector<int>& indeks, int i) 
{
      if(std::find(indeks.begin(), indeks.end(), i) != indeks.end() ) // i vec u skupu indeks koji odg. parcijalnom rjesenju...
         return INFEASIBLE;
 
      float num = 0.0;
      int den = S.size() - indeks.size();
      for(int j = 0; j < S.size(); ++j) {
          if(S[j] != S[i] &&  std::find(indeks.begin(), indeks.end(), j) == indeks.end())
          { // S_j ne smije biti u vec dodanom parcijalnom skupu 
              num += Intersection[i][j];
          } 
      }
      return num / den;
}
/**
param:
@S: skup svih skupova (instance) 
@C: parcijalno rjesenje (trenutno) 
@i: index skupa S_i koji se razmatra za dodavanje u parcijalni skup @C
**/
float greedy_criterion(int i, vector<set<Point_2>>& C) // take s_i from S
{     
    set<Point_2> s = S[i];
    int f_m = f_minus(C, s);  //cout << "f_m: " << f_m << endl;
    if(f_m == 0)
       return INFEASIBLE;

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
/**
param: 
@indeks: indeks svih skupova koji su u trenutnom parcijalnom rjesenju;
@i: ineks skupa za koji racunamo gridi vrijednost.

**/
float gridi_criterion_dragan(vector<int> &index, int i)
{
      if(std::find(index.begin(), index.end(), i) != index.end())
          return INFEASIBLE;
      index.push_back(i); // trebalo bi ovo optimizovati 

      int correct_total = cardinality_by_index(index);  // broj pokrivenih tacaka diskretizacije
      int incorrect_total =  cardinalityD - correct_total;
      
      int wPartial = 0; 
      for(int i: index)
      {
          wPartial += Cost[ i ]; 
      }

      double obj1part = ((double) wPartial) / wTotal;
      double obj2part = (float) incorrect_total / cardinalityD;
      double obj = obj1part + obj2part;
      
      // drop vertex i from indeks 
      vector<int>::iterator it = std::find(index.begin(), index.end(), i); 
      index.erase(it);

      return obj;
}


int min_greedy(vector<set<Point_2>>& C, vector<int>& indeks)
{         //cout << "min_greedy" << endl;
          float g_m = -1; int dodaj;
          for(int i = 0; i < S.size(); ++i)
          { 
              float g_mi;
              switch(greedy){
                  case 0: g_mi = greedy_criterion(i, C); break;
                  case 1: g_mi = greedy_criterion_1(i); break;
                  case 2: g_mi = greedy_criterion_2(indeks, i); break;
                  default: g_mi = gridi_criterion_dragan(indeks, i);
              }
              if(g_mi <= g_m and g_mi != INFEASIBLE and !findA(indeks, i)) 
              { 
                 dodaj = i;
                 g_m = g_mi;   
              }
          }//cout << "dodati....." << dodaj << endl;
          return dodaj;
}

float greedy_procedure()
{
     vector<set<Point_2>> C;
     vector<int> indeks;

     int f_S = f(S); //cout << "f_S: " << f_S << endl;
     int f_C = 0;

     while(f_C != f_S) 
     { 
           int index_set = min_greedy(C, indeks);
           //cout<<pol[index_set]<<endl;
           //cout<<index_set<<endl;
           if(!findA(indeks, index_set)){ //jos nije dodan
           
               C.push_back(S[index_set]);//cout << "dodaj ----> " << index_set << endl;
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

int main( int argc, char **argv ) {

    read_parameters(argc, argv);
    read_from_file(path); // fill Cost, Intersection and Surface

    // ----------------------- dodjela tezina (proporcionalno broju pokrivenih tacaka iz D(P) ------------------------------------
    switch(w_type){
 
          case 0: { // tezina proporcionalna sa velicinom skupa S[i]
                  for(auto& X: S){
                      cardinalityD += X.size();
                  }
                  for (size_t i = 0; i < n; i++){
                      float w_i = n * n * ( ((float) S[i].size()) / cardinalityD );
                      wTotal += w_i;
                      Cost.push_back(w_i);
                  }; break;}
          default: {
                    for(size_t i = 0; i < n; i++) // non-weighted version of the problem
                       Cost.push_back(1);
                   }
   }
    // ---------------------------------greedy------------------------------------
    auto start = high_resolution_clock::now();
    float s = greedy_procedure();
    auto stop = high_resolution_clock::now();
    // ------------------------------greedy - end------------------------------------

    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution: " << duration.count() << " microseconds" << endl;
    cout << "Dovoljan broj cuvara je: " << s << endl;
    
    return EXIT_SUCCESS;
}
