#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_polygon_visibility_2.h>
#include <CGAL/Triangular_expansion_visibility_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/General_polygon_set_2.h>
#include <CGAL/Gps_segment_traits_2.h>
#include <CGAL/centroid.h>
#include <list>
#include <minisat/core/Solver.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <bits/stdc++.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
typedef CGAL::Polygon_2<K>                                              Polygon;
typedef CGAL::Point_2<K>                                                Point;
typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef CGAL::Polygon_2<Kernel>                                         Polygon_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Edge_const_iterator                              Edge_const_iterator;
typedef Arrangement_2::Halfedge_const_handle                            Halfedge_const_handle;
typedef Arrangement_2::Halfedge_handle                                  Halfedge_handle;
typedef Arrangement_2::Halfedge_handle                                  Halfedge_handle;
typedef CGAL::Simple_polygon_visibility_2<Arrangement_2,CGAL::Tag_true> RSPV;
typedef CGAL::Polygon_with_holes_2<Kernel>                              Polygon_with_holes_2;
typedef std::list<Polygon_with_holes_2>                                 Pwh_list_2;
typedef CGAL::Gps_segment_traits_2<Kernel>                              Traits_22;
typedef CGAL::General_polygon_set_2<Traits_22>                          Polygon_set_2;
typedef Arrangement_2::Face_handle                                      Face_handle; 
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>          TEV;

CGAL::Cartesian_converter<Kernel,K> converter; 
CGAL::Cartesian_converter<K, Kernel> converter2; 

using namespace std;
using namespace std::chrono;

/** pomocne strukture za gridi metod **/
vector<float> Surface;
vector<int> Cost;
vector<vector<float>>Intersection;
#define INFEASIBLE 100000

/** kraj pomocnih str. za gridi metode **/

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

Point_2* read_point_of_polygon(string& filename, int *n)
{
  Point_2* points;
  *n = 0;

  string text;
  ifstream MyReadFile(filename);

  while(getline (MyReadFile, text))
  {
    vector<string> line = split(text, " ");
    stringstream geek(line[0]); 
    *n = 0; 
    geek >> *n;
    points = new Point_2[*n];
    
    int numerator = 0, denominator = 1;
    vector<string> number;
    for (int i = 1; i < line.size(); i+=2)
    {
      number = split(line[i], "/");
      numerator = 0;
      denominator = 1;

      stringstream geekN(number[0]); 
      geekN >> numerator;
      stringstream geekD(number[1]);
      geekD >> denominator;

      double po1 =(double)numerator/denominator;

      number = split(line[i+1], "/");
      numerator = 0;
      denominator = 1;

      stringstream geekN2(number[0]); 
      geekN2 >> numerator;
      stringstream geekD2(number[1]);
      geekD2 >> denominator;

      double po2 =(double)numerator/(double)denominator;

      points[i/2] = Point_2(po1, po2);
    }

  }
  MyReadFile.close(); 

  return points;
}

Polygon_2 read_polygon(string filename, int *n)
{
  Point_2* points = read_point_of_polygon(filename, n);
  return Polygon_2(points, points+*n);
}

mpfr_void discratePVI(vector<Point_2> &D, Polygon_2 &p, set<Point_2> &points)
{
    if(p.is_simple())
      for (size_t i = 0; i < D.size(); i++)
      {
        int bs = p.bounded_side(D[i]);

        if(bs==CGAL::ON_BOUNDED_SIDE || bs==CGAL::ON_BOUNDARY) //  || bs==CGAL::ON_BOUNDARY
          (points).insert(D[i]);
      }
}

Polygon_2 arrangement_to_polygon(Arrangement_2& arr)
{
  Edge_const_iterator eit = arr.edges_begin();
  vector<Point_2> points_2 = {eit->source()->point(), eit->target()->point()};

  vector<Edge_const_iterator> edges_it;
  for(++eit; eit != arr.edges_end(); ++eit)
    edges_it.push_back(eit);

  while(edges_it.size()>0)
  {
    for (size_t i = 0; i < edges_it.size(); i++)
    {
      eit = edges_it[i];

      if(points_2[0]==eit->source()->point())
      {
        points_2.insert(points_2.begin(), eit->target()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[0]==eit->target()->point())
      {
        points_2.insert(points_2.begin(), eit->source()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[points_2.size()-1]==eit->source()->point())
      {
        points_2.push_back(eit->target()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
        continue;
      }
      if(points_2[points_2.size()-1]==eit->target()->point())
      {
        points_2.push_back(eit->source()->point());
        edges_it.erase (edges_it.begin()+i);
        i--;
      }
    }
  }
  

  int n = arr.number_of_edges();
  Point_2* points = new Point_2[n];

  for (size_t i = 0; i < n; i++)
    points[i] = points_2[i];
  
  return Polygon_2(points, points + n);
}

vector<Polygon_2> visibility_polygon(Point_2* points, int& n)
{
  vector<Polygon_2> vp;
  
  std::vector<Segment_2> segments;
  for (size_t i = 0; i < n-1; i++)
    segments.push_back(Segment_2(points[i], points[i+1]));
  segments.push_back(Segment_2(points[n-1], points[0]));

  Arrangement_2 env;
  CGAL::insert_non_intersecting_curves(env,segments.begin(),segments.end());
  Arrangement_2 regular_output;

  for (size_t i = 0; i < n; i++)
  {
    Halfedge_const_handle he = env.halfedges_begin();
    int q_1 = i-1;
    if (q_1==-1)
      q_1 = n-1;
    
    while (he->source()->point() != points[q_1] || he->target()->point() != points[i])
      he++;

    TEV tev(env);
    Face_handle fh = tev.compute_visibility(points[i], he, regular_output);
    /*
    RSPV regular_visibility(env);
    regular_visibility.compute_visibility(points[i], he, regular_output);
    */

    vp.push_back(arrangement_to_polygon(regular_output));
  }

  return vp;
}

void read_visibility_set(string &filename, vector<set<Point_2>> &PV, int &n)
{
  string text;
  ifstream MyReadFile(filename);
  int i = 0;
  while(getline (MyReadFile, text))
  {
    if(i == n) break;

    vector<string> line = split(text, ";");
    for (string t: line)
    {
      vector<string> tacka = split(t, ",");
      int x = 0, y= 1;

      stringstream geekN(tacka[0]); 
      geekN >> x;
      stringstream geekD(tacka[1]);
      geekD >> y;

      PV[i].insert(Point_2((double)x, (double)y));
    }
    i++;
  }
  MyReadFile.close();
}

void write_test(const string &filename, string &tekst)
{
  std::ofstream outfile;

  outfile.open(filename, std::ios_base::app); // append instead of overwrite
  outfile << tekst; 
}

template<class Kernel, class Container>
Traits_2::FT area_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  Traits_2::FT res = 0;
  if (! pwh.is_unbounded()) {
    res = pwh.outer_boundary().area();
  } else
    return res;

  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) 
    res -= (*hit).area();

  return res;
}

Traits_2::FT area_set_polygons(Polygon_set_2 S)
{
  Traits_2::FT povrsina = 0;

  std::list<Polygon_with_holes_2> res;
  S.polygons_with_holes (std::back_inserter (res));

  Pwh_list_2::const_iterator it_res;
  for (it_res = res.begin(); it_res != res.end(); ++it_res) 
    povrsina += area_polygon_with_holes(*it_res);
      
  return povrsina;
}

// Markove funkcije  -  malo modifikovano

int f(vector<set<Point_2>> &C)
{
  std::set<Point_2> d;
       for(set<Point_2> di : C){
     
           for (Point_2 dii: di){
               //cout << "dii" << dii << endl;
               d.insert(dii);
           }
        
       }
       //cout <<"D" << d.size() << endl; 
       //for(int diii: d)
       //  cout << diii << "\t" <<endl;

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

float greedy_criterion1(int i)
{
    return Cost[i] / Surface[i];

}

/**   Gridi funkcija (6)  prema indeksiranju u radu
TODO: Milane, Za svaka dva cuvara i,j, potrebna mi je struktura koja racuna povrsinu poligona koji je presjek poligona vidljivosti cvora i te poligona vidljivosti cvora j
--> cuvaj ih po mogucnosti u nekoj matricnoj strukturi Intersection[i][j]

Param:
@i: indeks cvora za koji se racuna gridi vrijednost 
@indeks: trenutni parcijalni skup  --> samo indeksi skupova uzeti radi optimizacije
@S: skup svih tjemena
 **/

float greedy_criterion2(vector<set<Point_2>>& S, vector<int>& indeks, int i)   //vector<set<Point_2>>& C, int i)
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
float greedy_criterion(vector<set<Point_2>>& S, int i, vector<set<Point_2>>& C) // take s_i from S
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

int min_greedy(vector<set<Point_2>>& S, vector<set<Point_2>>& C, vector<int>& indeks)
{         //cout << "min_greedy" << endl;
          float g_m = 10000000; int dodaj;
          for(int i = 0; i < S.size(); ++i)
          {  // cout << "i " << i << endl;
              //float g_mi = greedy_criterion(S, i, C); //cout << "gmi: " << g_mi << endl;
              //float g_mi = greedy_criterion1(i);
              float g_mi = greedy_criterion2(S, indeks, i);
              if(g_mi < g_m and g_mi != INFEASIBLE and !findA(indeks, i)) 
              { 
                 dodaj = i;
                 g_m = g_mi;   
              }
          }//cout << "dodati....." << dodaj << endl;
          return dodaj;
}

float greedy_procedure(vector<set<Point_2>> &S)
{
     vector<set<Point_2>> C;
     vector<int> indeks;

     int f_S = f(S); //cout << "f_S: " << f_S << endl;
     int f_C = 0;

     while(f_C != f_S) 
     { 
           int index_set = min_greedy(S, C, indeks);
           
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

int main(int argc, char const *argv[])
{
  for(int file_number=50; file_number<=200; file_number+=2) // start = 8
  for(int file_order=1; file_order<=1; file_order++) // file_order<=file_number
  {
    const string category = "small";
    const string result_location = "test-small-greedy-criterion2.txt";

    string filename = "min-" + std::to_string(file_number) + "-" + std::to_string(file_order) + ".pol";
    string location = "instance/" + category + "/" + filename;
    string predprocesing = "predprocesing/" + category + "/" + filename;
  
    // ----------------------------read polygon---------------------------------
    int n;
    Point_2* points = read_point_of_polygon(location, &n);
    Polygon_2 p = Polygon_2(points, points+n);
    cout<<"Created polygon "<<filename<<endl;

    // -------------------------visibility points------------------------------
    cout<<"Creating visibility polygons...................";
    vector<Polygon_2> pv = visibility_polygon(points, n);
    // sve tacke ce biti orjentisane surotno kazaljki na satu
    for(size_t i = 0; i < n; i++)
      if(pv[i].orientation()==-1)
        pv[i].reverse_orientation();
    // kraj promjene orjentacije
    cout<<"end!!!"<<endl;

    //-------------------------------- area of poligons -------------------------
    cout<<"Calculating area of polygon...................";
    for (size_t i = 0; i < n; i++)
      Surface.push_back((float)converter(pv[i].area()));
    cout<<"end!!!"<<endl;
    
    //-------------------------------- area of intersection poligons -------------------------
    cout<<"Calculating area of digerence polygon...................";
    for (size_t i = 0; i < n; i++)
    {
      vector<float> area_row(n);
      Polygon_set_2 pvs(pv[i]);
    
      for (size_t j = 0; j < n; j++)
        if(i==j)
          area_row[j] = Surface[i];
        else
        {
          Polygon_set_2 pvs(pv[i]);
          pvs.intersection(pv[j]);
          area_row[j] = (float)converter(area_set_polygons(pvs));
        }
        Intersection.push_back(area_row);
    }
    cout<<"end!!!"<<endl;



    // --------------------------read visibility set----------------------------
    cout<<"Creating discretization for visibility polygons..................."<<endl;
    vector<set<Point_2>> pvD(n);
    read_visibility_set(predprocesing, pvD, n);
    cout<<"end!!!"<<endl;

    // ---------------------------------greedy------------------------------------
    // vector<int> cost;
    for (size_t i = 0; i < n; i++)
         Cost.push_back(1);

    auto start = high_resolution_clock::now();
    float s = greedy_procedure(pvD);
    auto stop = high_resolution_clock::now();
    // ------------------------------greedy - end------------------------------------

    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution>: " << duration.count() << " microseconds" << endl;
    string result_write = std::to_string(file_number) + "-" + std::to_string(file_order) + ";" + std::to_string((int)s) + ";" + std::to_string(duration.count()) + "\n";
    write_test(result_location, result_write);
    cout<<"Dovoljan broj cuvara je: "<<s<<endl;
  
  }
  return EXIT_SUCCESS;
}
