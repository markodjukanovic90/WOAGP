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

typedef Arrangement_2::Face_handle                              Face_handle;
// Define the used visibility class 
typedef CGAL::Triangular_expansion_visibility_2<Arrangement_2>  TEV;

CGAL::Cartesian_converter<Kernel,K> converter; 
CGAL::Cartesian_converter<K, Kernel> converter2; 

using namespace std;
using namespace std::chrono;

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

double* resolution(Polygon_2& p)
{
  double* dxy = new double[2];
  dxy[0] = 999999;
  dxy[1] = 999999;
  unsigned int n = p.size();
  double d;

  double absx = abs(CGAL::to_double(p[n-1].x()) - CGAL::to_double(p[0].x()));
  if (absx > 0)
    dxy[0] = absx;

  double absy = abs(CGAL::to_double(p[n-1].y()) - CGAL::to_double(p[0].y()));
  if (absy > 0)
    dxy[1] = absy;
  
  for (size_t i = 0; i < n-1; i++)
  {
    d = abs(CGAL::to_double(p[i].x()) - CGAL::to_double(p[i+1].x()));
    if (d != 0 and d<dxy[0])
      dxy[0] = d;

    d = abs(CGAL::to_double(p[i].y()) - CGAL::to_double(p[i+1].y()));
    if (d != 0 and d<dxy[1])
      dxy[1] = d;
  }
  
  return dxy;
}

vector<Point_2> discrate(Polygon_2& p)
{
  vector<Point_2> points;
  for (size_t i = 0; i < p.size(); i++)
    points.push_back(p[i]); 

  CGAL::Bbox_2 box = p.bbox();
  double x = box.xmin();
  double y = box.ymin();
  double xmax = box.xmax();
  double ymax = box.ymax();
  double* resolutionxy = resolution(p);

  while (x <= xmax)
  {
    y = box.ymin();
    while (y <= ymax)
    {
      
      if (p.is_simple())
      {
        Point_2 point(x,y);
        int bs = p.bounded_side(point);

        if(bs==CGAL::ON_BOUNDED_SIDE) //  || bs==CGAL::ON_BOUNDARY
          points.push_back(point);
      }
      y += resolutionxy[1];
    }
    x += resolutionxy[0];
  } 
  
  return points;
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




template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;
  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}

template<class Kernel, class Container>
void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = "; 
    print_polygon (pwh.outer_boundary());
  }
}

void print_set_poligon(Polygon_set_2 S)
{
  std::list<Polygon_with_holes_2> res;
  S.polygons_with_holes (std::back_inserter (res));

  Pwh_list_2::const_iterator it_res;
  for (it_res = res.begin(); it_res != res.end(); ++it_res) 
    print_polygon_with_holes(*it_res);
}

void write_test(string location, string tekst)
{
  std::ofstream outfile;

  outfile.open(location, std::ios_base::app); // append instead of overwrite
  outfile << tekst;
}

template<class Kernel, class Container>
Traits_2::FT area_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  Traits_2::FT res = 0;
  if (! pwh.is_unbounded()) {
    res = abs(pwh.outer_boundary().area());
  } else
    return res;
/*
  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) 
    res -= abs((*hit).area());
*/
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

int main(int argc, char const *argv[])
{
  for(int file_number=8; file_number<201; file_number+=2) // start = 8
  for(int file_order=1; file_order<=1; file_order++) // file_order<=file_number
  {
    const string category = "small";
    string filename = "min-" + std::to_string(file_number) + "-" + std::to_string(file_order);
    
    string location = "instance/" + category + "/" + filename + + ".pol";
    string output = "preprocessing/" + category + "/" + filename + "_D1.txt";
    string timeexecution = "executiontime-" + category + "_D1.txt";

    // ----------------------------read polygon---------------------------------
    auto start = high_resolution_clock::now();
    int n;
    Point_2* points = read_point_of_polygon(location, &n);
    Polygon_2 p = Polygon_2(points, points+n);
    cout<<"Created polygon "<<filename<<endl;

    // -------------------------visibility polygons------------------------------
    cout<<"Creating visibility polygons...................";
    vector<Polygon_2> pv = visibility_polygon(points, n);
    // sve tacke ce biti orjentisane surotno kazaljki na satu
    for(size_t i = 0; i < pv.size(); i++)
      if(pv[i].orientation()==-1)
        pv[i].reverse_orientation();
    // kraj promjene orjentacije
    cout<<"end!!!"<<endl;

    // -----------------------diskretization polygon------------------------------
    cout<<"Creating diskretization polygon P..............";
    vector<Point_2> D = discrate(p);
    cout<<"end!!!"<<endl;

    // -----------------------diskretization visibility polygons-------------------
    cout<<"Creating diskretization for visibility polygons...................";
    vector<set<Point_2>> pvD(p.size());
    for (size_t i = 0; i < pv.size(); i++)
        discratePVI(D, pv[i], pvD[i]);

    //-------------------------------- area of poligons -------------------------
    vector<float> Surface;
    cout<<"Calculating area of polygon...................";
    for (size_t i = 0; i < n; i++)
      Surface.push_back((float)converter(pv[i].area()));
    cout<<"end!!!"<<endl;
    
    //-------------------------------- area of intersection poligons -------------------------
    cout<<"Calculating area of digerence polygon..................."<<endl;
    vector<vector<float>>Intersection(n);
    for (size_t i = 0; i < n; i++)
      for (size_t j = i+1; j < n; j++)
        {
          Polygon_set_2 pvs(pv[i]);
          pvs.intersection(Polygon_set_2(pv[j]));
          //print_set_poligon(pvs);
          //cout<<p[i]<<p[j]<<"="<<(float)converter(area_set_polygons(pvs))<<endl;
          Intersection[i].push_back( (float)converter(area_set_polygons(pvs)) );
        }
  
    cout<<"end!!!"<<endl;

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution>: " << duration.count() << " microseconds" << endl;
    write_test(timeexecution, filename + ";" + std::to_string(duration.count()) + "\n");
    // ------------------------------------outut----------------------------------
    cout<<"Write output..................."<<endl;
    write_test(output, std::to_string(n) + "\n");
    for (size_t i = 0; i < n; i++)
    {
        string tjeme = "(" + std::to_string((int)(converter(p[i]).x())) + "," + 
          std::to_string((int)(converter(p[i]).y())) + ")";
        
        string vidljivost = "";
        for (Point_2 xy : pvD[i])
          vidljivost += "(" + std::to_string((int)(converter(xy).x())) + "," + 
            std::to_string((int)(converter(xy).y())) + ")\t";
        
        write_test(output, tjeme + ":" + vidljivost + "\n");
    }
    for (size_t i = 0; i < n; i++)
    {
        string tjeme = "(" + std::to_string((int)(converter(p[i]).x())) + "," + 
          std::to_string((int)(converter(p[i]).y())) + ")";
        
        write_test(output, tjeme + ":" + std::to_string(Surface[i]) + "\n");
    }

    for (size_t i = 0; i < n; i++)
    {
      for (size_t j = i+1; j < n; j++)
      {
        string tjemeI = "(" + std::to_string((int)(converter(p[i]).x())) + "," + 
          std::to_string((int)(converter(p[i]).y())) + ")";
        
        string tjemeJ = "(" + std::to_string((int)(converter(p[j]).x())) + "," + 
          std::to_string((int)(converter(p[j]).y())) + ")";
        
        write_test(output, tjemeI + "," + tjemeJ + ":" + std::to_string(Intersection[i][j-i-1]) + "\n");
      }
    }
    cout<<"end!!!"<<endl;
  }

  return EXIT_SUCCESS;
}