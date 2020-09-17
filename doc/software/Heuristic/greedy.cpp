#include <list>
// #include <minisat/core/Solver.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <chrono>
#include <bits/stdc++.h>
#include <utility>
#include <string>
#include <math.h>       /* sqrt */

using namespace std;
using namespace std::chrono;

typedef pair<float, float> Point_2;

int temp = 0;
/** pomocne strukture za gridi metod **/
map<Point_2,vector<Point_2>> Visibility;//for each point, we provide the list od Guards which see that point
map<Point_2,int> numberOfGuards;//for each point we keep the number of Guards in solution which see that point
set<Point_2> CoveredPoints;//set of points covered by the solution
vector<int> indeks; //solution
int alg = 0; // if alg = 1: we execute a Greedy+LS procedure
bool turn_ls = 0; // tunr on LS into Greedy+LS

/** pomocne strukture za gridi metod **/

vector<float> Surface; // for a vertex indexed by i, we return the surface of V(i)
vector<float> Cost; // at index i --> w_i
vector<vector<float>>Intersection; // Intersection[i][j] = Volume of V(i) intersected by V(j)
vector<set<Point_2>> S; // i-th position of S is a set of all points that are visible from vertex i
vector<Point_2> Vertices; // vector of vertices
vector<float> avg_visi; // average visibility vertex
#define INFEASIBLE 1000000
int cardinalityD; // |D(P)|
float wTotal = 0; // total weight 
int t_lim = 0; // time limit
int greedy = 0; // tip gridija => 0: Lovasz; 1: price-per-unit; 2: intersection-based; 3: Greedy by Dragan
std::string path;
std::string output = "";
int n = 0; // |V(P)|
int w_type = 0; // ti ptezine koju pozivamo; 0: tezina proporcionalna velicini vidljivosti svakog vrha; 1: --; 2: --
/** kraj pomocnih str. za gridi metode **/


inline int stoi(string &s) {

     return atoi(s.c_str());
}

inline double stof(string &s) {

     return atof(s.c_str());
}


ostream& operator<<(ostream&ostr, Point_2 p){
	ostr<<"("<<p.first<<","<<p.second<<")";
	return ostr;
}

bool findA(vector<int>& s, int a)
{

     for(int si : s)
         if(si == a)
            return true;

      return false;  
}

float dist(int i, int j)
{
     return sqrt( pow( Vertices[i].first - Vertices[j].first, 2) + pow( Vertices[ i ].second - Vertices[ j ].second, 2) );
}

float distP(Point_2 x, Point_2 y)
{
     return sqrt( pow( x.first - y.first, 2) + pow( x.second - y.second, 2) );
}

void updateCoveredPointsAdd(set<Point_2>& CovPoints,map<Point_2,int> & noOfG, int v){//we add all points which are covered by v, when v is added to solution
    for(Point_2 p : S[v]){
    	
       CovPoints.insert(p);
        if (noOfG.find(p) == noOfG.end()){
       		noOfG.insert(std::make_pair(p,1));//first guard which cover point p
       		//cout<<++temp<<"  First time: "<<p<<endl;
       	}
       	else{
       			//p was already covered
			   noOfG.find(p)->second++;
       		
       	}
       	
       
   }

}
void updateCoveredPointsRemove(set<Point_2>& CovPoints,map<Point_2,int> & noOfG, int v){//we remove all points which are covered by v, when v is removed from solution
    for(Point_2 p : S[v]){
    	
        if (noOfG.find(p) == noOfG.end()){
        	cout<<"------------------ERROR:"<<p<<"  "<<noOfG.find(p)->second<<endl;
        	cin.get();
        	return;
        }
       // cout<<"Tacka "<<p<<" pokrivena sa: "<< noOfG.find(p)->second<<endl;
        if(noOfG.find(p)->second == 1)
        {
        		//only the guard v see p in the solution
        		//remove from CovPoints and remove from noOfGuards
        		CovPoints.erase(p);
        		noOfG.erase(p);
        		//cout<<"Point "<<p<< "is no longer covered"<<endl;
        		
        		
		}
		else
		{
			noOfG.find(p)->second--;//one guard less see p
		}
			
       
   }

}


/*create map structure <key, value> 
key is a point from D(P), value is the vector of Polygon vertices which cover the key
*/
void fill_visibility(){
	//for each point from D(P) we create a vector of Polygon vertices which cover the point
	//start with S
	
	for (int i = 0; i < S.size();i++){
		//set S[i]
		for(auto ix: S[i]){
			//look for the point ix
			//return an iterator to the found element (or to the end() if the element was not found)
			std::map<Point_2,vector<Point_2>>::iterator it = Visibility.find(ix);
			if (it != Visibility.end())//if the point already exists in the map add the vertex i
    			it->second.push_back(Vertices[i]);
			else //we add a new element to Visibility> key = ix, value is a vector containing only the point Vertices[i]
				
				{
					vector<Point_2> vec;
					vec.push_back(Vertices[i]);
					Visibility.insert(std::make_pair(ix,vec));
				}
			
		}
					
	}
		
}
void print_visibility(){
	cout<<"Printing visibility..."<<endl;
	for(auto ix: Visibility){
		cout<<"Point: "<<ix.first<<" visible from: "<<endl; //print key
		for(auto iy:ix.second) //print vertices which cover the key
			cout<<iy<<"\t";
		cout<<endl;
		
	}
}

bool LS(float * obj){
	cout<<"LS executed"<<endl;
	if(0)//check if LS should be executed or not. At this moment, we go in LS
		return false;
	/* to do...*/
	vector<int> copy = indeks;//make a copy
	set<Point_2> CoveredPoints2 = CoveredPoints;
	map<Point_2,int> numberOfGuards2 = numberOfGuards;
	float obj_copy = *obj;
	bool global_succ = false;
	bool reset = false;
	int i = 0;
	while(i<copy.size()){
		int vertex = copy[i];
		cout<<"Cvor: "<<vertex<<endl;
		int neigh1 = (vertex+1) % n;
		int neigh2 = (vertex-1) % n; 
               if( neigh2 <  0) neigh2 = n - 1; 

		if(!findA(copy,neigh1))//neigh1 is not in the solution
	        {
			        cout<<"Test1: "<<neigh1<<":  "<<CoveredPoints2.size()<<endl;
				//cin.get();
				updateCoveredPointsRemove(CoveredPoints2,numberOfGuards2,vertex);
				
				//cout<<"Test2: "<<CoveredPoints2.size()<<endl;
				//cin.get();
				updateCoveredPointsAdd(CoveredPoints2,numberOfGuards2,neigh1);
				//updateCoveredPointsAdd(CoveredPoints2,numberOfGuards2,neigh2);
				//cout<<"Test3: "<<CoveredPoints2.size()<<endl;
				if(CoveredPoints2.size()>CoveredPoints.size()){
					//prelazimo odmah u novo rjesenje
					cout<<"uspjelo"<<endl;
					//cin.get();
					//prelazimo odmah u to rjesenje
					//indeks.erase(indeks.begin()+i);
					//indeks.push_back(neigh1);
					CoveredPoints = CoveredPoints2;
					numberOfGuards = numberOfGuards2;
					//prepravljamo i kopiju
					copy.erase(copy.begin()+i);
					copy.push_back(neigh1);
					indeks = copy;
					i = i-1;
					//i = copy.size();
					global_succ = true;
				}
	       }
	       i++;
	}
	return global_succ;	
}


void read_parameters(int argc, char **argv) {

     int iarg=1;
     while (iarg < argc) {
     if (strcmp(argv[iarg],"-f") == 0) path = (argv[++iarg]);
     else if(strcmp(argv[iarg],"-t") == 0) t_lim = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-greedy") == 0) greedy = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-w_type") == 0) w_type = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-alg") == 0) alg = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-turn_ls") == 0) turn_ls = atoi(argv[++iarg]);
     else if(strcmp(argv[iarg],"-l") == 0) output = argv[++iarg];
     else ++iarg;
   }
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

         Vertices.push_back(vertex);

        //izdvjanje skupa tacaka koje se vide iz tog vrha

        std:: string ostatak = line.substr(end+1, line.length() - end);

        std::string delimiter = "\t";

        size_t pos = 0;
        std::string token;
        set<Point_2> vertexSet;

        while ((pos = ostatak.find(delimiter)) != std::string::npos) {
            token = ostatak.substr(0, pos);
            // std::cout << token << std::endl;

            //obrada svakog tokena da bi se dobile koordinate tacaka


         std::string delim3 = ",";

         auto start3 = 0U;
         auto end3 = vrh.find(delim3);

         std:: string  x =  token.substr(start3+1, end3 - start3-1);
         std:: string  y =  token.substr(end3+1, token.length()-end3-2);

         Point_2 tacka;
         tacka.first = stof(x);
         tacka.second = stof(y);

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
  for (int i: index){
    for(Point_2 p : S[i])
       result.insert(p);
  }
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

int f(vector<int> &C)
{
      std::set<Point_2> d;
      for(int ix: C){
          for (Point_2 dii: S[ix]){
              // cout << "dii" << dii.first << endl;
               d.insert(dii);
          }
      }
      return (int)d.size();
}

int f_minus(vector<int>& C, int i) // f(C u {S_i}) - f(S)
{
    //cout <<" f_minus" << C.size() << " - " << ss.size() << endl;
    int count = 0; // ako i \in ss nije niti u jednom od skupova u S, onda count++;
    if(C.size() == 0 ){
       //cout << "before " << S[i].size() - count;
       return S[i].size();
    }
    // cout << "continue" << endl;
    for(auto ix: S[i]) 
    {
        //cout << ix.first << " " << ix.second << endl;
        bool presented = false;
        for(int id : C){
            if(S[id].count(ix) > 0){
               presented = true;
               break; 
            }
        }
        if(presented) 
           count++;
    } 
    return S[i].size() - count; 
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
      int den = S.size() - indeks.size(); cout<< "den: " << den <<" " <<  S.size()  << endl; 
      for(int j = 0; j < S.size(); ++j) {
          if( std::find(indeks.begin(), indeks.end(), j) == indeks.end() and i < j)
          {   // S_j ne smije biti u vec dodanom parcijalnom skupu 
              // cout << i << " " << j << endl;
              num += Intersection[i][j];
          } 
      }
      //cout << "obj: " << num /den << endl;
      return num / den;
}
/**
param:
@S: skup svih skupova (instance) 
@C: parcijalno rjesenje (trenutno) 
@i: index skupa S_i koji se razmatra za dodavanje u parcijalni skup @C
**/
float greedy_criterion(vector<int>& C, int i) // take s_i from S
{   //cout << "greedy_criterion " << endl;
    //set<Point_2> s;
    //s = S[i]; cout << s.size() << endl;
    int f_m = f_minus(C, i);  //cout << "f_m: " << f_m << endl;
    if(f_m == 0)
       return INFEASIBLE;

    float val = ((float)Cost[i]) / (f_minus(C, i));     //cout << "val: " << val << endl; 
    return val;
}

/**
param: 
@indeks: indeks svih skupova koji su u trenutnom parcijalnom rjesenju;
@i: ineks skupa za koji racunamo gridi vrijednost.

**/
float gridi_criterion_dragan(vector<int> &index, int i)
{

      if(index.size() > 0 and std::find( index.begin(), index.end(), i)  != index.end() )
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
      //cout << "obj: " << obj << endl;
      return obj;
}


int min_greedy(vector<int>& indeks)
{        
          cout << "min_greedy" << endl;
          float g_m = INFEASIBLE; 
          int dodaj = -1;
          for(int i = 0; i < S.size(); ++i)
          { 
              float g_mi;
              switch(greedy){

                  case 0: g_mi = greedy_criterion(indeks, i); break;
                  case 1: g_mi = greedy_criterion_1(i); break;
                  case 2: g_mi = greedy_criterion_2(indeks, i) ; break;
                  default: g_mi = gridi_criterion_dragan(indeks, i);
              }
              bool us = g_mi <= g_m; 
              bool us1 = g_mi != INFEASIBLE;

              if(us and us1 and !findA(indeks, i)) 
              { 
                 dodaj = i;
                 g_m = g_mi;   
              }
          } cout << "dodati....." << dodaj << endl;
          return dodaj;
}

float greedy_procedure()
{  
     vector<int> C;
     vector<int> indeks;
     vector<int> Sx; cout << "n: "<< n << endl;
     for(int i=0; i < n; ++i)
         Sx.push_back(i);
     cout << Sx.size() << endl;
     int f_S = f(Sx); cout << "f_S: " << f_S << endl;
     int f_C = 0;

     while(f_C != f_S) 
     { 
           int index_set = min_greedy(indeks); cout <<"index--------------> " << index_set << endl;
           //cout<<pol[index_set]<<endl;
           cout<<"index set: " << index_set<<endl;
           if(!findA(indeks, index_set)){ //jos nije dodan

               //C.push_back((set<int>) S[index_set] );//cout << "dodaj ----> " << index_set << endl;
               indeks.push_back(index_set);
           } 
           
           f_C = f(indeks); 
            
           cout << "f_C--------------->" << f_C << " |C|=" << C.size() << endl;
     }

     float greedy_val = 0;
     for(int ix: indeks){
         
        greedy_val += Cost[ ix ];
        
     }
     cout << "number of guards: " << indeks.size() << endl;
     return greedy_val; 
     //return C.size();
}

// izbacivanje cvora
// izaberemo random cvor i iz parcijalnog rjesenja C (copy)
//pogledamo njegovu okolinu --> cvorove v_i-k,...v_i-1, v_i+1,...,v_i+k
// uzmemo distancu D i uzmemo samo one cvorove koji su na udaljenosti vecoj od D od v_i (skup VC)
// random izbacimo neki iz skupa VC...

int izbaciCvor(){

     int K = 3;
     vector<int> Candidates; // kandidati za izbacivanje
     int i = random() % indeks.size(); // uzmemo i-ti vektor
     int cvor = indeks[ i ];
     int count1 = 0; 
     for(int i = cvor + 1; i < n; ++i) 
     {
         if(std::find(indeks.begin(), indeks.end(), i) != indeks.end()){ // already in indeks
            Candidates.push_back(i); 
         }
         count1++; 
          
         if(count1 == K -  1) 
            break;
     }
     if(count1 < K - 1) 
     {

        for(int i = 0; count1 == (K - 1); ++i)
        {
           if(std::find(indeks.begin(), indeks.end(), i) != indeks.end()){ // already in indeks
            Candidates.push_back(i); 
           }
           count1++;
           if(count1 == K - 1) 
              break;
        } 
     }
     count1 = 0; 
     for(int i = cvor - 1; i >= 0; --i) 
     {
         if(std::find(indeks.begin(), indeks.end(), i) != indeks.end()){ // already in indeks
            Candidates.push_back(i); 
         }
         count1++; 
           
     }
     if(count1 < K - 1) 
     {
        for(int i = n - 1; i >= 0; --i) 
        {
            if(std::find(indeks.begin(), indeks.end(), i) != indeks.end()){ // already in indeks
               Candidates.push_back(i); 
            }
            count1++;
            if(count1 == K - 1) 
               break;
        }
     }  
     // mjerimo udaljenost potencijalnih covrova koje hocemo izbaciti of cvora <cvor>
     // ako je neki do cvorova na udaljenosti D > 0 od cvora, izbacimo ga iz liste kandidata
     float Dist = 2.0; 

     for(vector<int>::iterator it =  Candidates.begin(); it != Candidates.end(); ++it )
     {
         if(dist(cvor, *it) > 2.0)
         {
            Candidates.erase(it); // izbaci <i> iz potencijalnih cvorova za izbacivanje
         } 
     }
     //biraj random cvor iz Candidates 
     int rx = rand() % Candidates.size();
     return Candidates[ rx ];
}

float greedy_LS()
{
	vector<int> C;

     vector<int> Sx; 
	 cout << "n: "<< n << endl;
     for(int i=0; i < n; ++i)
         Sx.push_back(i);
     cout << Sx.size() << endl;
     int f_S = f(Sx); 
	 cout << "f_S: " << f_S << endl;
	 int f_C = 0;
	 float obj_val = 0;
     while(f_C != f_S) 
     { 
           int index_set = min_greedy(indeks); 
           cout <<"index--------------> " << index_set << endl;
           //cout<<pol[index_set]<<endl;
           cout<<"index set: " << index_set<<endl;
           cout<<"Koliko tacaka pokriva: " << S[index_set].size() << endl;
           if(!findA(indeks, index_set)){ //jos nije dodan
           
               //C.push_back((set<int>) S[index_set] );//cout << "dodaj ----> " << index_set << endl;
               //add vertex to solution
		indeks.push_back(index_set);
               //update set of covered points
               updateCoveredPointsAdd(CoveredPoints,numberOfGuards,index_set);//add  
               //update obj_val
                obj_val += Cost[index_set];//add cost to
           } 
           
           f_C = f(indeks); 
            //f_C = CoveredPoints.size();
           cout << "Control: f_C--------------->" << f_C << endl;// " Covered=" << CoveredPoints.size() << endl;
           bool succ = (turn_ls == 1);
           while(succ){
           	succ = LS(&obj_val);//dragan - later def global
           	if(succ){
           		cout<<"LS succeeded"<<endl;
           		//break;
                }
           }
           
                 int k = 1;
		   vector<int> copy = indeks;
		   float obj_val_copy = obj_val;
		   bool success = true;
		   //cout<<"Duzina indeksa: "<<indeks.size()<<endl;
		   while(k<copy.size()/2)
		   {
		   		int i;
		   		//cout<<"izbacujemo:"<< k<<" elemenata."<<endl;
		   		for(i = 0; i < k; i++){//izbacimo k elemenata
		   			int r = rand() % copy.size(); //funkcija za izbacivanje? 
		   			cout<<i<<" Izbacujemo: "<<copy[r]<<endl;
		   			obj_val_copy -= Cost[copy[r]];
		   			copy.erase( copy.begin() + r );
		   			//cout<<copy[r]<<"umanjimo za: "<<Cost[copy[r]]<<endl;		   			
				}
				//cout<<"vrijednost obj_val_copy: "<<obj_val_copy<<endl;
				//vracamo k elemenata
				i = 0;
				while(i < k){
				
			           int index_set = min_greedy(copy); 
           			   cout<<"index set u popravci: " << index_set<<endl;
           			   if(!findA(copy, index_set)){ //jos nije dodan
           			      i++;
				      copy.push_back(index_set);
				      //cout<<index_set<<" uvecamo za: "<<Cost[index_set]<<endl;
                	              obj_val_copy += Cost[index_set];//add cost to
                                }
                }
                cout<<"vraceno:"<< k<<" elemenata. obj_val_copy: "<<obj_val_copy<<endl;
                //cout<<"obj_val: "<<obj_val<<endl;
                if(obj_val_copy<obj_val){
			cout<<"obj_val: "<<obj_val<<endl;
			cout<<"obj_val_copy: "<<obj_val_copy<<endl;
			indeks = copy;
			obj_val = obj_val_copy;
                	cin.get();
                }
                //obj_val_copy = obj_val;
				k++;	
           		
				 
 			}
        }
	/*provjera*/
        float check1 = 0;
	for (int i = 0; i < indeks.size();i++){
		check1 += Cost[i];
	}
        obj_val = check1; // daoan update
	cout<<"Ukupan kost: " << check1 << " Broj strazara: "<<indeks.size()<<endl;
	float check2 = f(indeks); 
    //f_C = CoveredPoints.size();
    cout << "Control: f_C--------------->" << check2 << " Polazno= " << f_S << endl;
    
    cout << "number of guards: " << indeks.size() << endl;
    return obj_val; 
     //return C.size();
	
}




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


bool equalP(Point_2 x, Point_2 y)
{
    if(x.first==y.first and x.second==y.second)
        return true;

    return false;
}

void avg_visi_vertex()
{
    for (size_t i = 0; i < Vertices.size(); i++)
    {
        int n = S[i].size()-1;
        float dists[n];
        size_t ind = 0;

        for(Point_2 point : S[i])
            if(!equalP(Vertices[i], point))
            {
                dists[ind] = distP(Vertices[i], point);
                ind++;
            }

        float avg = 0.0;
        for (size_t j = 0; j < n; j++)
            avg += dists[j];
         
        avg_visi.push_back(avg/n);
    }
}

void write_test(string tekst)
{
  std::ofstream outfile;

  outfile.open(output, std::ios_base::app); // append instead of overwrite
  outfile << tekst;
}


int main( int argc, char **argv ) {

    read_parameters(argc, argv);
    read_from_file(path); // fill Cost, Intersection and Surface
    cout << S.size() << " " << S[0].size() << " " << Surface[0] <<  endl;
    n = Surface.size();
    avg_visi_vertex(); //********************************* create average visibility vertex**************************************

    /*for(set<Point_2>& s: S)
    {
        for(Point_2 sx: s){
           cout << "(" << sx.first << ", " << sx.second << ")" << " ";
        }
        cout << endl;

    }*/

    // ----------------------- dodjela tezina (proporcionalno broju pokrivenih tacaka iz D(P) ------------------------------------
    switch(w_type){
 
          case 0: { // tezina proporcionalna sa velicinom skupa S[i]
                  for(auto& X: S){
                      cardinalityD += X.size();
                  }
                  for (size_t i = 0; i < n; i++){
                      float w_i = n * ( ((float) S[i].size()) / cardinalityD );
                      wTotal += w_i;
                      Cost.push_back(w_i);
                  }; break;}
          case 1: {     
                  float dist0  = (dist(0, n-1) + dist(0,1) ) / 2; 
                  Cost.push_back(dist0); wTotal += dist0;
                   for(int i = 1; i < n; ++i){ 
                      float w_i = ( (0.0 + dist(i-1, i) + dist(i, i+1) ) / 2.0 );
                      Cost.push_back( w_i );
                      wTotal += w_i;
                      //cout << "w_i= " << Cost[i] << endl;
                   }
                     break;}
          case 2: {
                    avg_visi_vertex();
                    for(float c :avg_visi)
                    {
                        Cost.push_back(c); 
                        wTotal += c;
                    } 
                    break;
                  } 
          case 3: { // random weights
                    for(size_t i = 0; i < n; i++) //  weights <= 10
                    {
                        Cost.push_back( rand() % 10);
                        wTotal += Cost[ i ];
                    }
                    break;
                   }
          default: {
                    for(size_t i = 0; i < n; i++) // non-weighted version of the problem
                       Cost.push_back(1);
                    wTotal = Cost.size();
                   }
                
   }   cardinalityD = Cost.size();
    // ---------------------------------greedy------------------------------------
    cout << "Run Greedy" << " with type: " <<  greedy << endl;
    auto start = high_resolution_clock::now();
    float s = 10000000;
    if(alg == 1) 
        s = greedy_LS();
    else 
        s = greedy_procedure();

    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start); 
    cout << "Time Execution: " << duration.count() << " microseconds" << endl;
    cout << "Result: " << s << endl;

    if(output.compare("") != 0){
        string name_polygon = split(split(path, "/")[4], "_")[0];
        //cout<<name_polygon<<"---"<<output<<endl;
        write_test(name_polygon + ";" + std::to_string(s) + ";" + std::to_string(duration.count()) + "\n");
    }
    
    
    return EXIT_SUCCESS;
}
