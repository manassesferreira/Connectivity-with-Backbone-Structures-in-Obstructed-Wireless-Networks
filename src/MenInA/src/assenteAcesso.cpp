#include "common.h"
#include <Rcpp.h>
using namespace Rcpp;

#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <cstdlib>

#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include <lemon/lgf_reader.h>
//#include <lemon/math.h> ///usr/share/R/include/R_ext/Constants.h:36:24: error: expected unqualified-id before numeric constant #define PI             M_PI

#include "test_tools.h"

using namespace std;
using namespace lemon;

GRAPH_TYPEDEFS(SmartGraph);

void checkMatching(const SmartGraph& graph,
                   const MaxMatching<SmartGraph>& mm) {
  int num = 0;

  IntNodeMap comp_index(graph);
  UnionFind<IntNodeMap> comp(comp_index);

  int barrier_num = 0;

  for (NodeIt n(graph); n != INVALID; ++n) {
    check(mm.status(n) == MaxMatching<SmartGraph>::EVEN ||
          mm.matching(n) != INVALID, "Wrong Gallai-Edmonds decomposition");
    if (mm.status(n) == MaxMatching<SmartGraph>::ODD) {
      ++barrier_num;
    } else {
      comp.insert(n);
    }
  }

  for (EdgeIt e(graph); e != INVALID; ++e) {
    if (mm.matching(e)) {
      check(e == mm.matching(graph.u(e)), "Wrong matching");
      check(e == mm.matching(graph.v(e)), "Wrong matching");
      ++num;
    }
    check(mm.status(graph.u(e)) != MaxMatching<SmartGraph>::EVEN ||
          mm.status(graph.v(e)) != MaxMatching<SmartGraph>::MATCHED,
          "Wrong Gallai-Edmonds decomposition");

    check(mm.status(graph.v(e)) != MaxMatching<SmartGraph>::EVEN ||
          mm.status(graph.u(e)) != MaxMatching<SmartGraph>::MATCHED,
          "Wrong Gallai-Edmonds decomposition");

    if (mm.status(graph.u(e)) != MaxMatching<SmartGraph>::ODD &&
        mm.status(graph.v(e)) != MaxMatching<SmartGraph>::ODD) {
      comp.join(graph.u(e), graph.v(e));
    }
  }

  std::set<int> comp_root;
  int odd_comp_num = 0;
  for (NodeIt n(graph); n != INVALID; ++n) {
    if (mm.status(n) != MaxMatching<SmartGraph>::ODD) {
      int root = comp.find(n);
      if (comp_root.find(root) == comp_root.end()) {
        comp_root.insert(root);
        if (comp.size(n) % 2 == 1) {
          ++odd_comp_num;
        }
      }
    }
  }

  check(mm.matchingSize() == num, "Wrong matching");
  check(2 * num == countNodes(graph) - (odd_comp_num - barrier_num),
         "Wrong matching");
  return;
}


// [[Rcpp::export]]
List assenteAcesso(List B, List C, List D, int g, double r, double eps) {
  std::list<int> _c, covered;
  std::list<double> _x, _y;

  NumericVector Bx = as<NumericVector>(B[0]);
  NumericVector By = as<NumericVector>(B[1]);
  NumericVector c1 = as<NumericVector>(B[2]);
  NumericVector c2 = as<NumericVector>(B[3]);
  NumericVector c3 = as<NumericVector>(B[4]);
  NumericVector c4 = as<NumericVector>(B[5]);

  for (int i = 0; i < c1.size(); ++i) {
    if (c4[i]<0 && c3[i]<0 && c2[i]>0){
      //std::cout << i << " " << c1[i] << " " << c2[i] << std::endl;
      _c.push_back(c1[i]);
      _c.push_back(c2[i]);
    }
  }
  int s = _c.size();

  _c.sort();
  _c.unique();

  //cout<<endl; cout<<endl;
  int cont=1;
  int *wisemap;
  wisemap = new int[_c.size()+1];
  std::list<int>::const_iterator c;
  for (c = _c.begin(); c != _c.end(); ++c) {
    //std::cout << *c << " ";
    wisemap[cont]=*c;
    cont++;
  }
  //cout<<endl; cout<<endl;

  cont = 0;
  int t=_c.size();
  int *u,*v;
  u = new int[s];
  v = new int[s];
  for (int i = 0; i < c1.size(); ++i) {
    if (c4[i]<0 && c3[i]<0 && c2[i]>0){
      int a = std::distance(wisemap, std::find(wisemap, wisemap + _c.size()+1, c1[i]));
      int b = std::distance(wisemap, std::find(wisemap, wisemap + _c.size()+1, c2[i]));
      //std::cout << i << " " << cont << " " << a << " " << b << std::endl;
      u[cont]=a;
      v[cont]=b;
      cont++;
    }
  }
  s=cont-1;
/*----------------------------------------------------------------------------*/

  SmartGraph from;
  SmartGraph::NodeMap<int> fnm(from);
  SmartGraph::ArcMap<int> fam(from);
  SmartGraph::EdgeMap<int> fem(from);
  SmartGraph::Node fn = INVALID;
  SmartGraph::Arc fa = INVALID;
  SmartGraph::Edge fe = INVALID;
  std::vector<SmartGraph::Node> fnv;

  for (int i = 0; i < t; ++i) {
    SmartGraph::Node node = from.addNode();
    fnv.push_back(node);
    fnm[node] = i;
    if (i == 0) fn = node;
  }

  int e1,e2;
  for (int i = 0; i < s; ++i) {
    e1=u[i]-1;
    e2=v[i]-1;
    SmartGraph::Edge edge = from.addEdge(fnv[e1], fnv[e2]);
    fem[edge] = e1 * e1 + e2 * e2;
    fam[from.direct(edge, true)] = e1 + e2 * e2;
    fam[from.direct(edge, false)] = e1 * e1 + e2;
    if (e1 == 0 && e2 == 0) fa = from.direct(edge, true);
    if (e1 == 0 && e2 == 0) fe = edge;
  }

  MaxMatching<SmartGraph> mm(from);
  mm.run();

  checkMatching(from, mm);
  /*cout <<"----------------------"<< endl;
  cout << "mm= " << mm.matchingSize() << endl;
  cout << "s= " << t - 2*mm.matchingSize() << endl;
  cout <<"----------------------"<< endl;*/

  int u1,v2,chave;
  bool found;
  if(mm.matchingSize()>0){
    for (EdgeIt e(from); e != INVALID; ++e) {
      if (mm.matching(e)) {

        u1 = wisemap[from.id(from.u(e))+1];
        v2 = wisemap[from.id(from.v(e))+1];

        chave=0; found=false;
        while(!found){
          if( c1[chave] == u1 && c2[chave] == v2){
            found=true;
          }else{
            chave++;
          }
        }

        /*std::cout << " " << from.id(from.u(e)) << "-" << from.id(from.v(e)) <<
           " " << wisemap[from.id(from.u(e))+1] << "-" << wisemap[from.id(from.v(e))+1] <<
           "\t\t" << c1[chave] << "-" << c2[chave] << "\t" << Bx[chave] << "-" << By[chave]
           << std::endl;*/

         _x.push_back(Bx[chave]);
         _y.push_back(By[chave]);
         covered.push_back(u1);
         covered.push_back(v2);
      }
    }
  }

/******************************************************************************/

  std::list<int>::const_iterator it;
  for (int i = 0; i < c1.size(); ++i) {
    if (c4[i]<0 && c3[i]<0 && c2[i]<0){

      found=false;
      it = covered.begin();
      while ( it != covered.end() && !found) {
          //std::cout << *it << std::endl;
        if (*it == c1[i]){
          found=true;
        }else{
          ++it;
        }
      }

      if(!found){
        //std::cout << i  << " " << c1[i] << " " << Bx[i] << " " << By[i] << std::endl;
        _x.push_back(Bx[i]);
        _y.push_back(By[i]);
      }
    }
  }

  delete [] wisemap;
  delete [] u;
  delete [] v;

  return Rcpp::List::create(
    Rcpp::Named("x") = _x, Rcpp::Named("y") = _y
  );
}
