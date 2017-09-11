#include "common.h"
#include "existeConectividade.h"
#include "segmentosDoCruzamento.h"
#include "weighted_quick_union.h"

#include <Rcpp.h>
using namespace Rcpp;

NumericVector getComponente(int Vertices, NumericVector u, NumericVector v) {
  NumericVector out(Vertices);
  int Edges = u.size();
  if(Edges>0){
    UF *uf;
    uf = new weightedQuickUnion(Vertices);
    int node1, node2;
    //cout<<Edges<<" "<<v.size()<<endl;
    cout << " union "  << endl;
    for(int i = 0; i < Edges; ++i) {
      //cout << i << ": " << u[i] << "," << v[i] << endl;
      node1=u[i]; node2=v[i];
      uf->_union(node1, node2);
    }
    cout << " find "  << endl;
    for (int i = 0; i < Vertices; i++) {
      //cout << i << ": " << uf->find(i) << endl;
      out[i] = uf->find(i);
/*
      for (int j = i+1; j < Vertices; j++) {
        if (uf->connected(i, j)){
          cout << "(" << i << "," << j << ") "<< uf->connected(i, j) << endl;
        }
      }
*/
    }
  }
  return out ;
}

// [[Rcpp::export]]
List computeComponente(List D, int g, double r, double eps) {
 cout << " in "  << endl;
 //cout << g << " : " << r << " : " << eps << endl;

 int Dispositivos = as<NumericVector>(D[0]).size();
 int Segmentos = 2 * g * ( g - 1 );
 int s;

 NumericVector seg = as<NumericVector>(D[0]);
 NumericVector pos = as<NumericVector>(D[1]);
 std::list<double> _x;
 std::list<double> _y;

 typedef struct {
   std::list<int> lista;
 } TipoArrayList;
 TipoArrayList *dispositivosNoSegmento;
 dispositivosNoSegmento=new TipoArrayList[Segmentos];

 for (int i = 0; i < Dispositivos; ++i) {
   s = (int)seg[i];
   //std::cout << i << " " << s << " " << pos[i] << std::endl;
   dispositivosNoSegmento[s].lista.push_back(i);
   if( s%(2*g-1) > (g-2) ){ //segmento vertical
     //std::cout << "\t" << (double)(s%(2*g-1)-(g-1)) << " " << pos[i] << std::endl;
     _x.push_back( (double)(s%(2*g-1)-(g-1)) );
     _y.push_back( pos[i] + s/(2*g-1) );
   }else{ //segmento horizontal
     //std::cout << "\t" << pos[i] << " " << (double)(floor(s/(2*g-1))) << std::endl;
     _x.push_back( pos[i] + s%(2*g-1) );
     _y.push_back( (double)(floor(s/(2*g-1))) );
   }
 }

 cout << " _u _v "  << endl;

 NumericVector _u;
 NumericVector _v;

 //nos segmentos
 for (int i =0; i < 2*(g-1)*g; ++i){
   std::list<int>::const_iterator aux,d1,d2;//mas mantenha o respeito
   for (d1 = dispositivosNoSegmento[i].lista.begin();
     d1 != dispositivosNoSegmento[i].lista.end(); ++d1) {
     aux = d1;
     for (d2 = ++aux; d2 != dispositivosNoSegmento[i].lista.end(); ++d2) {
       if ( fabs(pos[*d1] - pos[*d2]) <= r ){
         _u.push_back(*d1);
         _v.push_back(*d2);
       }
     }
   }
 }

//nos cruzamentos
 for (int i = 0; i < g*g; ++i) {
   //std::cout << i << " : ";
   int borda; //1-OS,2-S,3-LS,4-O,5-L,6-NO,7-N,8-NL
   int N, L, O, S;
   int xCruz = i%g;
   int yCruz = floor(i/g);

   //cout << "\t segmentosDoCruzamento "  << endl;
   segmentosDoCruzamento(g,xCruz,yCruz,&N,&L,&O,&S,&borda);

   std::list<int> dispositivosNoCruzamento;
   dispositivosNoCruzamento.insert(
     dispositivosNoCruzamento.end(),
     dispositivosNoSegmento[N].lista.begin(),
     dispositivosNoSegmento[N].lista.end()
   );
   dispositivosNoCruzamento.insert(
     dispositivosNoCruzamento.end(),
     dispositivosNoSegmento[L].lista.begin(),
     dispositivosNoSegmento[L].lista.end()
   );
   dispositivosNoCruzamento.insert(
     dispositivosNoCruzamento.end(),
     dispositivosNoSegmento[O].lista.begin(),
     dispositivosNoSegmento[O].lista.end()
   );
   dispositivosNoCruzamento.insert(
     dispositivosNoCruzamento.end(),
     dispositivosNoSegmento[S].lista.begin(),
     dispositivosNoSegmento[S].lista.end()
   );

   //std::copy(std::begin(dispositivosNoCruzamento), std::end(dispositivosNoCruzamento),
   //        std::ostream_iterator<int>(std::cout, " "));
   //std::cout << std::endl;
   //std::cout << std::endl;

   std::list<int>::const_iterator d1;
   std::list<double>::iterator d1x, d1y;

   std::list<int> proximosAoCruzamento;
   for (d1 = dispositivosNoCruzamento.begin();
     d1 != dispositivosNoCruzamento.end(); ++d1) {
     d1x = _x.begin(); d1y = _y.begin();
     std::advance(d1x, *d1); std::advance(d1y, *d1);

     if ( existeConectividade(g, borda, r, eps, *d1x, *d1y, (double)xCruz, (double)yCruz) ){
       proximosAoCruzamento.push_back(*d1);
     }
   }

   //std::cout << std::endl;
   //std::copy(std::begin(proximosAoCruzamento), std::end(proximosAoCruzamento),
   //          std::ostream_iterator<int>(std::cout, " "));
   //std::cout << std::endl;
   //std::cout << std::endl;

   //cout << "\t d1 d2 "  << endl;
   std::list<int>::const_iterator d2,aux;
   std::list<double>::iterator d2x, d2y;
   for (d1 = proximosAoCruzamento.begin();
     d1 != proximosAoCruzamento.end(); ++d1) {
     aux = d1;
     for (d2 = ++aux; d2 != proximosAoCruzamento.end(); ++d2) {
       //std::cout << *d1 <<  " " << *d2 << std::endl;

       d1x = _x.begin(); d1y = _y.begin(); std::advance(d1x, *d1); std::advance(d1y, *d1);
       //std::cout << "\t" << *d1x <<  " " << *d1y << std::endl;

       d2x = _x.begin(); d2y = _y.begin(); std::advance(d2x, *d2); std::advance(d2y, *d2);
       //std::cout << "\t" << *d2x <<  " " << *d2y << std::endl;

       if ( existeConectividade(g, borda, r, eps, *d1x, *d1y, *d2x, *d2y) ) {
         //std::cout << "\t" << *d1 <<  " " << *d2 << std::endl;
         _u.push_back( *d1 );
         _v.push_back( *d2 );
       }
     }
   }
   //std::cout << std::endl;
 }

 cout << " _C "  << endl;
 NumericVector _C = getComponente(Dispositivos, _u, _v); //using union-find method

 delete [] dispositivosNoSegmento;
 return Rcpp::List::create(Rcpp::Named("x") = _x,
                           Rcpp::Named("y") = _y,
                           Rcpp::Named("u") = _u,
                           Rcpp::Named("v") = _v,
                           Rcpp::Named("comp") = _C);
}
