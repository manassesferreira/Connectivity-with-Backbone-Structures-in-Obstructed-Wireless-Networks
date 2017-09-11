#include "common.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List definaDispositivo(int g, int mu) {
 int s=2*(g-1)*g; //o número de segmentos do grid

 std::list<int> seg; //em qual segmento está o dispositivo
 for (int i = 0; i < s; ++i) {
   for (int j = 0; j < mu; ++j) {
     seg.push_back(i);
   }
 }

 std::list<double> pos; //em qual posicao (0,1) está o dispositivo
 for (int i = 0; i < mu*s; ++i) {
   pos.push_back( ((double) rand() / (RAND_MAX)) );
 }

 return Rcpp::List::create(Rcpp::Named("seg") = seg,
                           Rcpp::Named("pos") = pos);
}
