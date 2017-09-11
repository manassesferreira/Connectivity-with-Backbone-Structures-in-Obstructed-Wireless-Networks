#include "common.h"
#include "existeConectividade.h"
#include "segmentosDoCruzamento.h"
#include <algorithm>    // std::max
#include <Rcpp.h>
using namespace Rcpp;

typedef struct {
  std::list<int> lista;
} TipoArrayList;

typedef struct {
  int d,c;
  double x,y;
} TipoDisp;
bool by_x ( const TipoDisp &a, const TipoDisp &b ) {
  return a.x < b.x;
}
bool by_y ( const TipoDisp &a, const TipoDisp &b ) {
  return a.y < b.y;
}

bool cabeUmAcesso(int g, double r, double eps, double xCruz, double yCruz,
  double Nx, double Ny, double Lx, double Ly, double *Bx, double *By){
  //std::cout << "\t" << Nx << " " << Ny << " - " << Lx << " " << Ly << std::endl;
  if( Ny-yCruz >= Lx-xCruz ){
    *Bx = Nx;
    *By = ( yCruz>Ny-r ? yCruz : Ny-r );
    //std::cout << "\t\t" << " N " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Lx, Ly) ) return true;
  }else{
    *Bx = ( xCruz>Lx-r ? xCruz : Lx-r );
    *By = Ly;
    //std::cout << "\t\t" << " L " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) ) return true;
  }
  return false;
}

bool cabeUmAcessoNLO(int g, double r, double eps, double xCruz, double yCruz,
  double Nx, double Ny, double Lx, double Ly, double Ox, double Oy, double *Bx, double *By){
  //std::cout << "\t" << Nx << " " << Ny << " - " << Lx << " " << Ly << " - " << Ox << " " << Oy << std::endl;
  if( Ny-yCruz >= Lx-xCruz && Ny-yCruz >= xCruz-Ox ){
    *Bx = Nx;
    *By = ( yCruz>Ny-r ? yCruz : Ny-r );
    //std::cout << "\t\t" << " N " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Lx, Ly) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) ) return true;
  }else if( Lx-xCruz >= Ny-yCruz && Lx-xCruz >= xCruz-Ox ){
    *Bx = ( xCruz>Lx-r ? xCruz : Lx-r );
    *By = Ly;
    //std::cout << "\t\t" << " L " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) ) return true;
  }else if( xCruz-Ox >= Ny-yCruz && xCruz-Ox >= Lx-xCruz  ){
    *Bx = ( r>xCruz-Ox ? xCruz : Ox+r );
    *By = Oy;
    //std::cout << "\t\t" << " O " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) ) return true;
  }
  return false;
}

bool cabeUmAcessoNLS(int g, double r, double eps, double xCruz, double yCruz,
  double Nx, double Ny, double Lx, double Ly, double Sx, double Sy, double *Bx, double *By){
  //std::cout << "\t" << Nx << " " << Ny << " - " << Lx << " " << Ly << " - " << Sx << " " << Sy << std::endl;
  if( Ny-yCruz >= Lx-xCruz && Ny-yCruz >= yCruz-Sy ){
    *Bx = Nx;
    *By = ( yCruz>Ny-r ? yCruz : Ny-r );
    //std::cout << "\t\t" << " N " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Lx, Ly) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( Lx-xCruz >= Ny-yCruz && Lx-xCruz >= yCruz-Sy ){
    *Bx = ( xCruz>Lx-r ? xCruz : Lx-r );
    *By = Ly;
    //std::cout << "\t\t" << " L " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( yCruz-Sy >= Ny-yCruz && yCruz-Sy >= Lx-xCruz  ){
    *Bx = Sx;
    *By = ( r>yCruz-Sy ? yCruz : Sy+r );
    //std::cout << "\t\t" << " S " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) ) return true;
  }
  return false;
}

bool cabeUmAcessoNOS(int g, double r, double eps, double xCruz, double yCruz,
  double Nx, double Ny, double Ox, double Oy, double Sx, double Sy, double *Bx, double *By){
  //std::cout << "\t" << Nx << " " << Ny << " - " << Ox << " " << Oy << " - " << Sx << " " << Sy << std::endl;
  if( Ny-yCruz >= xCruz-Ox && Ny-yCruz >= yCruz-Sy ){
    *Bx = Nx;
    *By = ( yCruz>Ny-r ? yCruz : Ny-r );
    //std::cout << "\t\t" << " N " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( xCruz-Ox >= Ny-yCruz && xCruz-Ox >= yCruz-Sy ){
    *Bx = ( r>xCruz-Ox ? xCruz : Ox+r );
    *By = Oy;
    //std::cout << "\t\t" << " O " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( yCruz-Sy >= Ny-yCruz && yCruz-Sy >= xCruz-Ox  ){
    *Bx = Sx;
    *By = ( r>yCruz-Sy ? yCruz : Sy+r );
    //std::cout << "\t\t" << " S " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Ox, Oy, *Bx, *By) ) return true;
  }
  return false;
}

bool cabeUmAcessoLOS(int g, double r, double eps, double xCruz, double yCruz,
  double Lx, double Ly, double Ox, double Oy, double Sx, double Sy, double *Bx, double *By){
  //std::cout << "\t" << Lx << " " << Ly << " - " << Ox << " " << Oy << " - " << Sx << " " << Sy << std::endl;
  if( Lx-xCruz >= xCruz-Ox && Lx-xCruz >= yCruz-Sy ){
    *Bx = ( xCruz>Lx-r ? xCruz : Lx-r );
    *By = Ly;
    //std::cout << "\t\t" << " L " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( xCruz-Ox >= Lx-xCruz && xCruz-Ox >= yCruz-Sy ){
    *Bx = ( r>xCruz-Ox ? xCruz : Ox+r );
    *By = Oy;
    //std::cout << "\t\t" << " O " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( yCruz-Sy >= Lx-xCruz && yCruz-Sy >= xCruz-Ox  ){
    *Bx = Sx;
    *By = ( r>yCruz-Sy ? yCruz : Sy+r );
    //std::cout << "\t\t" << " S " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Ox, Oy, *Bx, *By) ) return true;
  }
  return false;
}

bool cabeUmAcessoNLOS(int g, double r, double eps, double xCruz, double yCruz, double Nx, double Ny,
  double Lx, double Ly, double Ox, double Oy, double Sx, double Sy, double *Bx, double *By){
  //std::cout << "\t" << Nx << " " << Ny << " - " << Lx << " " << Ly << " - " << Ox << " " << Oy << " - " << Sx << " " << Sy << std::endl;
  if( Ny-yCruz >= Lx-xCruz && Ny-yCruz >= xCruz-Ox && Ny-yCruz >= yCruz-Sy ){
    *Bx = Nx;
    *By = ( yCruz>Ny-r ? yCruz : Ny-r );
    //std::cout << "\t\t" << " N " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, *Bx, *By, Lx, Ly) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( Lx-xCruz >= Ny-yCruz && Lx-xCruz >= xCruz-Ox && Lx-xCruz >= yCruz-Sy ){
    *Bx = ( xCruz>Lx-r ? xCruz : Lx-r );
    *By = Ly;
    //std::cout << "\t\t" << " L " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Ox, Oy) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( xCruz-Ox >= Ny-yCruz && xCruz-Ox >= Lx-xCruz && xCruz-Ox >= yCruz-Sy ){
    *Bx = ( r>xCruz-Ox ? xCruz : Ox+r );
    *By = Oy;
    //std::cout << "\t\t" << " O " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, *Bx, *By, Sx, Sy) ) return true;
  }else if( yCruz-Sy >= Ny-yCruz && yCruz-Sy >= Lx-xCruz && yCruz-Sy >= xCruz-Ox  ){
    *Bx = Sx;
    *By = ( r>yCruz-Sy ? yCruz : Sy+r );
    //std::cout << "\t\t" << " S " << *Bx << " " << *By << std::endl;
    if ( existeConectividade(g, 0, r, eps, Nx, Ny, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Lx, Ly, *Bx, *By) &&
      existeConectividade(g, 0, r, eps, Ox, Oy, *Bx, *By) ) return true;
  }
  return false;
}

int maisProximoAoCruzamentoDoSegmento(int g, char S, int s, int borda,
  TipoArrayList *dispositivosNoSegmento, NumericVector comp, NumericVector Cx, NumericVector Cy,
  int *Sd, int *Sc, double *Sx, double *Sy){
  TipoDisp *disps;
  disps=new TipoDisp[dispositivosNoSegmento[s].lista.size()];
  std::list<int>::const_iterator d1;
  int cont=0;
  for (d1 = dispositivosNoSegmento[s].lista.begin();
    d1 != dispositivosNoSegmento[s].lista.end(); ++d1) {
    disps[cont].d = *d1;
    disps[cont].c = (int)comp[*d1];
    disps[cont].x = (double)Cx[*d1];
    disps[cont].y = (double)Cy[*d1];
    cont++;
  }
  if( s%(2*g-1) > (g-2) ){ //segmento vertical
    std::sort(disps,disps+cont,by_y);
  }else{ //segmento horizontal
    std::sort(disps,disps+cont,by_x);
  }

  if (S == 'N' || S == 'L'){
    *Sc=disps[0].c;
    *Sd=disps[0].d;
    *Sx=disps[0].x;
    *Sy=disps[0].y;
  } else if (S == 'O' || S == 'S'){
    *Sc=disps[dispositivosNoSegmento[s].lista.size()-1].c;
    *Sd=disps[dispositivosNoSegmento[s].lista.size()-1].d;
    *Sx=disps[dispositivosNoSegmento[s].lista.size()-1].x;
    *Sy=disps[dispositivosNoSegmento[s].lista.size()-1].y;
  }

  switch ( borda ) {
    case 1: //OS
      if (*Sx > 1 ) *Sx = *Sx - (g-1);//O
      if (*Sy > 1 ) *Sy = *Sy - (g-1);//S
      break;
    case 2: //S
      if (*Sy > 1 ) *Sy = *Sy - (g-1);//S
      break;
    case 3: //LS
      if (*Sx <= 1 ) *Sx = *Sx + (g-1);//L
      if (*Sy > 1 ) *Sy = *Sy - (g-1);//S
      break;
    case 4: //O
      if (*Sx > 1 ) *Sx = *Sx - (g-1);//O
      break;
    case 5: //L
      if (*Sx <= 1 ) *Sx = *Sx + (g-1);//L
      break;
    case 6: //NO
      if (*Sy <= 1 ) *Sy = *Sy + (g-1);//N
      if (*Sx > 1 ) *Sx = *Sx - (g-1);//O
      break;
    case 7: //N
      if (*Sy <= 1 ) *Sy = *Sy + (g-1);//N
      break;
    case 8: //NL
      if (*Sy <= 1 ) *Sy = *Sy + (g-1);//N
      if (*Sx <= 1 ) *Sx = *Sx + (g-1);//L
      break;
  }

  delete [] disps;
}

// [[Rcpp::export]]
List busqueBase(List C, List D, int g, double r, double eps) {
  std::list<int> _c1, _c2, _c3, _c4, _c;
  std::list<double> _x, _y;
  int Dispositivos = as<NumericVector>(D[0]).size();
  int Segmentos = 2 * g * ( g - 1 );
  int s, c;

  NumericVector seg = as<NumericVector>(D[0]);
  NumericVector pos = as<NumericVector>(D[1]);

  NumericVector Cx = as<NumericVector>(C[0]);
  NumericVector Cy = as<NumericVector>(C[1]);
  NumericVector comp = as<NumericVector>(C[4]);

  TipoArrayList *dispositivosNoSegmento;
  dispositivosNoSegmento=new TipoArrayList[Segmentos];

  TipoArrayList *componentesNoSegmento;
  componentesNoSegmento=new TipoArrayList[Segmentos];

  for (int i = 0; i < Dispositivos; ++i) {
    s = (int)seg[i];
    c = (int)comp[i];
    //std::cout << i << " " << s << " " << c << std::endl;
    dispositivosNoSegmento[s].lista.push_back(i);
    componentesNoSegmento[s].lista.push_back(c);
    _c.push_back(c);
  }

  //dois a dois - nos segmentos
  for (int i =0; i < 2*(g-1)*g; ++i){

    std::list<int>::const_iterator c1;
    componentesNoSegmento[i].lista.sort();
    componentesNoSegmento[i].lista.unique();
    //std::cout << componentesNoSegmento[i].lista.size() << std::endl;

    if (componentesNoSegmento[i].lista.size() > 1 ){

      TipoDisp *disps;
      disps=new TipoDisp[dispositivosNoSegmento[i].lista.size()];

      std::list<int>::const_iterator d1;
      std::list<double>::iterator *Sx, *Sy;
      int cont=0;
      for (d1 = dispositivosNoSegmento[i].lista.begin();
        d1 != dispositivosNoSegmento[i].lista.end(); ++d1) {
          //std::cout << *d1 << " " << i << " " << (int)comp[*d1] << "\t" << Cx[*d1] << " " << Cy[*d1] << std::endl;

          disps[cont].d = *d1;
          disps[cont].c = (int)comp[*d1];
          disps[cont].x = (double)Cx[*d1];
          disps[cont].y = (double)Cy[*d1];
          cont++;
      }

      if( i%(2*g-1) > (g-2) ){ //segmento vertical
        std::sort(disps,disps+cont,by_y);
      }else{ //segmento horizontal
        std::sort(disps,disps+cont,by_x);
      }

      //std::cout << std::endl; std::cout << std::endl;
      int compAtual=disps[0].c;
      for (int j=1; j < dispositivosNoSegmento[i].lista.size(); ++j ){

        if(disps[j].c != compAtual){
          if( i%(2*g-1) > (g-2) ){ //segmento vertical
            if ( (disps[j].y - disps[j-1].y) <= 2*r){
              //std::cout << "\t" << disps[j].x << " " << (disps[j-1].y + disps[j].y ) /2 << std::endl;
              _x.push_back(disps[j].x);
              _y.push_back((disps[j-1].y + disps[j].y)/2);
              _c1.push_back(disps[j-1].c);
              _c2.push_back(disps[j].c);
              _c3.push_back(-2);
              _c4.push_back(-2);
            }
          }else{ //segmento horizontal
            if ( (disps[j].x - disps[j-1].x) <= 2*r){
              //std::cout << "\t" << (disps[j-1].x + disps[j].x ) /2 << " " << disps[j].y << std::endl;
              _x.push_back((disps[j-1].x + disps[j].x)/2);
              _y.push_back(disps[j].y);
              _c1.push_back(disps[j-1].c);
              _c2.push_back(disps[j].c);
              _c3.push_back(-2);
              _c4.push_back(-2);
            }
          }
          compAtual = disps[j].c;
        }
        //std::cout << disps[j].d << " " << disps[j].c << "\t" << disps[j].x << " " << disps[j].y << " " << std::endl;
      }
      delete [] disps;
    }
  }

  //nos cruzamentos
  for (int i =0; i < g*g; ++i){
    int borda; //1-OS,2-S,3-LS,4-O,5-L,6-NO,7-N,8-NL
    int N, L, O, S;
    int xCruz = i%g;
    int yCruz = floor(i/g);
    segmentosDoCruzamento(g,xCruz,yCruz,&N,&L,&O,&S,&borda);

    int Nd,Nc,Ld,Lc,Od,Oc,Sd,Sc;
    double Nx,Ny,Lx,Ly,Ox,Oy,Sx,Sy;
    maisProximoAoCruzamentoDoSegmento(g, 'N', N, borda, dispositivosNoSegmento,
      comp, Cx, Cy, &Nd, &Nc, &Nx, &Ny);
    //std::cout << Nd << " " << Nc << " " << Nx << " " << Ny << " " << std::endl;
    maisProximoAoCruzamentoDoSegmento(g, 'L', L, borda, dispositivosNoSegmento,
      comp, Cx, Cy, &Ld, &Lc, &Lx, &Ly);
    //std::cout << Ld << " " << Lc << " " << Lx << " " << Ly << " " << std::endl;
    maisProximoAoCruzamentoDoSegmento(g, 'O', O, borda, dispositivosNoSegmento,
      comp, Cx, Cy, &Od, &Oc, &Ox, &Oy);
    //std::cout << Od << " " << Oc << " " << Ox << " " << Oy << " " << std::endl;
    maisProximoAoCruzamentoDoSegmento(g, 'S', S, borda, dispositivosNoSegmento,
      comp, Cx, Cy, &Sd, &Sc, &Sx, &Sy);
    //std::cout << Sd << " " << Sc << " " << Sx << " " << Sy << " " << std::endl;


    std::list<int> componentesDosMaisProximos = {Nc,Lc,Oc,Sc};
    componentesDosMaisProximos.sort();
    componentesDosMaisProximos.unique();
    if (componentesDosMaisProximos.size() > 1 ){
    //dois a dois - nos cruzamentos
      bool NL=false,NO=false,NS=false,LO=false,LS=false,OS=false;
      double Bx,By;
      //std::cout << Nc << " " << Lc << " " << Oc << " " << Sc << " " << std::endl;

      if(Nc!=Lc && cabeUmAcesso(g, r, eps, xCruz, yCruz, Nx, Ny, Lx, Ly, &Bx, &By)){
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":NL\t" << Bx << " " << By << "\t" << Nc << " " << Lc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Nc);
          _c2.push_back(Lc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          NL=true;
        }
      }//NL

      if(Nc!=Oc && cabeUmAcesso(g, r, eps, xCruz, yCruz, Nx, Ny,
        xCruz+fabs(Ox-xCruz), Oy, &Bx, &By)){
        Bx=xCruz-fabs(Bx-xCruz);
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":NO\t" << Bx << " " << By << "\t" << Nc << " " << Oc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Nc);
          _c2.push_back(Oc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          NO=true;
        }
      }//NO

      if(Nc!=Sc && Ny-Sy<2*r){
        Bx=(Nx + Sx)/2;
        By=(Ny + Sy)/2;
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":NS\t" << Bx << " " << By << "\t" << Nc << " " << Sc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Nc);
          _c2.push_back(Sc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          NS=true;
        }
      }//NS

      if(Lc!=Oc && Lx-Ox<2*r){
        Bx=(Lx + Ox)/2;
        By=(Ly + Oy)/2;
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":LO\t" << Bx << " " << By << "\t" << Lc << " " << Oc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Lc);
          _c2.push_back(Oc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          LO=true;
        }
      }//LO

      if(Lc!=Sc && cabeUmAcesso(g, r, eps, xCruz, yCruz, Sx,
        yCruz+fabs(Sy-yCruz), Lx, Ly, &Bx, &By)){
        By=yCruz-fabs(By-yCruz);
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":LS\t" << Bx << " " << By << "\t" << Lc << " " << Sc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Lc);
          _c2.push_back(Sc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          LS=true;
        }
      }//LS

      if(Oc!=Sc && cabeUmAcesso(g, r, eps, xCruz, yCruz, Sx,
        yCruz+fabs(Sy-yCruz), xCruz+fabs(Ox-xCruz), Oy, &Bx, &By)){
        Bx=xCruz-fabs(Bx-xCruz);
        By=yCruz-fabs(By-yCruz);
        if(!(Bx>g-1||Bx<0)&&!(By>g-1||By<0)){
          //std::cout << i << ":OS\t" << Bx << " " << By << "\t" << Oc << " " << Sc << std::endl;
          _x.push_back(Bx);
          _y.push_back(By);
          _c1.push_back(Oc);
          _c2.push_back(Sc);
          _c3.push_back(-2);
          _c4.push_back(-2);
          OS=true;
        }
      }//OS


      if (componentesDosMaisProximos.size() > 2 ){
      //três a três - nos cruzamentos
        if(NL && NO && LO){
          if(cabeUmAcessoNLO(g, r, eps, xCruz, yCruz, Nx, Ny, Lx, Ly, Ox, Oy, &Bx, &By)){
            //std::cout << i << ":NLO\t" << Bx << " " << By << "\t" << Nc << " " << Lc << " " << Oc << std::endl;
            _x.push_back(Bx);
            _y.push_back(By);
            _c1.push_back(Nc);
            _c2.push_back(Lc);
            _c3.push_back(Oc);
            _c4.push_back(-3);
          }
        }
        if(NL && NS && LS){
          if(cabeUmAcessoNLS(g, r, eps, xCruz, yCruz, Nx, Ny, Lx, Ly, Sx, Sy, &Bx, &By)){
            //std::cout << i << ":NLS\t" << Bx << " " << By << "\t" << Nc << " " << Lc << " " << Sc << std::endl;
            _x.push_back(Bx);
            _y.push_back(By);
            _c1.push_back(Nc);
            _c2.push_back(Lc);
            _c3.push_back(Sc);
            _c4.push_back(-3);
          }
        }
        if(NO && NS && OS){
          if(cabeUmAcessoNOS(g, r, eps, xCruz, yCruz, Nx, Ny, Ox, Oy, Sx, Sy, &Bx, &By)){
            //std::cout << i << ":NOS\t" << Bx << " " << By << "\t" << Nc << " " << Oc << " " << Sc << std::endl;
            _x.push_back(Bx);
            _y.push_back(By);
            _c1.push_back(Nc);
            _c2.push_back(Oc);
            _c3.push_back(Sc);
            _c4.push_back(-3);
          }
        }
        if(LO && LS && OS){
          if(cabeUmAcessoLOS(g, r, eps, xCruz, yCruz, Lx, Ly, Ox, Oy, Sx, Sy, &Bx, &By)){
            //std::cout << i << ":LOS\t" << Bx << " " << By << "\t" << Lc << " " << Oc << " " << Sc << std::endl;
            _x.push_back(Bx);
            _y.push_back(By);
            _c1.push_back(Lc);
            _c2.push_back(Oc);
            _c3.push_back(Sc);
            _c4.push_back(-3);
          }
        }

        if (componentesDosMaisProximos.size() > 3 ){
        //quatro a quatro - nos cruzamentos
          if(NL && NO && NS && LO && LS && OS){
            if(cabeUmAcessoNLOS(g, r, eps, xCruz, yCruz, Nx, Ny, Lx, Ly, Ox, Oy, Sx, Sy, &Bx, &By)){
              //std::cout << i << ":NLOS\t" << Bx << " " << By << "\t" << Nc << " " << Lc << " " << Oc << " " << Sc << std::endl;
              _x.push_back(Bx);
              _y.push_back(By);
              _c1.push_back(Nc);
              _c2.push_back(Lc);
              _c3.push_back(Oc);
              _c4.push_back(Sc);
            }
          }
        }

      }

    }
  }

  //um a um - por componente
  _c.sort(); _c.unique();

  std::list<int>::const_iterator c1;
  for (c1 = _c.begin(); c1 != _c.end(); ++c1) {
    //std::cout << *c1 << " ";
    bool busqueUm=true; int d1=0;
    while (busqueUm && d1 < Dispositivos ){
      if ( comp[d1] == *c1 ) {
        _x.push_back(Cx[d1]);
        _y.push_back(Cy[d1]);
        _c1.push_back(*c1);
        _c2.push_back(-1);
        _c3.push_back(-1);
        _c4.push_back(-1);

        busqueUm = false;
      } else d1++;
    }
    //std::cout << d1 << " " << comp[d1] << std::endl;
  }


  delete [] componentesNoSegmento;
  delete [] dispositivosNoSegmento;

  return Rcpp::List::create(
    Rcpp::Named("x") = _x, Rcpp::Named("y") = _y,
    Rcpp::Named("c1") = _c1, Rcpp::Named("c2") = _c2, Rcpp::Named("c3") = _c3, Rcpp::Named("c4") = _c4
  );
}
