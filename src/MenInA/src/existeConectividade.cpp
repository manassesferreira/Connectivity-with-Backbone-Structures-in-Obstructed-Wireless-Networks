#include <math.h>
//#include <iostream>

void bordaPeriodica(int g, int borda, double *dx, double *dy, double d1x, double d1y, double d2x, double d2y){
  double d1xP=d1x, d1yP=d1y, d2xP=d2x, d2yP=d2y;
  switch ( borda ) {
    case 1: //OS
      if (d1x > 1 ) d1xP = d1x - (g-1);//O
      if (d2x > 1 ) d2xP = d2x - (g-1);//O
      if (d1y > 1 ) d1yP = d1y - (g-1);//S
      if (d2y > 1 ) d2yP = d2y - (g-1);//S
      break;
    case 2: //S
      if (d1y > 1 ) d1yP = d1y - (g-1);//S
      if (d2y > 1 ) d2yP = d2y - (g-1);//S
      break;
    case 3: //LS
      if (d1x <= 1 ) d1xP = d1x + (g-1);//L
      if (d2x <= 1 ) d2xP = d2x + (g-1);//L
      if (d1y > 1 ) d1yP = d1y - (g-1);//S
      if (d2y > 1 ) d2yP = d2y - (g-1);//S
      break;
    case 4: //O
      if (d1x > 1 ) d1xP = d1x - (g-1);//O
      if (d2x > 1 ) d2xP = d2x - (g-1);//O
      break;
    case 5: //L
      if (d1x <= 1 ) d1xP = d1x + (g-1);//L
      if (d2x <= 1 ) d2xP = d2x + (g-1);//L
      break;
    case 6: //NO
      if (d1y <= 1 ) d1yP = d1y + (g-1);//N
      if (d2y <= 1 ) d2yP = d2y + (g-1);//N
      if (d1x > 1 ) d1xP = d1x - (g-1);//O
      if (d2x > 1 ) d2xP = d2x - (g-1);//O
      break;
    case 7: //N
      if (d1y <= 1 ) d1yP = d1y + (g-1);//N
      if (d2y <= 1 ) d2yP = d2y + (g-1);//N
      break;
    case 8: //NL
      if (d1y <= 1 ) d1yP = d1y + (g-1);//N
      if (d2y <= 1 ) d2yP = d2y + (g-1);//N
      if (d1x <= 1 ) d1xP = d1x + (g-1);//L
      if (d2x <= 1 ) d2xP = d2x + (g-1);//L
      break;
  }
  *dx = fabs(d1xP - d2xP);
  *dy = fabs(d1yP - d2yP) ;
}

int existeConectividade(int g, int borda, double r, double eps, double d1x, double d1y, double d2x, double d2y){
  double dx, dy, wid = 2*eps;

  if (borda > 0) {
    bordaPeriodica(g, borda, &dx, &dy, d1x, d1y, d2x, d2y);
  }else{
    dx = fabs(d1x - d2x);
    dy = fabs(d1y - d2y);
  }

  //std::cout << "\t\t" << dx << " " << dy << std::endl;

  if( dx == 0 ){
    if ( dy > r ) return 0;
    else return 1;
  } else if( dy == 0 ){
    if ( dx > r ) return 0;
    else return 1;
  } else {
    if ( dx*dx + dy*dy > r*r ) return 0;
    else {
      if( ( dx < eps || dy < eps ) || ( dx < wid && dy < wid ) ||
          ( dx < wid && fabs(dy*(1-eps/dx)) < eps ) ||
          ( dy < wid && fabs(dx*(1-eps/dy)) < eps )
        ) return 1;
      else return 0;
    }
  }
}
