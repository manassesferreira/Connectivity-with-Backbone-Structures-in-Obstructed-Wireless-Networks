//usamos o contorno períodico na função abaixo para evitar efeitos de borda
//como por exemplo, cruzamentos com dois/três segmentos, ao invés, de quatro
//código para o tipo de borda 1-OS,2-S,3-LS,4-O,5-L,6-NO,7-N,8-NL
void segmentosDoCruzamento(int g, int xCruz, int yCruz, int *N, int *L, int *O, int *S, int *borda){
  int LCruz = (2*g-1)*yCruz + xCruz;
  int NCruz = LCruz + (g-1);
  int OCruz = LCruz - 1;
  int SCruz = LCruz - g;

  if(xCruz==0){
    OCruz=LCruz + (g-2);
    if(yCruz==0){ *borda=1;
      SCruz=(2*g-1)*(g-1)-g + xCruz;
    }else if(yCruz==(g-1)){ *borda=6;
      NCruz=xCruz + (g-1);
    }else{ *borda=4;
      //nothing
    }
  }else if(xCruz==(g-1)){
    LCruz=LCruz - (g-1);
    if(yCruz==0){ *borda=3;
      SCruz=(2*g-1)*(g-1)-g + xCruz;
    }else if(yCruz==(g-1)){  *borda=8;
      NCruz=xCruz + (g-1);
    }else{  *borda=5;
      //else
    }
  }else{
    if(yCruz==0){ *borda=2;
      SCruz=(2*g-1)*(g-1)-g + xCruz;
    }else if(yCruz==(g-1)){ *borda=7;
      NCruz=xCruz + (g-1);
    }else{ *borda=0;
      //matters
    }
  }
  //std::cout << NCruz << " " << LCruz << " " << OCruz << " " << SCruz << " " << std::endl;
  *N=NCruz; *L=LCruz; *O=OCruz; *S=SCruz;
}

