%Unidades en el SI, excepto el diámetro D(mm)

h0=15,ng=100,k=1.5,x=0.54,CVm=0.05,le=0.25,s=0.5,D=13.6,I0=0.02

vectorUnos=ones(ng,1);
matrizAcum=triu(ones(ng),0);

%Distancia al origen y cota de cada emisor
xR=transpose(matrizAcum)*(s.*vectorUnos);
zR=-I0.*xR;

%Distribución de presión para comenzar a calcular
h=h0.*vectorUnos;
hant=0.*vectorUnos;

%VariaciÃ³n manufactura
varManuf=1+CVm.*randn(ng,1);%Caudal de cada gotero/caudal meedio

  %CÃ¡lculo iterativo de la distribuciÃ³n de presiones en el ramal.
  while max(abs(h-hant))>1e-3;
    hant=h;
    q=(k.*h.^x).*varManuf;
    q=q.*(q>0);
    %Las pÃ©rdidas de carga se calculan con Blasius, viscosidad a 20ÂºC y longitud equivalente en las inserciones de ramal
    h=h0.*vectorUnos-zR-transpose(matrizAcum)*(0.465.*(matrizAcum*q).^1.75.*D.^-4.75.*(le+s));
  end