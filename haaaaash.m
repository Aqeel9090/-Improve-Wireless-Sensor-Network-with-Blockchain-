clc
clear all;
L=100;
% Xini=0.7;
Xini=rand;
alfa=4;
z=0,

%---------------hash generate-------------------------------------

x(1)=Xini;
for i=2:L
    x(i)=alfa*x(i-1)*(1-x(i-1));

    y(i)=(x(i)*100);
   z =[z,dec2bin(y(i))],

end

z (1:256)




