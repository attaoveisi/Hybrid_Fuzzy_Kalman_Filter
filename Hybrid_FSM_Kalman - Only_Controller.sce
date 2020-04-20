clc
clear

A=[-0.773996386153021,96.9795539479971,0,0,0,0;-96.4009703982576,-1.63547954896882,0,0,0,0;0,0,-2.34381015553514,467.827125148274,0,0;0,0,-467.843356672122,-2.35020309294106,0,0;0,0,0,0,-205.796226457431,502.302308896131;0,0,0,0,-502.497944160656,-205.977559293018];
dum1=size(A);
n=dum1(1,1);
B=[-0.230745598618469,-0.153904745996477;0.0612777559552695,0.0606173681667591;-25.0923163191213,13.5416143638357;4.01556833701487,-3.21559278869453;-4.90371370694887,20.5766878763535;-5.44248661253549,6.08611624326062];
dum2=size(B);
q=dum2(1,2);
C=[-23.1554338545303,9.56281796834669,0.435674513053965,-0.466194198179862,-0.575339270789856,-0.316369760481547];
dum3=size(C);
p=dum3(1,1);
H=[0.258916103898922;0.101956592759767;-1.80129126804772;-4.56185147008189;25.0926976512821;-19.7319673467244];
dum4=size(H);
s=dum4(1,2);
G1=H;
H1=ones(n,1)*1e-6;
dum5=size(H1);
r=dum5(1,2);

loadmatfile('Lk.mat');

loadmatfile('M_A.mat');
loadmatfile('N_A.mat');
loadmatfile('M_B.mat');
loadmatfile('N_B.mat');
loadmatfile('M_C.mat');
loadmatfile('N_C.mat');

dum6=sqrt(2e-2);
M_A=dum6*M_A;
M_B=dum6*M_B;
M_C=dum6*M_C;
N_A=dum6*N_A;
N_B=dum6*N_B;
N_C=dum6*N_C;

fm=1e-1;

//dum6=1e-3;
//M_A=[0.1,0,0,0.1,0,0]'*dum6;
//M_B=[0.1,0,0.1,0,0,0]'*dum6;
//M_C=0.1*dum6;
//N_A=[0,0,0.1,0,0,0.1];
//N_B=[0.1,0];
//N_C=[0,0,0.1,0,0,0.1];

Linf=C'*C*1e-3;
function [LME, LMI, OBJ]=HybridFSMKDRC(XLIST)
[X,Kh,eps6,eps7,eps8,eps9,gama_wc]= XLIST(:)
LME=list(X-X')
LMI=list(-([X*A'+A*X+B*Kh+Kh'*B'+eps6*fm^2*eye(n,n)+eps7*M_A*M_A'+eps8*M_B*M_B'+eps9*M_B*M_B'+Linf'*Linf,-B*Kh,G1,H1,X,X*N_A',Kh'*N_B',zeros(n,n);-Kh'*B',zeros(n,n+s+r+n+n+n),Kh'*N_B';G1',zeros(s,n),-gama_wc*eye(s,s),zeros(s,r+n+n+n+n);H1',zeros(r,n+s+r+n+n+n+n);X,zeros(n,n+s+r),-eps6*eye(n,n),zeros(n,n+n+n);N_A*X,zeros(n,n+s+r+n),-eps7*eye(n,n),zeros(n,n+n);N_B*Kh,zeros(n,n+s+r+n+n),-eps8*eye(n,n),zeros(n,n);zeros(n,n),N_B*Kh,zeros(n,s+r+n+n+n),-eps9*eye(n,n)]),X,eps6,eps7,eps8,eps9,-gama_wc+5e7)
OBJ=[]
endfunction

dum7=1e0;

X0=eye(n,n)*1e1;
Kh0=zeros(q,n);
eps6_0=1*dum7;
eps7_0=1*dum7;
eps8_0=1*dum7;
eps9_0=1*dum7;
gama_wc0=1e9;

Init_guess=list(X0,Kh0,eps6_0,eps7_0,eps8_0,eps9_0,gama_wc0);

Mbound=1e0;
abstol=5e-6;
nu=10;
maxiters=500;
reltol=1e-10;
options=[Mbound,abstol,nu,maxiters,reltol];

Ans_LMI=lmisolver(Init_guess,HybridFSMKDRC,options);
//Ans_LMI=lmisolver(Init_guess,HybridFSMKDRC);

X=Ans_LMI(1);
Kh=Ans_LMI(2);
eps6=Ans_LMI(3);
eps7=Ans_LMI(4);
eps8=Ans_LMI(5);
eps9=Ans_LMI(6);
gama_wc=Ans_LMI(7);
K=Kh*(X^-1);
