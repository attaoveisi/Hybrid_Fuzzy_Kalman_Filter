clc
clear

loadmatfile('A.mat');
loadmatfile('B.mat');
loadmatfile('C.mat');
loadmatfile('H.mat');

B=B/1e3;

dum1=size(A);
n=dum1(1,1);
dum2=size(B);
q=dum2(1,2);
dum3=size(C);
p=dum3(1,1);
dum4=size(H);
s=dum4(1,2);
G1=H*1;
H1=ones(n,1)*1e-3;
dum5=size(H1);
r=dum5(1,2);

loadmatfile('Lk.mat');

loadmatfile('M_A.mat');
loadmatfile('N_A.mat');
loadmatfile('M_B.mat');
loadmatfile('N_B.mat');
loadmatfile('M_C.mat');
loadmatfile('N_C.mat');

dum6=sqrt(1);
M_A=dum6*M_A;
M_B=dum6*M_B;
M_C=dum6*M_C;
N_A=dum6*N_A;
N_B=dum6*N_B;
N_C=dum6*N_C;

fm=1e-8;

//dum6=1e-3;
//M_A=[0.1,0,0,0.1,0,0]'*dum6;
//M_B=[0.1,0,0.1,0,0,0]'*dum6;
//M_C=0.1*dum6;
//N_A=[0,0,0.1,0,0,0.1];
//N_B=[0.1,0];
//N_C=[0,0,0.1,0,0,0.1];

Linf=C'*C;
function [LME, LMI, OBJ]=HybridFSMKDRC(XLIST)
[X,Kh,Ls,eps1,eps2,eps3,eps4,eps5,gama_we,eps6,eps7,eps8,eps9,gama_wc]= XLIST(:)
LME=list(X-X')
LMI=list(-([X*A'+A*X-X*C'*Lk'-Lk*C*X+eps1*M_B*M_B'+eps2*M_A*M_A'+eps3*Lk*M_C*M_C'*Lk'+eps4*M_B*M_B'+eps5*fm^2*eye(n,n)+X,zeros(n,n),G1,H1,-Ls,Kh'*N_B',zeros(n,n+n+n+n);zeros(n,n+n+s+r+p+n),X*N_A',X*N_C',Kh'*N_B',X;G1',zeros(s,n),-gama_we*eye(s,s),zeros(s,r+p+n+n+n+n+n);H1',zeros(r,n+s+r+p+n+n+n+n+n);-Ls',zeros(p,n+s+r+p+n+n+n+n+n);N_B*Kh,zeros(n,n+s+r+p),-eps1*eye(n,n),zeros(n,n+n+n+n);zeros(n,n),N_A*X,zeros(n,s+r+p+n),-eps2*eye(n,n),zeros(n,n+n+n);zeros(n,n),N_C*X,zeros(n,s+r+p+n+n),-eps3*eye(n,n),zeros(n,n+n);zeros(n,n),N_B*Kh,zeros(n,s+r+p+n+n+n),-eps4*eye(n,n),zeros(n,n);zeros(n,n),X,zeros(n,s+r+p+n+n+n+n),-eps5*eye(n,n)]),-([X*A'+A*X+B*Kh+Kh'*B'+eps6*fm^2*eye(n,n)+eps7*M_A*M_A'+eps8*M_B*M_B'+eps9*M_B*M_B'+Linf'*Linf,-B*Kh,G1,H1,X,X*N_A',Kh'*N_B',zeros(n,n);-Kh'*B',zeros(n,n+s+r+n+n+n),Kh'*N_B';G1',zeros(s,n),-gama_wc*eye(s,s),zeros(s,r+n+n+n+n);H1',zeros(r,n+s+r+n+n+n+n);X,zeros(n,n+s+r),-eps6*eye(n,n),zeros(n,n+n+n);N_A*X,zeros(n,n+s+r+n),-eps7*eye(n,n),zeros(n,n+n);N_B*Kh,zeros(n,n+s+r+n+n),-eps8*eye(n,n),zeros(n,n);zeros(n,n),N_B*Kh,zeros(n,s+r+n+n+n),-eps9*eye(n,n)]),X,eps1,eps2,eps3,eps4,eps5,gama_we,eps6,eps7,eps8,eps9,gama_wc)
OBJ=[]
endfunction

dum7=1e0;

X0=eye(n,n)*1e3;
Kh0=zeros(q,n);
Ls0=zeros(n,p);
eps1_0=1*dum7;
eps2_0=1*dum7;
eps3_0=1*dum7;
eps4_0=1*dum7;
eps5_0=1*dum7;
gama_we0=1e10;
eps6_0=1*dum7;
eps7_0=1*dum7;
eps8_0=1*dum7;
eps9_0=1*dum7;
gama_wc0=1e10;

Init_guess=list(X0,Kh0,Ls0,eps1_0,eps2_0,eps3_0,eps4_0,eps5_0,gama_we0,eps6_0,eps7_0,eps8_0,eps9_0,gama_wc0);

Mbound=1e5;
abstol=5e-6;
nu=2;
maxiters=5000;
reltol=1e-10;
options=[Mbound,abstol,nu,maxiters,reltol];

Ans_LMI=lmisolver(Init_guess,HybridFSMKDRC,options);
//Ans_LMI=lmisolver(Init_guess,HybridFSMKDRC);

X=Ans_LMI(1);
Kh=Ans_LMI(2);
Ls=Ans_LMI(3);
eps1=Ans_LMI(4);
eps2=Ans_LMI(5);
eps3=Ans_LMI(6);
eps4=Ans_LMI(7);
eps5=Ans_LMI(8);
gama_we=Ans_LMI(9);
eps6=Ans_LMI(10);
eps7=Ans_LMI(11);
eps8=Ans_LMI(12);
eps9=Ans_LMI(13);
gama_wc=Ans_LMI(14);
K=Kh*(X^-1);
