function result = R2Corr_EIF_Units(ModPar,IPar,flag)

% This program computes the approximated firing rate correction 2nd order
% of an AdEx neuron (with parameter a=0) in colored noise case.
%
% Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
%           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV]
%        2) IPar = vector with mean [IPar(1)] and std [IPar(2)] of Isyn 
%        3) flag = switch between <w> (1) and w-Dist approach (~=1)
%
% Output: firing rate correction 2nd order in Hz

gL=ModPar(2);
EL=ModPar(3);
sf=ModPar(4);
Vup=ModPar(5);
tcw=ModPar(6);     
b=ModPar(8);
Vr=ModPar(9);
Vth=ModPar(10);

mu0=IPar(1);
sig0=IPar(2);
r0=IPar(3);
avgW=tcw/1000*r0*b;
E=EL+mu0-avgW/gL;

F=@(r) (2/sig0^2) * (sf^2.*exp((r-Vth)./sf) - (r.^2)/2 + E*r);
df=@(r) exp((r-Vth)./sf) - 1;

if flag==1
    [P0,V]=CompPV_EIFw(ModPar,[]);
else
    [P0,V]=CompPV_EIF_Distw(ModPar,[]);
end
V=[V Vup]; P0=[P0 0];
Q0=@(z) interp1(V,P0,z);
dQ0=@(z) min(max(diff(P0)./mean(diff(V))),interp1(V(1:end-1),diff(P0)./mean(diff(V)),z));

M=Vr-20;
result=r0*integral2(@(V,x) df(x).*exp(F(V)-F(x)).*(Q0(x)./sf+dQ0(x)),M,Vup,@(V) V,Vup);


% (c) 2014 L. Hertaeg, D. Durstewitz and N. Brunel
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim