 function result = CompVarISI_EIF(ModPar,IPar)
 
% This program computes the variance of the ISI-distribution of an EIF
% neuron (for AdEx neuron: a=0 and w=b*tcw*f; average approach)
%
% Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
%           sf/mV, Vup/mV, w=b*tcw*f in pA, Vr/mV, Vth/mV]
%        2) IPar = vector with mean [IPar(1)] and std [IPar(2)] of Isyn 
%
% Output: variance of the ISI-distribution 

C=ModPar(1);
gL=ModPar(2);
tm=C/gL;
EL=ModPar(3);
sf=ModPar(4);
Vup=ModPar(5);
avgW=ModPar(6);
Vr=ModPar(7);
Vth=ModPar(8);

mu0=IPar(1);
mu0=mu0-avgW/gL;
sig0=IPar(2);
E=EL+mu0;

lb=Vr-20;
ATol=1e-4;
RTol=1e-2;

result = integral2(@InnerFunc,Vr,Vup,lb,@(V) V,'Method','tiled','AbsTol',ATol,'RelTol',RTol);

    function out = InnerFunc(V,u)  

        F=@(r) (2*E.*r - r.^2 + 2*sf^2.*exp((r-Vth)./sf))/sig0^2;
        I=arrayfun(@(u) integral(@(z) exp(F(z)-F(u)), lb, u, 'ArrayValued', true, 'AbsTol',ATol,'RelTol',RTol),u);

        out = (8*tm^2*exp(F(u)-F(V)).*I.^2)/sig0^4;
    end

 end

% (c) 2014 L. Hertaeg, D. Durstewitz and N. Brunel
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim