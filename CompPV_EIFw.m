function [P0,V]=CompPV_EIFw(para,VP)

    % This program computes the P(V)-distributon for the <w>-approach
    %
    % Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
    %           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV, 
    %           mean(Isyn)/pA, std(Isyn)/pA]
    %        2) VP = V-increment and V-regime of interest, VP(1)=bin,
    %                VP(2)=V:min
    %
    % Output: P(V) and V 

    C=para(1);
    gL=para(2);
    EL=para(3);
    sf=para(4);
    Vup=para(5);
    tcw=para(6);
    b=para(8);
    Vr=para(9);
    Vth=para(10);
    tm=C/gL;

    % define input regime
    m0=para(11)/gL;
    sig0=para(12)/sqrt(4*C*gL);

    % V-increment and V-regime of interest
    if isempty(VP)
        bin=0.1;
        V=Vr-20:bin:Vup;
    else
        bin=VP(1);
        V=VP(2):bin:Vup;
    end

    % define f0 
    OPTIONS=optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
    f0 = fminsearch(@SelfConstEqs_Q0,5,OPTIONS,para);
    avgW=b*tcw/1000*f0;
    m0=m0-avgW/gL;

    % compute P0 with EIF-FP equation
    P0=zeros(1,length(V));
    for i=1:length(V)
        P0(i)=integral(@(x) (2*tm/1000*f0/sig0.^2)*exp((2*(EL+m0)*(V(i)-x)+(x.^2-V(i).^2)-2*sf.^2*(exp((x-Vth)/sf)-exp((V(i)-Vth)/sf)))/sig0.^2),max(Vr,V(i)),Vup);
    end
end

%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

function q=SelfConstEqs_Q0(f,ModPar)

    f0_old=f;
    
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(ModPar(1:10));
    tm=Cm/gL/1000; tcw=tcw/1000;
    k=sqrt(tm/tcw);
    mu=ModPar(11)/gL;                   
    sig=ModPar(12)/sqrt(4*Cm*gL);       
    b=b/(k^2*gL);                       
    
    K=@(x) (2/sig^2)*(sf^2*exp((x-Vth)./sf) - (x.^2)/2 + (EL+mu-b*tm*f0_old)*x); 
    
    f0=sig^2./(2*tm*integral2(@(x,y) exp(K(x)-K(y)),Vr-50,Vup,@(x) max(x,Vr),Vup));

    if (f0>=0 && f0_old>=0)
        q = (f0_old-f0).^2;
    else
        q = 1e10;
    end

end

% (c) 2014 L. Hertaeg, D. Durstewitz and N. Brunel
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim