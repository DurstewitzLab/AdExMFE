function [P0,V]=CompPV_EIF_Distw(para,VP)

    % This program computes the P(V)-distributon for the w-Dist.-approach
    %
    % Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV,
    %           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV,
    %           mean(Isyn)/pA, std(Isyn)/pA]
    %        2) VP = V-increment and V-regime of interest, VP(1)=bin,
    %                VP(2)=V:min
    %
    % Output: P(V) and V

    % meaningful parameter names
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

    % define f_distr    
    OPTIONS = optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
    fr0 = fminsearch(@SelfConstEqs_Q0_Dist,f_est,OPTIONS,para);
    
    % hier definiere die overall Verteilung (also unabhängig von V) fuer w
    ModPar1(1:5) = para(1:5);      
    ModPar1(6) = para(6)*para(8)*fr0/1000;
    ModPar1(7:8) = para(9:10);
    
    if fr0>0
        VarISI=CompVarISI_EIF(ModPar1,[para(11)/gL para(12)/sqrt(4*C*gL)]);
        b1=(tcw/(fr0*VarISI/1000+tcw))^(1/((fr0/1000)^2*VarISI));
        if b1==1
            varW=0;
        else
            varW=b^2*tcw*fr0/2000 *((1+b1)/(1-b1) - 2*tcw*fr0/1000);
        end
    else
        varW=0;
    end
    
    wmax=b/(1-exp(-1000/(fr0*tcw)));
    wmin=max(0.001,b*exp(-1000/(fr0*tcw))/(1-exp(-1000/(fr0*tcw))));
    avgW=para(6)*para(8)*fr0/1000;
    Kn=gamcdf(wmax,avgW^2/varW,varW/avgW)-gamcdf(wmin,avgW^2/varW,varW/avgW); 
    Fpdf=@(z) gampdf(z,avgW^2/varW,varW/avgW)/Kn; 
    
    % compute f(w)
    z=wmin:(wmax-wmin)/18:wmax;
    out=zeros(1,length(z));
    for i=1:length(z)
        K=@(x) (2/sig0^2)*(sf^2*exp((x-Vth)./sf) - (x.^2)/2 + (EL+m0-z(i)/gL)*x); 
        out(i)=sig0^2./(2*tm/1000*integral2(@(x,y) exp(K(x)-K(y)),Vr-50,Vup,@(x) max(x,Vr),Vup));
    end

    % compute P with EIF-FP equation
    P0=zeros(1,length(V));
    for i=1:length(V)
        P0(i) = integral2(@(x,w) Fpdf(w).*(2*tm/1000.*interp1(z,out,w,'pchip')./sig0.^2).*exp((2*(EL+m0-w./gL).*(V(i)-x)+(x.^2-V(i)^2)-2*sf.^2*(exp((x-Vth)./sf)-exp((V(i)-Vth)/sf)))/sig0.^2),max(Vr,V(i)),Vup,wmin,wmax);
    end
end

%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%

function q=SelfConstEqs_Q0_Dist(f,ModPar)

    f0_old=f;
    
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(ModPar(1:10));
    tm=Cm/gL;
    m0=ModPar(11)/gL;
    sig0=ModPar(12)/sqrt(4*Cm*gL);
    ModPar1(1:5) = ModPar(1:5);      
    ModPar1(6) = ModPar(6)*ModPar(8)*f0_old/1000;
    ModPar1(7:8) = ModPar(9:10);
    
    if f0_old>0
        VarISI=CompVarISI_EIF(ModPar1,[ModPar(11)/gL ModPar(12)/sqrt(4*Cm*gL)]);
        b1=(tcw/(f0_old*VarISI/1000+tcw))^(1/((f0_old/1000)^2*VarISI));
        if b1==1
            varW=0;
        else
            varW=b^2*tcw*f0_old/2000 *((1+b1)/(1-b1) - 2*tcw*f0_old/1000);
        end
    else
        varW=0;
    end
    
    wmax=b/(1-exp(-1000/(f0_old*tcw)));
    wmin=max(0.001,b*exp(-1000/(f0_old*tcw))/(1-exp(-1000/(f0_old*tcw))));
    lb=Vr-100;
    avgW=ModPar(6)*ModPar(8)*f0_old/1000;
    K=gamcdf(wmax,avgW^2/varW,varW/avgW)-gamcdf(wmin,avgW^2/varW,varW/avgW); 
    Fpdf=@(z) gampdf(z,avgW^2/varW,varW/avgW)/K; 
        
    fun= @(z) 1000./integral2(@(V,x) (2*tm/sig0.^2)*exp((2*(EL+m0-z/gL)*(V-x)+(x.^2-V.^2)-2*sf.^2*(exp((x-Vth)/sf)-exp((V-Vth)/sf)))/sig0.^2),lb,Vup,@(V) max(V,Vr),Vup);
    f0=integral(@(z) Fpdf(z).*fun(z),wmin,wmax,'ArrayValued',true);
        
    if (f0>=0 && f0_old>=0)
        q = (f0_old-f0).^2;
    else
        q = 1e10;
    end

end

% (c) 2014 L. Hertaeg, D. Durstewitz and N. Brunel
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim