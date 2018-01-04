function fr0=Fr_AdExDist_SC(ModPar,f_guess)

    % This program computes the approximated firing rate (self-consistently)
    % of an AdEx neuron (with parameter a=0) by means of the distribution 
    % approach. 
    %
    % Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
    %           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV, 
    %           mean(Isyn)/pA, std(Isyn)/pA]
    %        2) f_guess = an initial guess for the firing rate in Hz
    %
    % Output: firing rate fr0/Hz      
    
    OPTIONS = optimset('TolFun',1e-2,'TolX',1e-2,'Display','off');
    fr0 = fminsearch(@SelfConstEqs_Q0_Dist,f_guess,OPTIONS,ModPar);
    
end

function q=SelfConstEqs_Q0_Dist(f,ModPar)

    f0_old=f;
    
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=(ModPar(1:10));
    tm=Cm/gL;
    m0=ModPar(11)/gL;
    sig0=ModPar(12)/sqrt(4*Cm*gL);
    ModPar1(1:5) = ModPar(1:5);      
    ModPar1(6) = ModPar(6)*ModPar(8)*f0_old/1000;
    ModPar1(7:8) = ModPar(9:10);
    
    if f0_old>0
        VarISI=CompVarISI_EIF(ModPar1,[m0 sig0]);
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
        
    fun= @(z) 1000./integral2(@(V,x) (2*tm/sig0.^2)*exp((2*(EL+m0-z/gL)*(V-x)+(x.^2-V.^2)-2*sf.^2*(exp((x-Vth)/sf)-exp((V-Vth)/sf)))/sig0.^2),lb,Vup,@(V) max(V,Vr),Vup,'AbsTol',1e-3,'RelTol',1e-3);
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