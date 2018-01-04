function f = FR_AdEx_2DFPeq_LCTL(ModPar,f_guess)

    % This function computes the firing rate derived from the zeroth and second
    % order contribution to the solution of the Fokker-Planck equation in the 
    % long correlation time limit.
    %
    % Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
    %           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV, 
    %           mean(Isyn)/pA, std(Isyn)/pA]
    %        2) f_guess = an initial guess for the firing rate in Hz
    %
    % Output: firing rate vector f=[zeroth order, second order correction,
    %         both together]

    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(ModPar(1:10));
    tm=Cm/gL/1000; tcw=tcw/1000;
    k=sqrt(tm/tcw);
    mu=ModPar(11)/gL;                   
    sig=ModPar(12)/sqrt(4*Cm*gL);       
    b=b/(k^2*gL);                       
    
    ATol=1e-2; 
    RTol=1e-1;
    Vlb=Vr-20;

    OPTIONS=optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
    fr0 = fminsearch(@SelfConstEqs_Q0,f_guess,OPTIONS,ModPar);

    if fr0>1e-8
    
        K=@(x) (2/sig^2)*(sf^2*exp((x-Vth)./sf) - (x.^2)/2 + (EL+mu-b*tm*fr0)*x);
        matlabpool(8)

        parfor i=1:4
            if i==1
                I(i)=integral3(@(V,y,z) exp(K(V)-K(z)),Vlb,Vup,@(V) V,Vup,@(V,y) max(y,Vr),Vup,'AbsTol',ATol,'RelTol',RTol,'Method','tiled'); 
            elseif i==2
                I(i)=integral(@(V) integral3(@(x,y,z) exp(K(V)-K(x)+K(y)-K(z)), V,Vup,Vlb,@(x) x,@(x,y) max(y,Vr),Vup,'AbsTol',ATol,'RelTol',RTol,'Method','tiled'),Vlb,Vup,'ArrayValued',true,'AbsTol',ATol,'RelTol',RTol); 
            elseif i==3
                I(i)=integral(@(V) integral3(@(x,y,z) exp(K(V)-K(z)), V,Vup,@(x) x,Vup,@(x,y) max(y,Vr),Vup,'AbsTol',ATol,'RelTol',RTol,'Method','tiled'),Vlb,Vup,'ArrayValued',true,'AbsTol',ATol,'RelTol',RTol);
            elseif i==4
                I(i)=integral(@(V) integral(@(a) integral3(@(b,c,d) exp(K(V)-K(b)+K(c)-K(d)),a,Vup,Vlb,@(b) b,@(b,c) max(c,Vr),Vup,'AbsTol',ATol,'RelTol',RTol,'Method','tiled'),V,Vup,'ArrayValued',true,'AbsTol',ATol,'RelTol',RTol),Vlb,Vup,'ArrayValued',true,'AbsTol',ATol,'RelTol',RTol); 
            end
        end

        matlabpool('close')
        
        KJQ=4*I(1)/sig^4; KJKQ=4*I(2)/sig^4; KJJQ=8*I(3)/sig^6; KJJKQ=8*I(4)/sig^6;
        s2_z0 = (fr0*tm*KJKQ - 1/(2*fr0*tm))/(1/(fr0*tm) + b*tm*fr0*KJQ);
        fr2_N = KJKQ - (KJQ*b*s2_z0)^2/(KJKQ - 1/(2*(fr0*tm)^2)) - b*s2_z0^2*KJQ/(fr0^2*tm^2*KJKQ - 1/2) - 1/(sqrt(2)*fr0*tm)^2;
        fr2_Z = b*s2_z0*fr0*tm*KJJKQ - b^2*s2_z0^2*fr0*tm*KJJQ + KJQ*b*s2_z0^3/(KJKQ - 1/(2*(fr0*tm)^2))*(b*KJQ - KJKQ/s2_z0 + b^2*fr0*tm*KJJQ - b*fr0*tm*KJJKQ/s2_z0);
        fr2 = fr2_Z/fr2_N;

    else
        fr2=0;
    end

    fr = fr0 + k^2 * fr2/tm;
    f = [fr0 fr2 fr];

end

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