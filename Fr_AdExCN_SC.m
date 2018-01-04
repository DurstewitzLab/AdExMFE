function [fr_SigRed,fr_R2]=Fr_AdExCN_SC(ModPar,f_guess,flag_1,flag_2)

    % This program computes the approximated firing rate (self-consistently)
    % of an AdEx neuron (with parameter a=0) in the colored noise case
    %
    % Input: 1) ModPar = vector with all parameters [C/pF, gL/nS, EL/mV, 
    %           sf/mV, Vup/mV, tcw/ms, a=0, b/pA, Vr/mV, Vth/mV, 
    %           mean(Isyn)/pA, std(Isyn)/pA, ts/ms]
    %        2) f_guess = an initial guess for the firing rate in Hz
    %        3) flag_1 = variance reduction method (0: both, 1: <w>, 2:
    %           w-Dist)
    %        4) flag_2 = ts-correction of FP equation (1: <w>, 2: w-Dist,
    %           3: all together)
    %
    % Output: fr_SigRed = firing rate for <w> [fr_SigRed(1)] and w-Dist.
    %                     [fr_SigRed(2)]
    %         fr_R2 = firing rate computed by means of corrections due to 
    %                 synaptic time constant ts unequal zero [<w>, w-Distr.
    %                 over fr2, w-Distr. over fr0 + k^2 * fr2]
   
    fr_SigRed=zeros(1,2);
    fr_R2=zeros(1,3);
    
    if flag_1==1 
        % fr_SigRed: <w>
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_SigRed(1) = fminsearch(@SelfConstEqs_func1,f_guess,OPTIONS,ModPar);
    elseif flag_1==2
        % fr_SigRed: w-Distr
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_SigRed(2) = fminsearch(@SelfConstEqs_func2,f_guess,OPTIONS,ModPar);
    elseif flag_1==0
        % fr_SigRed: <w>
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_SigRed(1) = fminsearch(@SelfConstEqs_func1,f_guess,OPTIONS,ModPar);
        
        % fr_SigRed: w-Distr
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_SigRed(2) = fminsearch(@SelfConstEqs_func2,fr_SigRed(1),OPTIONS,ModPar);
    else
        % nothing
    end
    
    if flag_2==1
        % fr_R2: <w>
        fr_R2(1)=Func3(ModPar,f_guess);
    elseif flag_2==2
        % fr_R2: w-Distr. over fr2
        fr_R2(2)=Func4(ModPar,f_guess);

        % fr_R2: w-Distr. over fr0 + k^2 * fr2
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_R2(3) = fminsearch(@SelfConstEqs_func5,fr_R2(2),OPTIONS,ModPar);
    elseif flag_2==3
        % fr_R2: <w>
        fr_R2(1)=Func3(ModPar,f_guess);
        
        % fr_R2: w-Distr. over fr2
        fr_R2(2)=Func4(ModPar,f_guess);

        % fr_R2: w-Distr. over fr0 + k^2 * fr2
        OPTIONS=optimset('TolFun',1e-3,'TolX',1e-3,'Display','off');
        fr_R2(3) = fminsearch(@SelfConstEqs_func5,f_guess,OPTIONS,ModPar);
    else
        % nothing
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  firing rates computed by EIF-<w> and EIF-F(w) approach + std-reduction
function q=SelfConstEqs_func1(f,ModPar)

    f0_old=f;
    
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(ModPar(1:10));
    tm=Cm/gL/1000; tcw=tcw/1000; ModPar(13)=ModPar(13)/1000;
    mu=ModPar(11)/gL;                                              
    sig=sqrt(ModPar(12)^2/(1+ModPar(13)/tm))/sqrt(4*Cm*gL);     
    
    K=@(x) (2/sig^2)*(sf^2*exp((x-Vth)./sf) - (x.^2)/2 + (EL+mu-b*tcw*f0_old/gL)*x); 
    
    f0=sig^2./(2*tm*integral2(@(x,y) exp(K(x)-K(y)),Vr-50,Vup,@(x) max(x,Vr),Vup));

    if (f0>=0 && f0_old>=0)
        q = (f0_old-f0).^2;
    else
        q = 1e10;
    end

end

function q=SelfConstEqs_func2(f,ModPar)

    f0_old=f;
    
    [Cm,gL,EL,sf,Vup,tcw,~,b,Vr,Vth]=names(ModPar(1:10));
    tm=Cm/gL;
    m0=ModPar(11)/gL;
    sig0=sqrt(ModPar(12)^2/(1+ModPar(13)/tm))/sqrt(4*Cm*gL);
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

% R2- corrections
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

function f0=Func3(ModPar,f_guess)
    
    [Cm,gL,~,~,~,~,~,~,~,~]=names(ModPar(1:10));
        
    sig0=ModPar(12)/sqrt(4*Cm*gL);
    mu0=ModPar(11)/gL;

    OPTIONS=optimset('TolFun',1e-8,'TolX',1e-8,'Display','off');
    f02 = fminsearch(@SelfConstEqs_Q0,f_guess,OPTIONS,ModPar);
    
    R2=R2Corr_EIF_Units(ModPar,[mu0 sig0 f02],1);
    f0=f02+(ModPar(13)*gL/Cm/2)*R2;

end


function f0=Func4(ModPar,f_guess)
   
    [Cm,gL,~,~,~,~,~,~,~,~]=names(ModPar(1:10));
    
    f02=Fr_AdExDist_SC(ModPar,f_guess);
        
    sig0=ModPar(12)/sqrt(4*Cm*gL);
    mu0=ModPar(11)/gL;
    R2=R2Corr_EIF_Units(ModPar,[mu0 sig0 f02],2);
    f0=f02+(ModPar(13)*gL/Cm/2)*R2;

end


function q=SelfConstEqs_func5(f,ModPar)

    f0_old=f;
    
    [Cm,gL,~,~,~,tcw,~,b,~,~]=names(ModPar(1:10));
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
    wmin=b*exp(-1000/(f0_old*tcw))/(1-exp(-1000/(f0_old*tcw)));
    
    f0=f0_EIFwDistr_ColoredNoise_R2Corr(ModPar,ModPar1,wmin,wmax,varW,ModPar(13));

    if (f0>=0 && f0_old>=0)
        q = (f0_old-f0).^2;
    else
        q = 1e10;
    end

end

function f0=f0_EIFwDistr_ColoredNoise_R2Corr(para0,para,wmin,wmax,varW,tc_syn)

    avgW=para(6);
    mu0=para0(11)/para(2); 
    sig0=para0(12)/sqrt(4*para(1)*para(2));

    if (wmax>0 && varW>0 && avgW>0)

        wi=wmin:(wmax-wmin)/15:wmax;
        for ii=1:length(wi)
            para(6)=wi(ii);
            fi(ii)=f0_EIFavgW0(para0,avgW);
            if fi(ii)>0
                R2=R2Corr_EIF_Units(para0,[mu0 sig0 fi(ii)],1);
                fi(ii)=fi(ii)+(tc_syn*para(2)/para(1)/2)*R2;
            end
        end

        w = wmin+1e-10:(wmax-wmin)/500:wmax;
        p = polyfit(wi,fi,10);
        fw = polyval(p,w);

        pdf_w = gampdf(w,avgW^2/varW,varW/avgW); 
        K = gamcdf(wmax,avgW^2/varW,varW/avgW)-gamcdf(wmin,avgW^2/varW,varW/avgW);
        C = trapz(w,pdf_w); % normalization due to numerical issues

        if wmin==0
            pdf_w = pdf_w*K/C;
            p0 = 1-K;
            para(6)=0;
            
            f_w0=f0_EIFavgW0(para0,avgW);
            if f_w0>0
                R2=R2Corr_EIF_Units(para0,[mu0 sig0 f_w0],1);
                f_w0=f_w0+(tc_syn*para(2)/para(1)/2)*R2;
            end
            f0 = p0*f_w0+trapz(w,pdf_w.*fw);
        else
            pdf_w = pdf_w/C;
            f0 = trapz(w,pdf_w.*fw);
        end

    else
        para(6)=0;
        f0 = f0_EIFavgW0(para0,avgW);
        if f0>0
            R2=R2Corr_EIF_Units(para0,[mu0 sig0 f0],1);
            f0=f0+(tc_syn*para(2)/para(1)/2)*R2;
        end
    end
end

function f=f0_EIFavgW0(para,avgW)

    [C,gL,EL,sf,Vup,~,~,~,Vr,Vth]=names(para(1:10));
    tm=C/gL;

    m0=para(11)/gL-avgW/gL;
    sig0=para(12)/sqrt(4*C*gL);

    lb=Vr-20;

    t1=quad2d(@(V,x) (2*tm/sig0.^2)*exp((2*(EL+m0)*(V-x)+(x.^2-V.^2)-2*sf.^2*(exp((x-Vth)/sf)-exp((V-Vth)/sf)))/sig0.^2),lb,Vr,Vr,Vup);
    t2=quad2d(@(V,x) (2*tm/sig0.^2)*exp((2*(EL+m0)*(V-x)+(x.^2-V.^2)-2*sf.^2*(exp((x-Vth)/sf)-exp((V-Vth)/sf)))/sig0.^2),Vr,Vup,@(V) V,Vup);

    f=1000./(t1+t2);

end

% (c) 2014 L. Hertaeg, D. Durstewitz and N. Brunel
% Central Institute of Mental Health, Mannheim University of Heidelberg 
% and BCCN Heidelberg-Mannheim