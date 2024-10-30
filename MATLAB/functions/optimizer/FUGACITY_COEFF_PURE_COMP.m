
%function [PHI,Comp_factor_Z]=FUGACITY_COEFF_PURE_COMP(TYPE,PHASE,ALFA_function,T,P,Tc, Pc,acentric,ALPHA_constants)
function [PHI,Comp_factor_Z]=FUGACITY_COEFF_PURE_COMP(TYPE,PHASE,T,P,a,b,VECTOR_PROPERTIES_PURE_COMP)

% pure comp props
Tc=VECTOR_PROPERTIES_PURE_COMP(1);
Pc=VECTOR_PROPERTIES_PURE_COMP(2);


% R in kPa.m3/kmol.K
R=8.314;
N=1; % total moles
%guess P

%TYPE: SRK or PR
% Select parameters for the cubic EOS: PR=Peng-Robinson,% SRK=Sove-Redlick_kwong
% calculation of a critical and b
if TYPE == string('PR')
    delta_1=1+sqrt(2);
    delta_2=1-sqrt(2);
elseif TYPE == string('SRK')||TYPE == "RK"
    delta_1=1;
    delta_2=0;    
end
%{
%{_____________________________________________________________________________
% calculation of alpha function
if TYPE == string('PR')
    
    if ALFA_function == string ('SOAVE')
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == string ('TWU')
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
elseif TYPE == string('SRK')
    
    if ALFA_function == string ('SOAVE')
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == string ('TWU')
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
end
%calculation of a=ac*alpha
a=alpha.*ac;
%}
%_______________________________________________________________________________
% Calculation of the Cubic EOS Roots 
AA=(P*a)/((R*T)^2);
BB=b*P/(R*T);
Poly_EOS=[1, (BB*(delta_1+delta_2-1)-1), ((BB^2)*(delta_1*delta_2-delta_2-delta_1)-BB*(delta_2+delta_1)+AA),(-(BB^3)*delta_1*delta_2-(BB^2)*delta_1*delta_2-AA*BB)];

%running on GPU
% gpu_Poly_EOS=gpuArray(Poly_EOS);
% roots_EOS=transpose(roots(gpu_Poly_EOS));

roots_EOS=transpose(roots(Poly_EOS));

%chose only real roots and sort roots from max to min
%real_roots_EOS=roots_EOS(imag(roots_EOS)==0);
real_roots_EOS=real(roots_EOS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
control_img_roots=imag(roots_EOS);

if T==Tc && P==Pc
    Zmean=mean(real(roots_EOS));
    %Zmean=real_roots_EOS; % Twaha Mohamed  wen Jun 22 The value of the compressibility factor at the critical point should be the average of all real roots
    real_roots_EOS=Zmean;    
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sort roots from max to min
real_roots_non_negative=real_roots_EOS(real_roots_EOS>0);
real_roots_EOS=sort(real_roots_non_negative,'descend');

%check if Z>bP/RT
real_roots_EOS=real_roots_EOS(real_roots_EOS>(b*P/(R*T)));
end

% case 1: only one root
if length(real_roots_EOS)==1
Z=real_roots_EOS;

%case 2: only 2 roots are different
elseif length(real_roots_EOS)==2 
Zl=min(real_roots_EOS); %liquid
Zh=max(real_roots_EOS); %vapor
Z=[Zl,Zh];
%case 3: 3 real roots all different
elseif length(real_roots_EOS)==3
Zl=min(real_roots_EOS); %liquid
Zh=max(real_roots_EOS); % vapor
Z=[Zl,Zh];
end

%__________________________________________________________________________

%chose the phase under analysis: liq or vapor
%if length(Z)=2, then, system is splitted into two phases 
if length(Z)~=1
    
    if PHASE==string('LIQUID')
        Zliq=Zl;
        % Calculation of fugacity coefficient
        V=(Zliq*R*T)/P;
        ar_RT=-(log(1-(b/V)))-((a/(R*T*b*(delta_1-delta_2)))*log((1+(delta_1*b/V))/(1+(delta_2*b/V))));
        PHI=exp(ar_RT+Zliq-1-log(Zliq));
        Comp_factor_Z=Zl;
    elseif PHASE==string('VAPOR')
        Zvap=Zh;
        V=(Zvap*R*T)/P;
        ar_RT=-(log(1-(b/V)))-((a/(R*T*b*(delta_1-delta_2)))*log((1+(delta_1*b/V))/(1+(delta_2*b/V))));
        PHI=exp(ar_RT+Zvap-1-log(Zvap));
        Comp_factor_Z=Zh;
    end
     
else
%if length(Z)=1, then, system is as one phase   
Z=real_roots_EOS;
V=(Z*R*T)/P;
ar_RT=-(log(1-(b/V)))-((a/(R*T*b*(delta_1-delta_2)))*log((1+(delta_1*b/V))/(1+(delta_2*b/V))));
PHI=exp(ar_RT+Z-1-log(Z));
Comp_factor_Z=Z;

end

%___________________________________________________________________________



end %end function