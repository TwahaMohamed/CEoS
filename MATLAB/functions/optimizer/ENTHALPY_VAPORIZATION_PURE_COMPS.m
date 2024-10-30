
% TYPE: type of EOS, SRK or PR
% calculates first enthalpy of the liquind and then enthalpy of vapor

function [H_vap]=ENTHALPY_VAPORIZATION_PURE_COMPS(TYPE,ALFA_function,Zliq,Zvap, T,P,VECTOR_PROPERTIES_PURE_COMP)
% R in kPa.m3/kmol.K
R=8.314;
N=1; % total moles

%props pure comp
Tc=VECTOR_PROPERTIES_PURE_COMP(1);
Pc=VECTOR_PROPERTIES_PURE_COMP(2);
ALPHA_constants(1)=VECTOR_PROPERTIES_PURE_COMP(3);
ALPHA_constants(2)=VECTOR_PROPERTIES_PURE_COMP(4);
ALPHA_constants(3)=VECTOR_PROPERTIES_PURE_COMP(5);
%CP_constants(1)=VECTOR_PROPERTIES_PURE_COMP(6) ;
%CP_constants(2)=VECTOR_PROPERTIES_PURE_COMP(7);
%CP_constants(3)=VECTOR_PROPERTIES_PURE_COMP(8);
%CP_constants(4)=VECTOR_PROPERTIES_PURE_COMP(9);
%CP_constants(5)=VECTOR_PROPERTIES_PURE_COMP(10);


%TYPE: SRK or PR
% Select parameters for the cubic EOS: PR=Peng-Robinson,% SRK=Sove-Redlick_kwong
% calculation of a critical and b
if TYPE == "PR"
    ac=0.45724*(R^2)*((Tc)^2)/(Pc);
    b=0.07780*R*Tc/Pc;
    delta_1=1+sqrt(2);
    delta_2=1-sqrt(2);
elseif TYPE == "SRK"||TYPE == "RK"
    ac=0.42748*(R^2)*((Tc)^2)/(Pc);
    b=0.08664*R*Tc/Pc;
    delta_1=1;
    delta_2=0;    
end
%_______________________________________________________________________________
% calculation of alpha function
if TYPE == "PR"
    
    if ALFA_function == "SOAVE-GN"|| ALFA_function == "SOAVE-AP"||ALFA_function == "SOAVE-CS"||ALFA_function == "Soave_P"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == "TWU-CS"||ALFA_function == "Twu"
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
elseif TYPE == "SRK"||TYPE == "RK"
    
    if ALFA_function == "SOAVE-GN"|| ALFA_function == "SOAVE-AP"|| ALFA_function == "SOAVE-CS"||ALFA_function == "Soave_P"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == "TWU-CS"||ALFA_function == "Twu"
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
end
%calculation of a=ac*alpha
a=alpha.*ac;

%___________________________________________________________________________
% derivatives for enthalpy calculation
Zeta=[Zliq,Zvap];
% Derivatives of alpha function: SOAVE or TWU
ALPHA_constants=[ALPHA_constants(1), ALPHA_constants(2),ALPHA_constants(3)];
[D_alpha_dT,D_alpha_dT_2] = alpha_fx_derivatives(T, Tc,ALPHA_constants,ALFA_function);

for i=1:2
Z=Zeta(i);

V=Z*R*T/P;
   
% Helmholtz residual energy function, F (also g and f)
B=b;
D=a;
g=log(1-(b/V));
f=(1/(R*B*(delta_1-delta_2)))*log((1+(delta_1*b/V))/(1+(delta_2*b/V)));
 
FX = (-N*g)-((D/T)*f);  %(Eq. 60, Page 88 Michelsen Book)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivatives of Cubic EOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%1) Derivatives of B
%1.1) Bi
 Bi=b;

% 1.2) Bij
 Bij=0;

%**************************************************************************
%2) Derivatives of D (depend on the Mixing Rule for Parameter D)
% 2.1) Di
Di=2*a;

%2.2) Calculation of DiT
DiT=2*ac*(D_alpha_dT);

%2.3) Calculation of DTT

DTT=ac*(D_alpha_dT_2);

%2.4) Calculation of DT

DT=0.5*DiT;

%2.5) Calculation of Dij
Dij=0;

%************************************************************************** 
%3) First Order Derivatives of F, g and f
%3.1) fv
fV=-1/(R*(V+delta_1*B)*(V+delta_2*B));
%3.2)fB
fB=-(f+V*fV)/B;
%3.3)gB
gB=-1/(V-B);
%3.4) gV
gV=B/(V*(V-B));
%3.5) FD
FD=-(f/T);
%3.6) FB
FB=-(N*gB)-((D/T)*fB);
%3.7) FV
FV=-(N*gV)-((D/T)*fV);
%3.8) FT
FT=(D/(T^2))*f;
%3.9) Fn
Fn=-g;

%**************************************************************************
% 4) Second Order Partial Derivatives of F, g and f

%4.1) fVV
fVV= (1/(R*B*(delta_1-delta_2)))*((-1/((V+delta_1*B)^2))+(1/((V+delta_2*B)^2)));
% 4.2) fBV
fBV=-(2*fV+V*fVV)/B;
% 4.3) fBB
fBB=-(2*fB+V*fBV)/B;
% 4.4) gBB
gBB=-1/((V-B)^2);
% 4.5) gBV
gBV=1/((V-B)^2);
% 4.6) gVV
gVV=(-1/((V-B)^2))+(1/(V^2));
%4.7) FVV
FVV=-N*gVV-((D/T)*fVV);
% 4.8) FTV
FTV=(D/(T^2))*fV;
% 4.9) FBD
FBD=-fB/T;
% 4.10) FDV
FDV=-fV/T;
% 4.11) FBB
FBB=-(N*gBB)-((D/T)*fBB);
% 4.12) FBV
FBV=-(N*gBV)-((D/T)*fBV);
% 4.13) FDT
FDT=f/(T^2);
% 4.14) FBT
FBT=(D*fB)/(T^2);
% 4.15) FTT
FTT=-2*(FT/T);
% 4.16) FnB
FnB=-gB;
% 4.17) FnV
FnV=-gV;

%**************************************************************************
% FIRST AND SECOND ORDER PARTIAL DERIVATIVES OF THE F FUNCTION
% These derivatives are used for the calculation of thermodynamic properties
%**************************************************************************
%First Order Derivatives
%Equation 66 Chapter 3 Michelsen
dF_dni=Fn;

%Equation 68 Chapter 3 Michelsen
dF_dV=FV;

%Equation 67 Chapter 3 Michelsen
dF_dT=FT+(FD*DT);

% Second Order Derivatives
%Equation 69 Chapter 3 Michelsen
d2F_dnidnj=0;

%Equation 70 Chapter 3 Michelsen
d2F_dnidT=0;

%Equation 71 Chapter 3 Michelsen
d2F_dnidV=FnV;


%Equation 72 Chapter 3 Michelsen
d2F_dT2=FTT+2*FDT*DT+FD*DTT;

%Equation 73 Chapter 3 Michelsen
d2F_dTdV=FTV+FDV*DT;

%Equation 74 Chapter 3 Michelsen
d2F_dV2=FVV;


% end derivatives
%**************************************************************************
%calculation Enathalpy of liquid/vapor

%calculation of Ar
Ar=(R*T)*FX;
% calculation Sr Equation 17 Chapter 2 Michelsen
Sr=R*((-T*dF_dT)-FX);
% calculation residual enthalpy of the phase: liquid or vapor Equation 19 Chapter 2 Michelsen
Hr(i)= Ar+(T*Sr)+(P*V)-(N*R*T);
end

%calculation enthalpy of vaporization= vapor-liquid (kJ/kmol)
H_vap=Hr(2)-Hr(1);

%enthalpy of vaporization (kJ/mol)
H_vap=H_vap/1000;

end %end function