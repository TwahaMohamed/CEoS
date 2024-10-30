
% TYPE: type of EOS, SRK or PR
% ALFA_function: type of function for alpha: SOAVE or TWU
% AF_constants: vector containing the parameters for the alpha function
% CP_constant is a vector with 5 elements: A, B, C and D (from Yaws' Tables)
%CP calculated in KJ/Kmol.K



function [CP]=HEAT_CAPACITY_PURE_COMPS(TYPE,ALFA_function,Z, T,P,VECTOR_PROPERTIES_PURE_COMP)
% R in kPa.m3/kmol.K
R=8.314;
N=1; % total moles

%props pure comp
Tc=VECTOR_PROPERTIES_PURE_COMP(1);
Pc=VECTOR_PROPERTIES_PURE_COMP(2);
ALPHA_constants(1)=VECTOR_PROPERTIES_PURE_COMP(3);
ALPHA_constants(2)=VECTOR_PROPERTIES_PURE_COMP(4);
ALPHA_constants(3)=VECTOR_PROPERTIES_PURE_COMP(5);
CP_constants(1)=VECTOR_PROPERTIES_PURE_COMP(6) ;
CP_constants(2)=VECTOR_PROPERTIES_PURE_COMP(7);
CP_constants(3)=VECTOR_PROPERTIES_PURE_COMP(8);
CP_constants(4)=VECTOR_PROPERTIES_PURE_COMP(9);
CP_constants(5)=VECTOR_PROPERTIES_PURE_COMP(10);


%TYPE: SRK or PR
% Select parameters for the cubic EOS: PR=Peng-Robinson,% SRK=Sove-Redlick_kwong
% calculation of a critical and b
if TYPE == string('PR')
    ac=0.45724*(R^2)*((Tc)^2)/(Pc);
    b=0.07780*R*Tc/Pc;
    delta_1=1+sqrt(2);
    delta_2=1-sqrt(2);
elseif TYPE == string('SRK') || TYPE =="RK"
    ac=0.42748*(R^2)*((Tc)^2)/(Pc);
    b=0.08664*R*Tc/Pc;
    delta_1=1;
    delta_2=0;    
end
%_______________________________________________________________________________
% calculation of alpha function
if TYPE == string('PR')
    
    if ALFA_function == "SOAVE-GN"||ALFA_function == "SOAVE-AP"|| ALFA_function == "Soave_P"||ALFA_function == "SOAVE-CS"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == "TWU-CS"  || ALFA_function == "Twu"
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
elseif TYPE == string('SRK') || TYPE =="RK"
    
    if ALFA_function ==  "SOAVE-GN"||ALFA_function == "SOAVE-AP" || ALFA_function == "Soave_P"||ALFA_function == "SOAVE-CS"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == "TWU-CS" || ALFA_function == "Twu"
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
end
%calculation of a=ac*alpha
a=alpha.*ac;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% derivatives of Cubic EOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derivatives of alpha function: SOAVE or TWU
ALPHA_constants=[ALPHA_constants(1), ALPHA_constants(2),ALPHA_constants(3)];
[d_alpha_dT,d2_alpha_dT_2] = alpha_fx_derivatives(T, Tc,ALPHA_constants,ALFA_function);

% calculation  of Volume of the Phase
V=Z*R*T/P;
    
% Helmholtz residual energy function, F (also g and f)
B=b;
D=a;
g=log(1-(b/V));
f=(1/(R*B*(delta_1-delta_2)))*log((1+(delta_1*b/V))/(1+(delta_2*b/V)));
%**************************************************************************
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
DiT=2*ac*(d_alpha_dT);

%2.3) Calculation of DTT

DTT=ac*(d2_alpha_dT_2);

%2.4) Calculation of DT  Equation 104 Chapter 3 Michelsen

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
%4.7) FVV Equation 94 Chapter 3 Michelsen
FVV=-N*gVV-((D/T)*fVV); 
% 4.8) FTV Equation 93 Chapter 3 Michelsen
FTV=(D/(T^2))*fV;
% 4.9) FBD Equation 91 Chapter 3 Michelsen
FBD=-fB/T;
% 4.10) FDV
FDV=-fV/T;
% 4.11) FBB
FBB=-(N*gBB)-((D/T)*fBB);
% 4.12) FBV
FBV=-(N*gBV)-((D/T)*fBV);
% 4.13) FDT
FDT=f/(T^2);
% 4.14) FBT Equation 88 Chapter 3 Michelsen
FBT=(D*fB)/(T^2);
% 4.15) FTT Equation 86 Chapter 3 Michelsen
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
%calculation heat capacity
% (Cv_residual/R) Equation 18 Chapter 2 Michelsen
Cv_res_R=-((T^2)*d2F_dT2)-(2*T*dF_dT);

%(dP/dT) Equation 10 Chapter 2 Michelsen
dP_dT=(-R*T*d2F_dTdV)+(P/T);

%(dP/dV)Equation 9 Chapter 2 Michelsen
dP_dV=(-R*T*d2F_dV2)-(N*R*T/(V^2));

%((Cp_residual-Cv_residual)/R) Equation 19 Chapter 2 Michelsen
Cp_res_Cv_res_R=((-T/R)*((dP_dT)^2)/(dP_dV))-N;

% (Cp_residual)
Cp_res=((Cp_res_Cv_res_R)+(Cv_res_R))*R;

%calculation ideal gas Cp
Cp_ideal_gas=(CP_constants(1)+CP_constants(2)*(T)+CP_constants(3)*(T^2)+CP_constants(4)*(T^3)+CP_constants(5)*(T^4));
% heat capacity
CP=Cp_ideal_gas+Cp_res;

% % Heat capacity analysis checks 
% %_________________________________________________________________________________________________________________________________________
% Cp_res2 =((((-T/R)*(((-R*T*(FTV+FDV*(ac*d_alpha_dT)))+(P/T))^2)/((-R*T*d2F_dV2)-(N*R*T/(V^2))))-N)+(-((T^2)*(FTT+2*FDT*(ac*d_alpha_dT)+FD*(ac*d2_alpha_dT_2)))-(2*T*(FT+(FD*(ac*d_alpha_dT))))))*R;
% 
% AA = d_alpha_dT*(d_alpha_dT*(R*T*FDV*ac)^2-2*(P/T)*(R*T*FDV*ac)+2*(R*T*FTV)*(R*T*FDV*ac))+((P/T)^2-2*(P/T)*(R*T*FTV)+(R*T*FTV)^2);
% 
% %AA = (P/T)^2-2*(P/T)*(R*T*FTV)-2*(P/T)*d_alpha_dT*(R*T*FDV*ac)+(R*T*FTV)^2+2*(R*T*FTV)*d_alpha_dT*(R*T*FDV*ac)+(d_alpha_dT*(R*T*FDV*ac))^2
% 
% %Y1 =(P/T);
% %Y2 =(R*T*FTV);
% %Y3 =(R*T*FDV*ac*d_alpha_dT);
% %AA = Y1^2 - 2*Y1*Y2 - 2*Y1*Y3 + Y2^2 +2*Y2*Y3 +Y3^2;
% 
% %AA= (Y1-Y2-Y3)^2;
% 
% %AA =((-R*T*(FTV+FDV*(ac*d_alpha_dT)))+(P/T))^2;
% 
% AA_actual = (dP_dT)^2;
% AA_Test = AA_actual-AA;
% 
% BB =(d2F_dV2*(-R*T)-((N*R*T)/(V^2)));
% BB_actual = ((-R*T*d2F_dV2)-(N*R*T/(V^2)));
% BB_Test = BB_actual-BB;
% 
% CC_1 = ((d_alpha_dT*(d_alpha_dT*(-T/R)*(R*T*FDV*ac)^2-2*(P/T)*(-T/R)*(R*T*FDV*ac)+2*(-T/R)*(R*T*FTV)*(R*T*FDV*ac))+(-T/R)*((P/T)^2-2*(P/T)*(R*T*FTV)+(R*T*FTV)^2))/dP_dV)-N;
% CC_actual = (-T/R)*(((-R*T*(FTV+FDV*(ac*d_alpha_dT)))+(P/T))^2)/dP_dV-N;
% CC_Test = CC_actual-CC_1;
% 
% A1 = ((-T/R)*(R*T*FDV*ac)^2)/dP_dV;
% A2 = (-2*(P/T)*(-T/R)*(R*T*FDV*ac)+2*(-T/R)*(R*T*FTV)*(R*T*FDV*ac))/dP_dV;
% A3 = ((-T/R)*((P/T)^2-2*(P/T)*(R*T*FTV)+(R*T*FTV)^2))/dP_dV;
% 
% CC_2 = (d_alpha_dT*(d_alpha_dT*A1+A2)+A3)-N;
% CC_Test_2=CC_actual-CC_2;
% 
% DD=-(T^2)*FTT-d_alpha_dT*(2*FDT*ac)*(T^2)-d2_alpha_dT_2*(FD*ac)*(T^2);
% DD_actual =-((T^2)*(FTT+2*FDT*(ac*d_alpha_dT)+FD*(ac*d2_alpha_dT_2)));
% DD_Test =DD_actual- DD;
% 
% 
% EE = - (2*T*FT)-d_alpha_dT*(2*T*FD*ac);
% EE_actual = -(2*T*(FT+(FD*(ac*d_alpha_dT))));
% EE_Test = EE_actual-EE;
% 
% 
% %Cp_res3 =(d_alpha_dT^2*X1+d_alpha_dT*X2+X3- N - (T^2)*FTT-d_alpha_dT*(2*FDT*ac)*(T^2)-d2_alpha_dT_2*(FD*ac)*(T^2) - (2*T*FT)-d_alpha_dT*(2*T*FD*ac))*R;
% %Cp_res3 = ( -d2_alpha_dT_2*(FD*ac)*(T^2) + (d_alpha_dT^2*X1 + d_alpha_dT*X2 - d_alpha_dT*(2*FDT*ac)*(T^2) -d_alpha_dT*(2*T*FD*ac))+(X3- N - (T^2)*FTT - (2*T*FT)) )*R;
% Cp_res3 =-d2_alpha_dT_2*(R*FD*ac)*(T^2) + d_alpha_dT*(d_alpha_dT*(A1*R)+(A2*R-(2*FDT*ac*R)*(T^2)-(2*T*FD*ac*R)))+R*(A3- N - (T^2)*FTT - (2*T*FT));
% Cp_res3_test = Cp_res-Cp_res3;
% 
% 
% X1 = (R*FD*ac)*(T^2);
% X2 = (((-T/R)*(R*T*FDV*ac)^2)/dP_dV)*R;
% X2_test= A1*R-X2;
% X3 = ((-2*(P/T)*(-T/R)*(R*T*FDV*ac)+2*(-T/R)*(R*T*FTV)*(R*T*FDV*ac))/dP_dV)*R-(2*FDT*ac*R)*(T^2)-(2*T*FD*ac*R);
% X3_test= (A2*R-(2*FDT*ac*R)*(T^2)-(2*T*FD*ac*R))-X3;
% X4 = R*(A3- N - (T^2)*FTT - (2*T*FT));
% Cp_res4= -d2_alpha_dT_2*X1 + d_alpha_dT*(d_alpha_dT*X2+X3)+X4;
% 
% Cp_res4_test = Cp_res- Cp_res4;


end %end function