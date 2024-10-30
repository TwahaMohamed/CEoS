% This function calculate the cubic EOS parameters for pure components
%ALFA_function: type of alpha function: SOAVE or TWU
%ALPHA_constants: constants for alpha function. for SOAVE only 1 parameter vector. for TWU, 3 parameter vector
% type is the EOS: SRK or PR

function [a,b,alpha]=CUBIC_EOS_PURE_COMPS_PARAMETERS(TYPE,ALFA_function,T,VECTOR_PROPERTIES_PURE_COMP)

%Tc,Pc,ALPHA_constants
Tc=VECTOR_PROPERTIES_PURE_COMP(1);
Pc=VECTOR_PROPERTIES_PURE_COMP(2);
ALPHA_constants(1)=VECTOR_PROPERTIES_PURE_COMP(3);
ALPHA_constants(2)=VECTOR_PROPERTIES_PURE_COMP(4);
ALPHA_constants(3)=VECTOR_PROPERTIES_PURE_COMP(5);


% R in kPa.m3/kmol.K
R=8.314;
N=1; % total moles

%TYPE: SRK or PR
% Select parameters for the cubic EOS: PR=Peng-Robinson,% SRK=Sove-Redlich-Kwong
% calculation of a critical and b
if TYPE == "PR"
    ac=0.45724*(R^2)*((Tc)^2)/(Pc);
    b=0.07780*R*Tc/Pc;
    delta_1=1+sqrt(2);
    delta_2=1-sqrt(2);
elseif TYPE == "SRK"||TYPE =="RK"
    ac=0.42748*(R^2)*((Tc)^2)/(Pc);
    b=0.08664*R*Tc/Pc;
    delta_1=1;
    delta_2=0;    
end
%_______________________________________________________________________________
% calculation of alpha function
if TYPE == "PR"
    
    if ALFA_function == "SOAVE-GN"||ALFA_function == "SOAVE-AP"||ALFA_function =="Soave_P"||ALFA_function =="SOAVE-CS"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc)))^2;
    elseif ALFA_function == "TWU-CS"||ALFA_function == "Twu"
    twu_1= ALPHA_constants(1);
    twu_2=ALPHA_constants(2);
    twu_3=ALPHA_constants(3);
    Tr=T/Tc;
    alpha=(Tr^(twu_1))*exp(twu_2*(1-(Tr^twu_3)));    
    end
    
elseif TYPE == "SRK"||TYPE =="RK"
    
    if ALFA_function == "SOAVE-GN"||ALFA_function == "SOAVE-AP"||ALFA_function =="Soave_P"||ALFA_function =="SOAVE-CS"
    m=ALPHA_constants(1);
    alpha=(1+m*(1-sqrt(T/Tc))).^2;
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

end