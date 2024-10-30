function [Psat, Zliq, Zvap]=SATURATION_PRESSURE_PURE_COMP(TYPE,Alpha_Function,T,P,VECTOR_PROPERTIES_PURE_COMP)

%this function calculates the vapor pressure of pure components using a
%cubic EOS
Tc = VECTOR_PROPERTIES_PURE_COMP(1);
tolerance=1;
criteria=1e-20;

% calculation EOS parameters
[a,b,~]=CUBIC_EOS_PURE_COMPS_PARAMETERS(TYPE,Alpha_Function,T,VECTOR_PROPERTIES_PURE_COMP);
thershold_value=100000;
loop_counter=0;

while  tolerance > criteria
    P_old=P;    
    %Calculation fugacity coeff's    
    [PHI_LIQUID,Zliq]=FUGACITY_COEFF_PURE_COMP(TYPE,'LIQUID',T,P,a,b,VECTOR_PROPERTIES_PURE_COMP);
    [PHI_VAPOR,Zvap]=FUGACITY_COEFF_PURE_COMP(TYPE,'VAPOR',T,P,a,b,VECTOR_PROPERTIES_PURE_COMP);
    
    if loop_counter==thershold_value

        randomNumber = 0.1899*log(10.484*(T/Tc));
        P=randomNumber;

    else 
    
        tolerance=(PHI_LIQUID-PHI_VAPOR)^2;   
        P=P_old*( PHI_LIQUID/ PHI_VAPOR); 
        loop_counter=loop_counter+1;
    end 

end 
 
Psat=P; 


end