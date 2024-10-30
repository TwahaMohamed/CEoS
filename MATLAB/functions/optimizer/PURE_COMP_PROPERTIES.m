% this function organizes the pure component parameters into a single vector


function [VECTOR_PROPERTIES_PURE_COMP]= PURE_COMP_PROPERTIES(Tc,Pc,ALPHA_constants,CP_constants)

%this is a row vector

% element 1: Tc
VECTOR_PROPERTIES_PURE_COMP(1)= Tc;

% element 2: Pc
VECTOR_PROPERTIES_PURE_COMP(2)= Pc;

%elements 3,4 & 5 constants for alpha function
VECTOR_PROPERTIES_PURE_COMP(3)= ALPHA_constants(1);
VECTOR_PROPERTIES_PURE_COMP(4)= ALPHA_constants(2);
VECTOR_PROPERTIES_PURE_COMP(5)= ALPHA_constants(3);

%elements 6,7,8,9 & 10 constants for ideal gas heat capacity
VECTOR_PROPERTIES_PURE_COMP(6)= CP_constants(1);
VECTOR_PROPERTIES_PURE_COMP(7)= CP_constants(2);
VECTOR_PROPERTIES_PURE_COMP(8)= CP_constants(3);
VECTOR_PROPERTIES_PURE_COMP(9)= CP_constants(4);
VECTOR_PROPERTIES_PURE_COMP(10)= CP_constants(5);


end