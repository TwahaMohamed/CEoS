function [D_alpha_dT,D_alpha_dT_2] = alpha_fx_derivatives(T, Tc,ALPHA_constants,alpha_fun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation  derivatives of the Twu function 
%_____________________________________________________________________________________________
    % checking to see what alpha function was selected   
    if alpha_fun =="TWU-CS" || alpha_fun == "Twu"
            twu_1= ALPHA_constants(1);
            twu_2=ALPHA_constants(2);
            twu_3=ALPHA_constants(3);
            % Calculation the first part if the Twu fucntion For the first derivative
            part_1a=(T/Tc)^twu_1;
            
            % Calculation the second part if the Twwu fucntion For the first derivative
            part_2a=twu_2*twu_3*((T/Tc)^twu_3)-twu_1;
            
            % Calculation the Third part if the Twu fucntion For the first derivative
            part_3a=exp(twu_2*(1-(T/Tc)^(twu_3)));
            
            % Calculation the first derivative of the Twu function
            D_alpha_dT =-(part_1a*part_2a*part_3a)/T;
            
        %end
        
%Calculation the second derivative of the Twu function 
%_____________________________________________________________________________________________
         
             % Calculation the first part if the Twu fucntion For the second derivative
            part_1b=(T/Tc)^(twu_1);
            
            % Calculation the second part if the Twu fucntion For the second derivative
            part_2b=(twu_2^2)*(twu_3^2)*((T/Tc)^(2*twu_3));
            
            % Calculation the third part if the Twu fucntion For the second derivative
            part_3b=((1-2*twu_1)*twu_2*twu_3-twu_2*twu_3^2)*((T/Tc)^twu_3)+twu_1^2-twu_1;
            
            % Calculation the forth part if the Twu fucntion For the second derivative
            part_4b=exp(twu_2*(1-(T/Tc)^twu_3));
            
           
            % Calculation the second derivative of the Twu function  
            D_alpha_dT_2= (part_1b*(part_2b+part_3b)*part_4b)/(T^2);
            
    end 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Calculation  derivatives of the SOAVE function 
   if alpha_fun =="SOAVE-GN" ||alpha_fun == "SOAVE-AP" || alpha_fun == "Soave"||alpha_fun == "SOAVE-CS"||alpha_fun == "Soave_P"
        m=ALPHA_constants(1);
        alpha_soave=(1+m*(1-sqrt(T/Tc)))^2;
        % Calculation the first derivative of the Soave function   
        D_alpha_dT =-(sqrt(alpha_soave))*m*(1/sqrt(Tc))*(T^-0.5);
        
        % Calculation the first part if the Soave fucntion For the second derivative
        part_1=0.5*(alpha_soave^-0.5)*D_alpha_dT*(T^-0.5);
        
        % Calculation the second part if the Soave fucntion For the second derivative
        part_2=(-0.5)*sqrt(alpha_soave)*(T^-1.5);
        
        % Calculation the second derivative of the Twu function  
        D_alpha_dT_2= -(m/sqrt(Tc))*(part_1+part_2);
   end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    

end