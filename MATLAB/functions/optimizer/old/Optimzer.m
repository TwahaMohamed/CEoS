    function [alphaConstants,fval]=Optimzer(ComponentInfo,ThermodynamicModelInfo,OptimizationOptionsVector)


    initialAlphaValues=ThermodynamicModelInfo{3};
    % % Construct the dynamic file path using genpath to add all subfolders
    % appFolderPath = fileparts(mfilename('fullpath'));
    % addpath(genpath(appFolderPath));            


    % Obtaining critical properties and ideal properties of component
    [Tc,Pc,CP_constants,mw,w,Tb]=CEoS_database_manager(ComponentInfo{1}.CAS(1));


    HeatCapacity_Temperature =rmmissing(ComponentInfo{1}.cp_temperature_K);
    HeatCapacity_ReduceTemperature = rmmissing(ComponentInfo{1}.cp_temperature_K)./Tc;
    HeatCapacity_experimentalValues=rmmissing(ComponentInfo{1}.cp_J_per_K_mol);

    Temperature = { rmmissing(ComponentInfo{1}.saturation_temperature_K),...
                    rmmissing(ComponentInfo{1}.enthalpy_vaporization_temperature_K)...
                    HeatCapacity_Temperature(HeatCapacity_ReduceTemperature<=OptimizationOptionsVector{3})};

    experimentalValues ={ rmmissing(ComponentInfo{1}.saturation_pressure_kPa),...
                          rmmissing(ComponentInfo{1}.enthalpy_vaporization_kJ_per_mol)...
                          HeatCapacity_experimentalValues(HeatCapacity_ReduceTemperature<=OptimizationOptionsVector{3})};

    options = optimoptions('fmincon','UseParallel',true,'Display','iter');
    nonlcon=@(alphaConstants)LeGuennec_constraints(alphaConstants,ThermodynamicModelInfo{2},OptimizationOptionsVector{1});
    fun= @(alphaConstants)ObjectiveFunction(Temperature,experimentalValues,ThermodynamicModelInfo,OptimizationOptionsVector,Tc,Pc,alphaConstants,CP_constants);
    % [alphaConstants,fval]=fmincon(fun,initialAlphaValues,a,b,a_eq,b_eq,lb,ub,nonlcon,options);
    problem=createOptimProblem('fmincon','x0',initialAlphaValues,'objective',fun,'nonlcon',nonlcon,'options',options);

    [alphaConstants,fval ] = run(GlobalSearch,problem);


end

function [OF]= ObjectiveFunction(Temperature,experimentalValues,ThermodynamicModelInfo,OptimizationOptionsVector,Tc,Pc,alphaConstants,CP_constants)
    % Construct the dynamic file path using genpath to add all subfolders
    appFolderPath = fileparts(mfilename('fullpath'));
    addpath(genpath(appFolderPath));         
    %Tags for volume translation
    optimizing_c=false;
    c_opt=0;

    % creating a vector containing pure comp properties; used to solve EOS
    [VECTOR_PROPERTIES_PURE_COMP]= PURE_COMP_PROPERTIES(Tc,Pc,alphaConstants,CP_constants);
    VECTOR_PROPERTIES_PURE_COMP(13:14)=[c_opt,optimizing_c];

    % 
    % [Psat_calc,~,~,~]   = thermo_prop_calc(ThermodynamicModelInfo{1},ThermodynamicModelInfo{2},"Psat",Temperature{1},ComponentInfo{1}.CAS(1),VECTOR_PROPERTIES_PURE_COMP);
    % [Venth_calc,~,~,~]  = thermo_prop_calc(ThermodynamicModelInfo{1},ThermodynamicModelInfo{2},"Venth",Temperature{2},ComponentInfo{1}.CAS(1),VECTOR_PROPERTIES_PURE_COMP);
    % [Cp_calc,~,~,~]     = thermo_prop_calc(ThermodynamicModelInfo{1},ThermodynamicModelInfo{2},"Cp",Temperature{2},ComponentInfo{1}.CAS(1),VECTOR_PROPERTIES_PURE_COMP);

    CEoS=ThermodynamicModelInfo{1};
    alphaFunction =ThermodynamicModelInfo{2};
    evaluationProperties=OptimizationOptionsVector{4};
    Psat_temp  = Temperature{1};
    Venth_temp = Temperature{2};
    Cp_temp    = Temperature{3};

    if evaluationProperties(1)
        parfor i=1:length(Psat_temp)
            % Calculating the saturation pressure for the current temperature value
            [Psat_calc(i,1), ~, ~]=SATURATION_PRESSURE_PURE_COMP(CEoS,alphaFunction,Psat_temp(i),0.0001,VECTOR_PROPERTIES_PURE_COMP);
        end
    else
        Psat_calc =[];
    end

    if evaluationProperties(2)
        parfor j=1:length(Venth_temp)
            % Calculating the saturation pressure for the current temperature value
            [Psat_Venth, Zliq_Venth, Zvap_Venth]=SATURATION_PRESSURE_PURE_COMP(CEoS,alphaFunction,Venth_temp(j),0.0001,VECTOR_PROPERTIES_PURE_COMP);


            % Calculating the enthalpy of vaporization for the current temperature value  
            [Venth_calc(j,1)]=ENTHALPY_VAPORIZATION_PURE_COMPS(CEoS,alphaFunction,Zliq_Venth,Zvap_Venth, Venth_temp(j),Psat_Venth,VECTOR_PROPERTIES_PURE_COMP);

        end
    else
        Venth_calc=[];
    end

    if evaluationProperties(3) 
    parfor k=1:length(Cp_temp)
        % Calculating the saturation pressure for the current temperature value
        [Psat_Cp, Zliq_Cp, ~]=SATURATION_PRESSURE_PURE_COMP(CEoS,alphaFunction,Cp_temp(k),0.0001,VECTOR_PROPERTIES_PURE_COMP);

            % Calculating the heat capacity for the current temperature value  
        [Cp_calc(k,1)]=HEAT_CAPACITY_PURE_COMPS(CEoS,alphaFunction,Zliq_Cp, Cp_temp(k),Psat_Cp,VECTOR_PROPERTIES_PURE_COMP);
    end
    else
        Cp_calc=[];
    end


    if evaluationProperties(1)
        [~,Psat_agg_error]  = error_calc(OptimizationOptionsVector{5}, experimentalValues{1},Psat_calc);
    else
        Psat_agg_error = 0;
    end

    if evaluationProperties(2)
        [~,Venth_agg_error] = error_calc(OptimizationOptionsVector{5}, experimentalValues{2},Venth_calc);
    else
        Venth_agg_error = 0;
    end

    if evaluationProperties(3)
        [~,Cp_agg_error]    = error_calc(OptimizationOptionsVector{5}, experimentalValues{3},Cp_calc);
    else
        Cp_agg_error = 0;
    end

    Psat_Wf =OptimizationOptionsVector{1,2}(1);
    Venth_Wf=OptimizationOptionsVector{1,2}(1);
    Cp_Wf   =OptimizationOptionsVector{1,2}(1);

    OF = (Psat_Wf*Psat_agg_error+Venth_Wf*Venth_agg_error+Cp_Wf*Cp_agg_error)^2;
    agg_error =[Psat_agg_error,Venth_agg_error,Cp_agg_error];
end

%******* NON LINEAR CONSTRAINS (FROM GUENNEC PAPER)*************************
function [c,ceq]=LeGuennec_constraints(alphaConstants,alphaFunction,apply_LeGuennec_constraints)
    if apply_LeGuennec_constraints == true
        if alphaFunction == "Twu"
            twu_1= alphaConstants(1);
            twu_2=alphaConstants(2);
            twu_3=alphaConstants(3);

            X=-3*(twu_3+twu_1-1);
            Y=(twu_3^2)+(3*twu_1*twu_3)-3*twu_3+3*(twu_1^2)-6*twu_1+2;
            Z=-twu_1*((twu_1^2)-3*twu_1+2);

            c(1)= twu_1;
            c(2)=-twu_2*twu_3;
            c(3)=-(1-2*twu_1+2*sqrt(twu_1*(twu_1-1)))+twu_3;
            c(4)= -(4*(Y^3)+4*Z*(X^3)+27*(Z^2)-18*X*Y*Z-((X^2)*(Y^2)));

            ceq=[];            
        elseif alphaFunction == "Soave"
            m= alphaConstants(1);

            c=-m*(1+m);
            ceq=[];
        end
    else
        c=[];
        ceq=[];
    end 

end