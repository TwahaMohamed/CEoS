function [screenedData,detectedOutlines,limits] = removeOutliers(Tbl,NumTrees,tau,binSize,k,plotAll,names,var_names)
%% REMOVEOUTLIERS Detect Outliers Using Quantile Regression
%% Grow Quantile Random Forest
% Grow a bag  regression trees using TreeBagger.
    % Function detects outliers in data set Using Quantile Regression. 
    % The detected outliers are then removed from the data set based 
    % on the defined limits from the Quantile Regression.
    
    % Tbl: A 2d table that contains all of the data for the x, and y axis
    % NumTrees: the number of decision trees created to 100 is a good default
    % value (used to avoid overfiting)
    % tau: provied the estemanted values needed to creat a quantile
    % regression ( 0.25 0.5 0.75) are a good default values
    
    % binSize: determins the number of section the inlet data will be split
    % into 50 is a good default value
    % plotAll: must be "true" or "false" , defined "plotAll” as "true" if you
    % want the function to show all plots including optional plots, or define
    % as false if you only the function to  only show required plots.
    % i- allow for use in loops
    Tbl.Properties.VariableNames = {'t','y'};
    Mdl = TreeBagger(NumTrees,Tbl,'y','Method','regression');
%% 
% Mdl is a TreeBagger ensemble.
%% Predict Conditional Quartiles and Interquartile Ranges
% Using quantile regression, estimate the conditional quartiles of "tau" equally 
% spaced values within the range of t.
    
    predT = linspace(min(Tbl.t),max(Tbl.t),binSize)';
    quartiles = quantilePredict(Mdl,predT,'Quantile',tau);
%% 
% quartiles is a n-by-3 matrix of conditional quartiles. Rows correspond to 
% the observations in t, and columns correspond to the probabilities in tau.
% 
% On the scatter plot of the data, plot the conditional mean and median responses.
    %predicting the averages 
    meanY = predict(Mdl,predT);
    %Optinal Ploting ( Mean,Quartiles, and data plot)
    if plotAll == true
        figure(randi([1 1000]));
    
        plot(Tbl.t,Tbl.y,'LineWidth',2);
        hold on
        plot(predT,[quartiles(:,2) ,meanY],'LineWidth',2);
        title(names+" Median and Mean response" )
        legend('Data','Median response','Mean response',...
            'Location','NorthWest');
        hold off
    end
%% 
% 
% 
% Compute the conditional , , and .
    iqr = quartiles(:,3) - quartiles(:,1);
%     k = 1.5;
    f1 = quartiles(:,1) - k*iqr;
    f2 = quartiles(:,3) + k*iqr;
%% 
% k = 1.5 means that all observations less than f1 or greater than f2 are considered 
% outliers, but this threshold does not disambiguate from extreme outliers. A 
% k of 3 identifies extreme outliers.
%% Compare Observations to Fences (Optional  Plot)
% Plot the observations and the fences.
    if plotAll == true
        figure(randi([1 1000]));
        plot(Tbl.t,Tbl.y,'.');
        
        hold on
        
        plot(predT,[f1 f2]);
        legend('Data','F_1','F_2','Location','NorthWest');
        axis tight
        title(names+' Outlier Detection Using Quantile Regression');
        hold off
    end 
%% 
% All simulated outliers fall outside , and some observations are outside this 
% interval as well.
% 
% Copyright 2012 The MathWorks, Inc.
%% Removal of Outliers
% Filling in points for the limits of  so as to match points alon the x-axis
    %Linearly interpolates based on the existing points in "f1" and "perdT" 
    % given the desired estimation points from "Tbl.t"
    f1_interp = interp1(predT,f1,Tbl.t);
    
    %Linearly interpolates based on the existing points in "f1" and "perdT" 
    % given the desired estimation points from "Tbl.t"
    f2_interp = interp1(predT,f2,Tbl.t);
%% Optional ploting ( data and Linearly interpolates limits)
    if plotAll == true
        %Creating a figure to plot on
        figure(randi([1 1000]));
        %Creating a plot of the data in it's current sate
        plot(Tbl.t,Tbl.y,'.');
        %allows for more plots to be added to current figure
        hold on
        % adding a plot of the lower limit 
        plot(Tbl.t(:,1),f1_interp,".",predT,f1,"o")
        % adding a plot of the upper limit 
        plot(Tbl.t(:,1),f2_interp,".",predT,f2,"o")
        
        pause(0.1);
         % adding a legend to the plot
        legend("Data Point","Lower Limit Point","f1","Upper Limit Point","f2","location","best")
        title(names+" limits plot")
        %seting the x label name
        xlabel(var_names(1))
    
        %setting the y label name
        ylabel(var_names(2))
        %prevents more plots from being added to current figure
        hold off
    end
%% 
% 
    % setting the loop variable "tempLoop" to equal the number of rows in
    % "Tbl" 
    tempLoop =size(Tbl,1);
    
    % reallocating the size of "rowDel" to improve speed
    rowDel = zeros(tempLoop);
    %Comparing the interpolated values of “f1_interp” and “f2_interp” to 
    % see if the data set has values that are outside of these limits.
    for i =1:tempLoop
            %Comparing the interpolated values of "f1_interp" at the 
            % current temperature to see if it Lises outside of the lower limit 
            if  Tbl.y(i) < f1_interp(i)
                % I the case that the value does lie outside the lower
                % limit we store this value in "rowDel"
                rowDel(i) = i;
                
            %Comparing the interpolated values of "f2_interp" at the 
            % current temperature to see if it Lises outside of the upper limi             
            elseif Tbl.y(i) > f2_interp(i)
    
                % I the case that the value does lie outside the upper
                % limit we store this value in "rowDel"
                rowDel(i) = i;
                
            end 
    end
    % Removing all the zero which represent rows containing data that are
    % within the limits of "f1" and "f2"
    rowDel(rowDel==0)=[];
    
    % Filtering "Tbl" for values that have been marked as outliers and
    % storing this values in "Outliers"
    Outliers = Tbl(rowDel,:);
    % Filtering "Tbl" for values that have been marked as outliers and
    % deleting these values from "Tbl"    
    Tbl(rowDel,:) =[];
    
%% Optional ploting (data and oultiers)
    if plotAll == true
        %Creating a figure to plot on
        figure(randi([1 1000]));
        %Creating a plot of the data with no outliers
        plot(Tbl.t,Tbl.y,'.');
       
        %allows for more plots to be added to current figure
        hold on
        % adding a plot of the outliers 
        plot(Outliers.t,Outliers.y,'r*');
        title(names+" data with oulier plot")
        %seting the x label name
        xlabel(var_names(1))
    
        %setting the y label name
        ylabel(var_names(2))
        %prevents more plots from being added to current figure
        hold off
    end
%% Required Plotting (data, limits , and outliers) not used in GUI
    % %Creating a figure to plot on
    % figure(randi([1 1000]));
    % 
    % %Creating a plot of the data with no outliers
    % plot(Tbl.t,Tbl.y,'.');
    % 
    % %allows for more plots to be added to current figure
    % hold on 
    % 
    % % adding a plot of the limits
    % plot(predT,[f1 f2]);
    % 
    % % adding a title
    % title(names+' Outlier Detection Using Quantile Regression');
    % 
    % %seting the x label name
    % xlabel(var_names(1))
    % 
    % %setting the y label name
    % ylabel(var_names(2))
    % 
    % % adding a plot of the outliers 
    % plot(Outliers.t,Outliers.y,'r*');
    % 
    % % adding a legend to the plot
    % legend('Data','F_1','F_2',"Outliers",'Location','best');
    % 
    % %prevents more plots from being added to current figure
    % hold off
    % 
    % % setting axis as "tight" maximizing used plot space
    % axis tight
    % 
    % figure(round(rand()*1000));
    % %Creating a plot of the data with no outliers
    % plot(Tbl.t,Tbl.y,'.');
    % 
    % % setting axis as "tight" maximizing used plot space
    % axis tight
    % 
    % % adding a title
    % title(names+' Outlier Removed');
    % 
    % %seting the x label name
    % xlabel(var_names(1))
    % 
    % %setting the y label name
    % ylabel(var_names(2))
    % 
    % %setting function output to be "Tbl"
    % 
    % Tbl_Output=Tbl;
    
%% Export Setup
    screenedData = Tbl;
    detectedOutlines =Outliers;
    predT_table = table(predT,VariableNames={'Temperature_K'});
    f1_table = table(f1,VariableNames={'upper_Fence'});
    f2_table = table(f2,VariableNames={'lower_Fence'});
    limits= [predT_table,f1_table,f2_table];
end