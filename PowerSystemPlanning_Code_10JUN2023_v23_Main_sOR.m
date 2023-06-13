clear; clc;
format shortG
import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
diary off;                                  % In case there is an error an it returns before reaching the end in the previous iteration
% Keep Track of time
initial = datetime;

% Directory
Database_name  = "Benchmark6bus";           %% Change name depending on which information we want to work with.
Route = "C:\Users\newsu\Documents\MATLAB" +"\\"+ Database_name;           %% The master directory will be the place where this code is stored.
cd(Route);                                  %% Change the directory to aim to the corresponding database.

% load demand into a variable
demand_base = readmatrix("Demand.xlsx");
opt_year    = "2050";                       %% opt_year must be a string.

%% Global information
% Be careful with globals: https://matlab.fandom.com/wiki/FAQ#Are_global_variables_bad.3F Pay especial attention to point 3.
global crit_time blocks rate Curt_penalty qt_days fst_study_horizon last_study_horizon stages scenarios
crit_time       = 0.1;       %% duration of first demand block
blocks          = 5;     
rate            = 0.1;          %% interest rate for projects. [0-1].
Curt_penalty    = 500;          %% we are declaring that curtailment has a value of 500 USD/MWh
qt_days         = 7;            % We define how long will be the vector for days to be tested depending all of the criteria
fst_study_horizon = 2021;
last_study_horizon = str2double(opt_year);  %% This can definitely be more elegant.
stages          = 12;                        %% stages should be seen as how do we divide each year. Now we are only dividing the year in months, or not dividing it at all.
scenarios       = 1;            % Number of renewable scenarios

% We make sure that the directory in which we will gonna graph is empty. To do so, in case information exists, we delete completely the directory. 
% In case it doesn't exist. We create the directory.
Reports_Directory = Route + "\\Reports";
if not(isfolder(Reports_Directory))
    mkdir(Reports_Directory);
elseif isfolder(Reports_Directory)
    rmdir(Reports_Directory,'s');
    mkdir(Reports_Directory);
end

% Start recording the output to a log file
logFilename = Route + '\Reports\LogOutput.txt';
diary(logFilename);

%% cluster demand
base_LDC = clustering_block(demand_base);       % Maybe we can think as stages as month or year. So we could have the option of making this yearly or annually.
%% Create LDC Curve for base year 
Curve_LDC_base = [];        % Curve_LDC_base is a vector with the cluster powers
for t=1:stages*blocks
    Curve_LDC_base = [Curve_LDC_base;ones(round(base_LDC.Hours(t)),1)*base_LDC.Cluster(t)];
end

% Plot the LDC vs chronological demand for year 2019
graph_LDC("2019",Curve_LDC_base,demand_base); % We use the function because we intend to use it several times

%% Apply demand growth and obtain a set of new set of demands for the following years
% using information from the document: Generation Indicative Expansion Plan
% 2022-2031 Honduran Independent System Operator

% Invoke function to have tables with LDC and Smooth Demand
[Global_LDC_Demand,Global_Chronological_demand] = global_demands(Curve_LDC_base,demand_base);

% Graph for year to be optimized
graph_LDC(opt_year,Global_LDC_Demand.(opt_year),Global_Chronological_demand.(opt_year));

%% Get the hour-block map and the equivalent generation for all the blocks
% Hour-block map
hbm = hourblockmap(demand_base,base_LDC,opt_year);

% Equivalent generation for renewables for each block
[renew_gen_block, renewables_gen] = RENEW_GEN_BLK(hbm);

%% Create the Long Term Load and store it in a Table for optimization
% Obtain the LDC distribution for the last year of study
opt_LDC = clustering_block(Global_Chronological_demand.(opt_year));
% Create the GEP_Load
GEP_Load = calc_GEP_Load(Global_LDC_Demand,opt_LDC);

%% Renewable Commitment quota
Renew_commitment = calc_renew_commitment(opt_year);  %% SUGGESTION: Perhaps we just need milestone years according to plan

%% Iterations counter
iterations       = 0;
% Load general information from Database
DBB_info = call_BDD_info(Database_name + ".xlsx");
%% Flexibility and Capacity indicators
% Firm Capacity Reserve
sOR              = 0.10*ones(1,last_study_horizon-fst_study_horizon+1);    % Reserve Initial used for the initial problem
% Ramping
ramp_initial    = ramp_calc([],iterations,DBB_info);
ramp_capacity   = ramp_initial*ones(1,last_study_horizon-fst_study_horizon+1);  % ramp capacity for every year within the study horizon
% Define Renew_Cap
Renew_Cap = renew_cap_calc([],[],iterations,DBB_info,GEP_Load,[],[]);

% Information for peak_Yearly_Installed capacity graph
flag_peakLoad_InstalledCapacity = 0;    % Create a flag that will be used for the Peak Demand per year vs Peak Installed Capacity per Iteration
Global_Feasibility = 0;                 % Create a flag that will keep track of feasibility in the problem
flag_UCP_Stage = "First";               % Create a flag that will be used to check if we are checking Firm Power and Constrain Renewables
yearly_Installed_Capacity = [];         % Initialize the vector at zero
legend_description = [];                % Initialize the vector as empty

%% Create optimizer objects
yalmip('clear');
construct_initial_time = datetime;
GEP_optimizer = GEP_optimizer_creator(GEP_Load,Renew_commitment,renew_gen_block,DBB_info,opt_year);    % Higher hierarchy GEP optimizer object
construct_final_time = datetime;
GEP_construct_time = construct_final_time - construct_initial_time;
disp(GEP_construct_time)
UCP_optimizer = UCP_optimizer_creator(DBB_info);         % Verifier UCP optimizer object

%% Main loop
while true                  % Infinite loop that will only be broken if we get to the finish line.
    iterations = iterations + 1;
    disp("This is iteration: " + num2str(iterations));
    % Display Time    
    final = datetime;
    disp(datetime);

    %% First Stage: Generation Expansion Problem
    disp("First Stage: Generation Expansion Problem")
    close all hidden;        % Attempt to free up memory

    % Evaluate the GEP Optimizer:
    [Investment,Output_G_exist,Output_G_cand,Cost,GEP_optimizer] = GEP_evaluate(GEP_optimizer,sOR,ramp_capacity,DBB_info,Renew_Cap,iterations);
    
    % Graph Findings
    % Graph of estimates of Operation
    graph_GEP(Output_G_exist,Output_G_cand,iterations,"off",GEP_Load,DBB_info,sOR,Investment);   %% Graph a Report with the information of the recent findings.
    % Graph of Firm power evolution over time
    [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,DBB_info,sOR);


    %% Second stage: Unit-Commitment Stage    
    disp(flag_UCP_Stage + " Stage UCP verification")
    for t = 1:last_study_horizon-fst_study_horizon+1
        for s = 1:scenarios
            % Signal start of for loop
            disp("Year currently under study: "+num2str(fst_study_horizon+t-1) + "- Scenario: " + s);            %% Flag to keep track of code advancements.
            
            % Create demand:
            % 1. Identify days for demand:
            test_days = ident_days_to_test(Global_Chronological_demand,num2str(fst_study_horizon+t-1),renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),Investment(:,t),length(demand_base),DBB_info);  
    
            % 2. Use identified days to obtain the demand:
            [Week_Chrono_LNL,Week_Chrono_MNL,Week_Chrono_LRG,Week_Chrono_MRG,Week_Chrono_MDM,Week_Chrono_RAN] = Load_Chrono_to_test(num2str(fst_study_horizon+t-1),Global_Chronological_demand,test_days);
    
            %% Unit Commitment Problem
            % General Feasibility variables to control the loop between GEP and UC
%             Feasibility_FirmPower = 0;
%             Feasibility_Ramping = 0;
            Feasibility = 0;
            Aggregate_Curtailment = 0;
            % 1. Test the UC Problem for Most Demand
            close all hidden;        % Attempt to free up memory        
            % Evaluate the UCP optimizer
            [Output_G_cand,Output_G_exist,MDM_marker,UCP_optimizer,MDM_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.most_Demand,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MDM,UCP_optimizer,"Most Demand",DBB_info);
            % Decide feasibility in case feasible, graph.
%             Feasibility_FirmPower = Feasibility_FirmPower + MDM_marker;
            Feasibility = Feasibility + MDM_marker;                         Global_Feasibility = Global_Feasibility + MDM_marker;
            Aggregate_Curtailment = Aggregate_Curtailment + MDM_Curtailment;
            if MDM_marker == 0
                graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Demand","off",DBB_info,Investment,s);
            else
                disp("Most Demand is not feasible.")
            end
            
            % 2. Test the UC Problem for Most Net Load
            close all hidden;        % Attempt to free up memory       
            % Evaluate the UCP optimizer
            [Output_G_cand,Output_G_exist,MNL_marker,UCP_optimizer,MNL_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.most_NetLoad,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MNL,UCP_optimizer,"Most Net Load",DBB_info);
            % Decide feasibility in case feasible, graph.
%             Feasibility_FirmPower = Feasibility_FirmPower + MNL_marker;
            Feasibility = Feasibility + MNL_marker;                         Global_Feasibility = Global_Feasibility + MNL_marker;
            Aggregate_Curtailment = Aggregate_Curtailment + MNL_Curtailment;
            if MNL_marker == 0
                graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Net Load","off",DBB_info,Investment,s);
            else
                disp("Most Net Load is not feasible.")
            end
    
            % 3. Test the UC Problem for Least Renewable Generation
            close all hidden;        % Attempt to free up memory
            % Evaluate the UCP optimizer
            [Output_G_cand,Output_G_exist,LRG_marker,UCP_optimizer,LRG_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.least_Renew_gen,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_LRG,UCP_optimizer,"Least Renewable Generation",DBB_info);
            % Decide feasibility in case feasible, graph.
%             Feasibility_FirmPower = Feasibility_FirmPower + LRG_marker;
            Feasibility = Feasibility + LRG_marker;                         Global_Feasibility = Global_Feasibility + LRG_marker;
            Aggregate_Curtailment = Aggregate_Curtailment + LRG_Curtailment;
            if LRG_marker == 0
                graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Least Renewable Generation","off",DBB_info,Investment,s);
            else
                disp("Least Renewable Generation is not feasible.")
            end
            
            % 4. Test the UC Problem for Least Net Load
            close all hidden;        % Attempt to free up memory        
            % Evaluate the UCP optimizer
            [Output_G_cand,Output_G_exist,LNL_marker,UCP_optimizer,LNL_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.least_NetLoad,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_LNL,UCP_optimizer,"Least Net Load",DBB_info);
            % Decide feasibility in case feasible, graph.
%             Feasibility_Ramping = Feasibility_Ramping + LNL_marker;
            Feasibility = Feasibility + LNL_marker;                         Global_Feasibility = Global_Feasibility + LNL_marker;
            Aggregate_Curtailment = Aggregate_Curtailment + LNL_Curtailment;
            if LNL_marker == 0
                graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Least Net Load","off",DBB_info,Investment,s);
            else
                disp("Least Net Load is not feasible.")
            end
    
            % 5. Test the UC Problem for Most Renewable Generation
            close all hidden;        % Attempt to free up memory        
            % Evaluate the UCP optimizer
            [Output_G_cand,Output_G_exist,MRG_marker,UCP_optimizer,MRG_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.most_Renew_gen,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_MRG,UCP_optimizer,"Most Renewable Generation",DBB_info);
            % Decide feasibility in case feasible, graph.
%             Feasibility_Ramping = Feasibility_Ramping + MRG_marker;
            Feasibility = Feasibility + MRG_marker;                         Global_Feasibility = Global_Feasibility + MRG_marker;
            Aggregate_Curtailment = Aggregate_Curtailment + MRG_Curtailment;
            if MRG_marker == 0
                graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Most Renewable Generation","off",DBB_info,Investment,s);
            else
                disp("Most Renewable Generation is not feasible.")
            end
    
%             % 6. Test the UC Problem for a Random Week
%             close all hidden;        % Attempt to free up memory        
%             % Evaluate the UCP optimizer
%             [Output_G_cand,Output_G_exist,RAN_marker,UCP_optimizer,RAN_Curtailment] = UCP_evaluate(renewables_gen.("S"+s).("Y"+num2str(fst_study_horizon+t-1)),test_days.Random,num2str(fst_study_horizon+t-1),Investment(:,t),Week_Chrono_RAN,UCP_optimizer,"Random days",DBB_info);
%             % Decide feasibility in case feasible, graph.
%             Aggregate_Curtailment = Aggregate_Curtailment + RAN_Curtailment;
%             if RAN_marker == 0
%                 graph_UCP(Output_G_exist,Output_G_cand,sOR(t),num2str(fst_study_horizon+t-1),"Random","off",DBB_info,Investment,s);
%             else
%                 disp("Random days is not feasible.")
%             end
    
            % For user readings:
            % disp("The average curtailment for the study weeks is: " + string(v.format(Aggregate_Curtailment/sum([MDM_marker;MNL_marker;LRG_marker;LNL_marker;MRG_marker;RAN_marker] == 0))) + newline)
            disp("The average curtailment for the study weeks is: " + string(v.format(Aggregate_Curtailment/sum([MDM_marker;MNL_marker;LRG_marker;LNL_marker;MRG_marker] == 0))) + newline)
        end

        % Output Information for User
        disp("Year "+ num2str(fst_study_horizon+t-1) + " finished." +newline);            %% Flag to keep track of code advancements.
    end

    %% Feedback
    if Feasibility ~= 0         % If one or more of the previous UCP are not feasible, we will increase the amount of reserve, to try to make them feasible within the next iterations. 
        disp("System Operating Reserve is incremented from "+ num2str(sOR(t)*100) +"% to "+ num2str((sOR(t)+.1)*100)  +"%."+newline)
        sOR = sOR + 0.1;        % SUGGESTION: Device a way to find the optimal amount of reserve. Maybe after it becomes feasible, apply a bisection method to find the best sOR within a given 1%. Perhaps for the next while, apply the same logic, but in increments of 1%.
        Feasibility = 0;
    else
        flag_peakLoad_InstalledCapacity = flag_peakLoad_InstalledCapacity+1;        % Change the flag for we to graph and exit
        [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,DBB_info,sOR);
        break           % We attempt to break the while loop
    end
end

% Stop recording the output
diary off;

final = datetime;
% Save the data to a mat file for review even after MatLab interface has been closed.
% clear DBB_info;         % This is an nonessential variable, used only to handle in an object all the database information. However, this are know parameters stored initially in Excel files, so it is useless to store them in our final mat file. Other variables are nonessential, and could as well be deleted.
save(Reports_Directory +"\\"+ Database_name + ".mat");

% Interact with user to let him know it is over
for d = 1:3
    pause(2)
    beep
end

disp("Algorithm is done. Running clock time is: " + string(final-initial));

%% Functions 

function master_LDC = clustering_block(demand_base)     %% FIX: stages should be able to be 1 for annually or 12 for monthly
    % We will base most of the reasoning in here in the file: Pruebas Herramienta para determinar BH_ 12 MAR 2019
    global crit_time blocks stages
    
    % Divide the entries into its given months
    % 1. Create a vector of dates that will be associated to the demand entries. Given that it is not important the year, we will do it with 1900 because it is not a leap year
    % The first entry is the 1900/JAN/1 0:00:00
    time_stamp(1,1) = datetime(1900,1,1,0,0,0);
    % The rest are just add an hour as defined 1/24
    for i = 2:length(demand_base)
        time_stamp(i,1) = time_stamp(i-1) + 1/24;
    end

    % 2. Make a vector with each month
    month_stamp = month(time_stamp);
    
    % 3. Attempt to apply the clustering technique to each individual month
    master_LDC = [];
    for m = 1:stages
        % a. Define the average duration time for the rest of the clusters other than the critical time
        rest_duraci = (1-crit_time)/(blocks-1);
        % b. Extract the demand associated with the current we are currently
        if stages == 1
            current_demand = demand_base;
        else
            current_demand = demand_base(find(month_stamp == m));
        end
        % c. Sort from greatest to smallest
        sorted_demand = sort(current_demand, 'descend');
        % d. Define Centroids and their corresponding lengths
        Lengths     = zeros(blocks,1);
        Centroids   = zeros(blocks,1);
        % e. Initialization:
        % e.1 Peak Block definition
        Lengths(1)  = round(length(current_demand)*crit_time);    % Length or duration of block 1
        Centroids(1)= mean(sorted_demand(1:Lengths(1)));                % We are stating that centroid 1 is the mean from the first data till the length 1.
        % e.1 Hot start of the duration of the rest of the centroids
        for i = 2:blocks
            if i == blocks  % In case we are at the last block
                Lengths(i)  = Lengths(i-1) + round(length(current_demand)*rest_duraci)-1; % We are accumulating the reference of the duration length
                Centroids(i)= mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));                   % we average those said duration lengths
            else            % In case for the rest of the blocks
                Lengths(i)   = Lengths(i-1) + round(length(current_demand)*rest_duraci);
                Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
            end
        end
        % f. Define a table of demand and distances to all Centroids
        Table = [];
        for i = 1:blocks
            Table(:,i) = abs(current_demand - Centroids(i));  % Distance to each centroid
        end
        [~,Duration] = min(Table,[],2);         % Assign flag to assign Centroid
        [duraci,~] = groupcounts(Duration);     % Count the duration of each Block
        % g. Update Lengths and Centroids
        Lengths(1) = duraci(1);
        Centroids(1) = mean(sorted_demand(1:Lengths(1)));
        % h. for loop to update blocks 2 till the end
        for i = 2:blocks
            Lengths(i)   = Lengths(i-1) + duraci(i);
            Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
        end
        % i. Create Load Duration Curve information output
        Power = zeros(blocks,1);
        % peak demand will be preserved
        Power(1) = max(sorted_demand);
        % All the rest of the blocks = centroids
        for i = 2:blocks
            Power(i) = Centroids(i);
        end
        % j. Correct First duration, as it is a fixed duration
        peak_excess = round((duraci(1) - round(length(demand_base)*crit_time))/(blocks-1));      % variable that will determine if we have an excess of number hours to adjust at peak demand
        duraci(1) = round(length(demand_base)*crit_time);     % correct it back to the fixed duration block
        duraci(2:end) = duraci(2:end) + peak_excess;          % we will adjust with the duration of the rest of the blocks equally
        % k. Correct the last power, to make the energy of the month stay the same
        Power(blocks) = abs((sum(sorted_demand)-Power(1:blocks-1)'*duraci(1:blocks-1)))/duraci(blocks);
        LDC = [duraci./length(current_demand),Power];
        % l. Start a counter for how many iterations it takes to reach convergence
        it = 1;
        % k. We will have a convergence criteria of iteration k - (k-1) < error, so we define both a vector k and k-1>
        LDC_k = LDC;
        LDC_k_1 = zeros(size(LDC_k,1),size(LDC_k,2));

        % m. While loop to start the clustering algorithm
        while sum(abs(LDC_k_1(:,2)-LDC_k(:,2))) >= blocks % the sum of all differences must be less than 1/20 MW times blocks
            LDC_k = LDC_k_1;
            % Update table of demand and distances to all Centroids
            for i =1:blocks
                Table(:,i) = abs(sorted_demand-Centroids(i));
            end
            
            % Assign flag to assign Centroid
            [~,Duration] = min(Table,[],2);
            % Count the duration of each Block
            [duraci,~] = groupcounts(Duration);
    
            % Update Lengths and Centroids
            Lengths(1) = duraci(1);
            Centroids(1) = mean(sorted_demand(1:Lengths(1)));
            % for loop to update blocks 2 till the finish
            for i = 2:blocks
                Lengths(i)   = Lengths(i-1) + duraci(i);
                Centroids(i) = mean(sorted_demand(Lengths(i-1)+1:Lengths(i)));
            end
    
            % Create Load Duration Curve information output
            Power = zeros(blocks,1);
            % peak demand will be preserved
            Power(1) = max(sorted_demand);
            % All the rest of the blocks = centroids except the last
            for i=2:blocks-1
                Power(i) = Centroids(i);
            end
        
            % Create Load Duration Curve information output
            Power = zeros(blocks,1);
            % peak demand will be preserved
            Power(1) = max(sorted_demand);
            % All the rest of the blocks = centroids except the last
            for i=2:blocks
                Power(i) = Centroids(i);
            end
    
            % Correct First duration
            peak_excess = round((duraci(1) - round(length(current_demand)*crit_time))/(blocks-1));      % variable that will determine if we have an excess of number hours to adjust at peak demand
            duraci(1) = round(length(current_demand)*crit_time); % correct it back to the fixed duration block
            duraci(2:end) = duraci(2:end) + peak_excess;          % we will adjust with the duration of the rest of the blocks equally
    
            % Correct last Power
            Power(blocks) = abs(sum(sorted_demand)-Power(1:blocks-1)'*duraci(1:blocks-1))/duraci(blocks);
            LDC_k_1 = [duraci./length(current_demand),Power];
            
            % Break criteria in case we have too many iterations
            if it >= 100
                break
            else
                it = it +1;
            end
        end
        % Make the last iteration equal to the current iteration
        LDC_k = LDC_k_1;
        % Convert the percentage load duration to absolute units
        LDC_k(:,1) = LDC_k(:,1)*length(current_demand);
        if sum(LDC_k(:,1)) < length(current_demand)
            LDC_k(1,1) = LDC_k(1,1) + (length(current_demand) - sum(LDC_k(:,1)));
        elseif length(current_demand) < sum(LDC_k(:,1))
            LDC_k(round(median([1:blocks])),1) = LDC_k(round(median([1:blocks])),1) + (length(current_demand) - sum(LDC_k(:,1)));
        end
        % Append the output to the master data file
        master_LDC = [master_LDC;LDC_k,[1:blocks]'];
    end
    master_LDC = array2table(master_LDC,'VariableNames',["Hours","Cluster","Block"]);
end

function [Global_LDC_Info,Global_Chronological_Info] = global_demands(Curve_LDC,demand_base)
    demand_growth = readmatrix("DemandGrowth.csv");         % load growth into a variable

    % Create Table with Info
    Global_LDC_Info = table(Curve_LDC);
    for ii = 1:length(demand_growth)
        if ii == 1      % For the first year that is not 2019, we multiply against the first growth factor
            current_demand = Curve_LDC*demand_growth(ii,2);
        else            % for the rest of the years, we multiply against the previous years we just calculated
            current_demand = current_demand*demand_growth(ii,2);
        end
        % We concatenate this into our Table
        Global_LDC_Info.(ii+1) = current_demand;
    end
    
    % We now will change the headers in our table
    % Initialize the table headers with year 2019
    headers = ["2019"];
    % For the rest we extract from the variable demand growth first column
    headers = [headers;num2str(demand_growth(:,1))];
    
    Global_LDC_Info.Properties.VariableNames = headers;
    
    % Now we repeat the procedure but for the chronological curve
    Global_Chronological_Info = table(demand_base);
    for ii = 1:length(demand_growth)
        % For the first year that is not 2019, we multiply against the first
        % growth factor
        if ii == 1
            current_demand = demand_base*demand_growth(ii,2);
        % for the rest of the years, we multiply against the previous years we
        % just calculated
        else
            current_demand = current_demand*demand_growth(ii,2);
        end
        % We concatenate this into our Table
        Global_Chronological_Info.(ii+1) = current_demand;
    end
    
    % We apply the same headers to the chronological demand
    Global_Chronological_Info.Properties.VariableNames = headers;
end

function graph_LDC(year,Curve_LDC,chronological_demand)
    global stages       % FIX: this figure should be able to graph properly if we change the stage from 1 to 12.
    if stages == 12     % In case we want to group the year in months, we will group the chronological demand into its parts.
        % Divide the entries into its given months
        % 1. Create a vector of dates that will be associated to the demand entries. Given that it is not important the year, we will do it with 1900 because it is not a leap year
        % The first entry is the 1900/JAN/1 0:00:00
        time_stamp(1,1) = datetime(1900,1,1,0,0,0);
        % The rest are just add an hour as defined 1/24
        for i = 2:length(chronological_demand)
            time_stamp(i,1) = time_stamp(i-1) + 1/24;
        end
    
        % 2. Make a vector with each month
        month_stamp = month(time_stamp);

        % 3. Append all the demand by month
        OrderedDemand = [];
        for m = 1:stages
            current_demand = sort(chronological_demand(month_stamp == m),"descend");
            OrderedDemand = [OrderedDemand;current_demand];
        end
    else
        OrderedDemand = sort(chronological_demand,"descend");
    end
    % Fig_Naam = "LDC vs Clustered LDC: " + year;
    Fig_Naam = "LDC vs Clustered LDC";
    fig = figure("Name",Fig_Naam,'NumberTitle','off','Visible',"off");
    
    hold on
    title(Fig_Naam,'FontSize',10);
    if max(max(OrderedDemand),max(Curve_LDC)) >= 20000
        OrderedDemand = OrderedDemand/1000;
        Curve_LDC = Curve_LDC/1000;
        ylabel("Power [GW]");
    else
        ylabel("Power [MW]");
    end
    plot(Curve_LDC);                            % We plot the Load Duration Curve
    plot(OrderedDemand);                        % We plot the sorted Chronological Demand
    ytickformat('%,4.0f');
    xlabel('Time (hours)');
    axis([0 inf 0 max(Curve_LDC)*1.025]);
    
    hold off
    legend('Clustered Load Duration Curve','Load Duration Curve',"Location","best");

    % Export the fig file
    Route = pwd + "\\Reports";     % Directory name
    % [imageData, alpha] = export_fig:          we disregard imageData and alpha because there is no need for it.
    [~,~] = export_fig(Route +"\LDC vs Clustered LDC " + year  +".pdf",[fig],'-append');           % In this step we export the figures to pdf.
    
    %% Close all Excel files possibly still open
    fclose('all');      close all hidden;       close all force;
end

function hourblock_map = hourblockmap(demand_base,base_LDC,opt_year)    % We will use the marker from base_LDC.Hours in which we have the duration of the block for that month.
    global blocks stages
    
% -------------------------------------------------------------%
% ----------------PREFACE--------------------------------------%
% -------------------------------------------------------------%
    switch stages
        case 12
            day_duration      = [31,28,31,30,31,30,31,31,30,31,30,31];
            month_hour_markers = zeros(12,2);
            for ii = 1:length(day_duration)
                if ii == 1                           % There is no hour zero or previous information, so we have to fill it out prior
                    month_hour_markers(ii,1) = 1;   % Initial hour of the month
                    month_hour_markers(ii,2) = month_hour_markers(ii,1) + day_duration(ii)*24 -1;       % Final hour of the month
                else
                    month_hour_markers(ii,1) = month_hour_markers(ii-1,2) + 1;   % Initial hour of the month
                    month_hour_markers(ii,2) = month_hour_markers(ii,1) + day_duration(ii)*24 -1;       % Final hour of the month
                end
            end
        case 1
            day_duration = 365;
            month_hour_markers = [1,8760]; 
        otherwise
            disp("Stages is not monthly nor annually.");
            return
    end
% -------------------------------------------------------------%
% ----------------END PREFACE----------------------------------%
% -------------------------------------------------------------%

    % 1. Create a vector for duration aggregated. This means that that specific block goes from the previous number up to the number stored for that block. Ex: 1-> 613   2-> 3986. So block 2 goes from 614-3986
    duraci_agg = zeros(stages*blocks,2);            % We can add a switch to change this. Or perhaps change it to the sum(base_LDC.Hours)
    for mm = 1:stages
        duraci_agg((mm-1)*5+1,1) = base_LDC.Hours((mm-1)*5+1,1);
        for bb = 2:blocks
            duraci_agg((mm-1)*5+bb,1) = duraci_agg((mm-1)*5+bb-1,1) + base_LDC.Hours((mm-1)*5+bb,1);
        end
    end
    
    % We will use the second column of the duraci_agg matrix to flag the blocks with the corresponding block number
    duraci_agg(:,2) = base_LDC.Block;
    
    % 2. We will create reusable vectors for comparators to assign the corresponding flag from LDC to hour
    aggregate_blk_flag = [];
    for mm = 1:stages
        blk_flag = zeros(24*day_duration(mm),1);    % this vector has the same length as hours in the stage
        blk_id = 1;                                 % We will use a marker to know in which hour of the month we are in. This will let us know in which block are still in.
        for ii = 1:length(blk_flag)
            blk_flag(ii) = duraci_agg(blk_id,2);    % We will be writing the number of the flag for increasing hours
            if ii == duraci_agg(blk_id,1)           % Until we identify we reach the last hour of that block
                blk_id = blk_id + 1;                % Then we will try to look for the following block
            end
        end
        aggregate_blk_flag = [aggregate_blk_flag;blk_flag];
    end
    
    % 4. Make a comparator for LDC and chronological order  
    % Create a time stamps vector
    % First make a first date
    year = str2double(opt_year);
    t = datetime([year,1,1]);               % Don't know why we cannot condense this into a single line
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(length(demand_base),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');
    month_stamp = month(time_stamps);

    % Now we will make the stagely comparator
    hourblock_map = [];      % The matrix we aim to fill out
    for mm = 1:stages
        % We define the current demand and its corresponding time stamps
        if stages == 12
            current_demand = demand_base(find(month_stamp == mm));
            current_time_stamps = time_stamps(find(month_stamp == mm));
        elseif stages == 1
            current_demand = demand_base;
            current_time_stamps = time_stamps;
        end
        current_blk_flag = aggregate_blk_flag(month_hour_markers(mm,1):month_hour_markers(mm,2));
        % We create an ordered matrix considering chronological time stamps and its base demand.
        A = table(current_time_stamps,current_demand);
        % Sort the table upon the demand
        A = sortrows(A,'current_demand','descend');
        % take that sorted table and the time stampos to attach it to the block flags
        B = table(A.current_time_stamps,current_blk_flag);
        % Just rename variables
        B.Properties.VariableNames = ["time_stamps","block_flag"];
        % Sort again by date in chronological manner
        C = sortrows(B,"time_stamps","ascend");

        % Append the product to our hourblock matrix
        hourblock_map= [hourblock_map;C];
    end 
end

function [renew_gen_block, renewables_gen] = RENEW_GEN_BLK(hbm)
    global blocks stages scenarios last_study_horizon fst_study_horizon

    % Get the sheet names from the Excel file
    [~, sheets] = xlsfinfo("RenewableScenario.xlsx");
    
    for s = 1:scenarios
        for t = 1:last_study_horizon - fst_study_horizon + 1
            % 0. Select random information from Renewable Scenario File
            randomSheetIndex = randi(numel(sheets));
            randomSheetName = sheets{randomSheetIndex};
            % 1. Take data as input
            renewables_gen_helper = readtable("RenewableScenario.xlsx",'ReadVariableNames',true,"UseExcel",false, 'Sheet',randomSheetName);      % The file has the variable names in the first line, so it is needed the ReadVariablesNames == true.
            RenewableScenariosName = renewables_gen_helper.Properties.VariableNames(2:end);                             % not taking into account the first one because that one is the time_stamp. Just considering the RenewableScenarios.
        
            % 2. Add hourblock flag to the table to make the average if regarding this criteria.
            renewables_gen_helper = addvars(renewables_gen_helper,hbm.block_flag,'After','time_stamp','NewVariableNames',['block_flag']);
            renewables_gen_helper = addvars(renewables_gen_helper,month(renewables_gen_helper.time_stamp),'After','block_flag','NewVariableNames',['month']);    % We add a new criteria, regarding the month    
        
            % 3. Make average if for all renewables
            renew_gen_block_helper = zeros(stages*blocks,width(RenewableScenariosName)+1);       % Matrix dimension is {stages*blocks, renewable sceanarios (Amount of Renewable Scenarios), we will also add the blocks}
            % Fill out the blocks column tag
            renew_gen_block_helper(:,1) = repmat([1:blocks]',stages,1);
        
            % Fill out the information for the average if
            for rr = 1:width(RenewableScenariosName)        % for loop counting around the amount of renewable scenarios we have.
                current_renewable_source = renewables_gen_helper.(string(RenewableScenariosName(rr)));       % we will need to add a helper variable to take out the current renewable source observed
                for mm = 1:stages             % For loop around stages (could be months or annually)
                    for bb = 1:blocks                   
                        renew_gen_block_helper((mm-1)*blocks+bb,rr+1) = mean(current_renewable_source(and(renewables_gen_helper.block_flag == bb, renewables_gen_helper.month == mm)));
                    end                         % Reads: average(renewable_source(for the entries == block || == month))
                end
            end
        
            % Convert renewewables generation block to table
            renew_gen_block_helper = array2table(renew_gen_block_helper);
            renew_gen_block_helper.Properties.VariableNames = [['Block';RenewableScenariosName']];
        
            % Append information to struct
            Year_Name = "Y" + num2str(fst_study_horizon + t -1);
            Scenario_Name = "S" + s;
            % First for the Hourly File:
            renewables_gen.(Scenario_Name).(Year_Name) = renewables_gen_helper;
            % Now for the Block File:
            renew_gen_block.(Scenario_Name).(Year_Name) = renew_gen_block_helper;
        end
    end
end

function GEP_optimizer = GEP_optimizer_creator(GEP_Load,Renew_commitment,renew_gen_block,DBB_info,opt_year)
    global blocks rate Curt_penalty fst_study_horizon stages scenarios last_study_horizon

    %% Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Fuel_price_data = DBB_info.Fuel_price_data;
    Flags = DBB_info.Flags;

    %% General Information for the whole mathematical problem. This data shouldn't change across the years.
    % Note: The renewable information might change in the future. As of now, the problem will be deterministic and the generation of
    % renewable will be the same across all years. However, to make it stochastic, we will need to create a function that chooses how
    % renewables will generate energy across each year.

    % Existing Generator Bounds
    Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
    
    % Candidate Generator Bounds
    Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit; ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
    
    % Max Installation Bounds: Given that decision variable is integer, we set
    % a limit on the total install capacity coming from each source
    MaxInstall = [Thermal_cand_data.MaxInstall; EcoThermal_cand_data.MaxInstall; Renewables_cand_data.MaxInstall; ESS_cand_data.MaxInstall];
    
    % Load
    PLoad = GEP_Load{:,2:end}; % Power Instances will be equal to the representing blocks
    
    % Renewable generation quota
    % Peaking Information: According to information from paper couple long-term
    % and short-term
    vCap_exist= Exist_Upper;    % the Capacity of Installed plants is the same as the Upper Limit of Existing plants
    % Peak Factor coefficient factor contribution for existing generators
    % Renewables Peak Factor coefficient factor
    Exist_renewables_capacity_factor = ones(stages*blocks,height(Renewables_exist_data));
    for s = 1:scenarios
        for t = 1: last_study_horizon - fst_study_horizon +1
            for i = 1:height(Renewables_exist_data)
                if t == 1
                    Exist_renewables_capacity_factor(:,i) = [renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_exist_data.RenewableScenario{i})];
                else
                    Exist_renewables_capacity_factor(:,i) = mean([Exist_renewables_capacity_factor(:,i),renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_exist_data.RenewableScenario{i})],2);
                end
            end
        end
    end
    % Take out the capacity factor of renewables for peak demand
    Exist_renewables_pPF = zeros(1,height(Renewables_exist_data));
    for i = 1:height(Renewables_exist_data)
        for mm = 1:stages
            if mm == 1
                Exist_renewables_pPF(1,i) = Exist_renewables_capacity_factor((mm-1)*blocks+1,i);
            else
                Exist_renewables_pPF(1,i) = mean([Exist_renewables_pPF(1,i),Exist_renewables_capacity_factor((mm-1)*blocks+1,i)]);
            end
        end
    end
    % Append all peak contribution factors
    pPF_exist = [Thermal_exist_data.HistoricalAvailability;EcoThermal_exist_data.HistoricalAvailability;Exist_renewables_pPF']; % Peak factor of all technologies in peak demand.
    
    % Peak Factor coefficient factor contribution for candidate generators
    Cand_renewables_capacity_factor = ones(stages*blocks,height(Renewables_cand_data));
    for s = 1:scenarios
        for t = 1: last_study_horizon - fst_study_horizon +1
            for i = 1:height(Renewables_cand_data)
                if t == 1
                    Cand_renewables_capacity_factor(:,i) = [renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_cand_data.RenewableScenario{i})];
                else
                    Cand_renewables_capacity_factor(:,i) = mean([Cand_renewables_capacity_factor(:,i),renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_cand_data.RenewableScenario{i})],2);
                end
            end
        end
    end
    % Take out the capacity factor of renewables for peak demand
    Cand_renewables_pPF = zeros(1,height(Renewables_cand_data));
    for i = 1:height(Renewables_cand_data)
        for mm = 1:stages
            if mm == 1
                Cand_renewables_pPF(1,i) = Cand_renewables_capacity_factor((mm-1)*blocks+1,i);
            else
                Cand_renewables_pPF(1,i) = mean([Cand_renewables_pPF(1,i),Cand_renewables_capacity_factor((mm-1)*blocks+1,i)]);
            end
        end
    end
    % Append all peak contribution factors
    pPF_cand = [Thermal_cand_data.HistoricalAvailability;EcoThermal_cand_data.HistoricalAvailability;Cand_renewables_pPF';ones(height(ESS_cand_data),1)]; % Peak factor of all technologies in peak demand.        We define the factor of Batteries as able to deliver all of their installed capacity.
    
    %% Mathematical Formulation
    % Decision variables
    % Investment
    Exist_xit      = binvar(Flags.quant_exist,width(PLoad),'full');                    % Binary variables to keep track the installed technologies.
    Exist_yit      = intvar(Flags.quant_exist,width(PLoad),'full');                    % Integer variable to keep track of the installed technologies.
    xit            = intvar(Flags.quant_cand,width(PLoad),'full');                     % Integer variables to define the installed capacity needed. Variable is 2-D: [quantity of type of technology decision, years]. The decision variables doesn't change per block.
    yit            = intvar(Flags.quant_cand,width(PLoad),'full');                     % Integer variable used to measure Investment. Not a decision variable. Dependent variable.
    wit            = intvar(Flags.quant_cand,width(PLoad),'full');                     % Integer variables used to measure Amortization variables 
    alpha_it       = intvar(Flags.quant_cand,width(PLoad),'full');                     % Integer variable to measure the early retirement. 
    % Traditional Operation
    Exist_git      = sdpvar(Flags.quant_exist,stages*blocks,width(PLoad),scenarios,'full');             % continuous variable to define generation output from exist generation. Variable is 3-D, for same reason as Cand_git.
    Cand_git       = sdpvar(Flags.flag_renewables_cand,stages*blocks,width(PLoad),scenarios,'full');    % continuous variable to define generation output from future generation. Variables is 3-D: [quantity of candidate generators, block, year]. There should be an independent value for each gen i, block b, and year t.
    % Battery operation
    q_Charge_itb   = sdpvar(height(ESS_cand_data),stages*blocks,width(PLoad),scenarios,'full');         % continuous variable to define discharge from candidate batteries.
    q_Discharge_itb = sdpvar(height(ESS_cand_data),stages*blocks,width(PLoad),scenarios,'full');        % continuous variable to define discharge from candidate batteries.
    d_itb          = binvar(height(ESS_cand_data),stages*blocks,width(PLoad),scenarios,'full');         % binary variable to constrain batteries to only charge or only discharge
    % Decision variables for curtailment
    Curt_exist_git = sdpvar(Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,stages*blocks,width(GEP_Load)-1,scenarios,'full');             % continuous variable to define curtailment from renewables. Variable is 3-D, for same reason as Exist_git.
    Curt_cand_git = sdpvar(Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,stages*blocks,width(GEP_Load)-1,scenarios,'full');             % continuous variable to define curtailment from renewables. Variable is 3-D, for same reason as Exist_git.

    % Nomenclature:
    % i: Generators. Will be used both for aggregated generators in Existing or Candidate decision variable object.
    % b: blocks.
    % t: years.
    
    % Parameters: parameters that make certain constraints change over time.
    sOR             = sdpvar(1,width(PLoad),'full');       % System Operating Reserve parameter
    ramp_capacity   = sdpvar(1,width(PLoad),'full');       % ramping capacity for each year.
    Renew_Cap       = sdpvar(1,width(PLoad),'full');       % Renewable Energy Cap for each year to prevent too much curtailment.

    %% Constraints
    % General
    Constraints = [];
    % Generator Bounds: Min and Max Constraints
    % Existing Generator Bounds
    for s = 1:scenarios
        for t = 1:width(PLoad)
                for b = 1:stages*blocks     
                    Constraints = Constraints + [(0 <= Exist_git(:,b,t,s) <= Exist_xit(:,t).*Exist_Upper) : 'Existing Generator Bounds'];       % given that we are not considering on and off conditions for the machines, the lower bound is zero.
                end
        end
    end
    % Candidate Generator Bounds: All technologies including batteries
    for s = 1:scenarios
        for t = 1:width(PLoad)
            for b = 1:stages*blocks                % vector of zeros size as long as Git_cand + Batteries <= [Git_cand;Batteries]              <= decisionvariable*UpperLimit
                Constraints = Constraints + [(0 <= Cand_git(1:Flags.flag_renewables_cand,b,t,s) <= xit(1:Flags.flag_renewables_cand,t).*Cand_Upper(1:Flags.flag_renewables_cand)): 'Candidate Generator Bounds']; % given that we are not considering on and off conditions for the machines, the lower bound is zero.    % FIX: We should change the CandUpper variable
                % Given that iterates over the same loops, we will add the battery charge and discharge constraint in here
                Constraints = Constraints + [(0 <= q_Discharge_itb(:,b,t,s) <= d_itb(:,b,t,s).*xit(Flags.flag_renewables_cand+1:end,t).*ESS_cand_data.UpperLimit): 'Battery discharge bounds'];    % Discharge constraint
                Constraints = Constraints + [(0 <= q_Charge_itb(:,b,t,s) <= ((1-d_itb(:,b,t,s)).*xit(Flags.flag_renewables_cand+1:end,t)).*ESS_cand_data.UpperLimit): 'Battery Charge Bounds'];       % charge constraint
            end                            % vector of zeros size as big as Batteries <=  Batteries charge <= decision variable*UpperLimit
        end
    end
    % Power Balance Constraint
    for s = 1:scenarios
        for t = 1:width(PLoad)
            for b = 1:stages*blocks
                Constraints = Constraints + [(sum(Exist_git(:,b,t,s)) + sum(Cand_git(:,b,t,s)) + sum(q_Discharge_itb(:,b,t,s)) == PLoad(b,t) + sum(q_Charge_itb(:,b,t,s))): 'Power Balance Constraint'];          % Fix batteries in here
            end     % This constraint reads: Sum(Exist_Git) + Sum(Candi_Git)              + batteries discharge         = Load for that specific block + batteries charge
        end
    end
    % Max Installation constraint: perhaps we cannot build more than a certain amount of technology given area constraints
    for t = 1:width(PLoad)
        Constraints = Constraints + [[xit(:,t).*Cand_Upper <= MaxInstall] : 'Max Install'];
    end

    % Peaking Equation: Firm Power Security
    for t = 1:width(PLoad)
        Constraints = Constraints + [[(pPF_exist'*(Exist_xit(:,t).*vCap_exist) + (xit(:,t).*pPF_cand)'*Cand_Upper) >= max(PLoad(:,t))*(1+sOR(t))]: 'Peaking Equation'];          %sOR: system Operating Reserve
    end                             % Reserve * (pPF_exist*InstalledCap + decision variable*pPF_cand*MaxCapacity) >= Peak demand at each year.

    % Ramping Capacity:
    for t = 1:width(PLoad)     % Ramp capacity at each year must be at least the RampUp of Existing Thermal, EcoThermal + RampUp of Candidate Thermal,EcoThermal + Installed Capacity of Batteries. We consider the installed capacity of Batteries because they can act immediately.
        Constraints = Constraints + [(ramp_capacity(t) <= Thermal_exist_data.RampUp'*Exist_xit(1:Flags.flag_thermal_exist,t) + EcoThermal_exist_data.RampUp'*Exist_xit(Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist,t) + Thermal_cand_data.RampUp'*xit(1:Flags.flag_thermal_cand,t) + EcoThermal_cand_data.RampUp'*xit(Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand,t) + ESS_cand_data.UpperLimit'*xit(Flags.flag_renewables_cand+1:Flags.flag_ESS_cand,t)): 'Ramping constraint'];
    end

    % Cap on the amount of energy that can come from Renewables
    for s = 1:scenarios
        for t = 1:width(PLoad)     % sum(Exist_Upper*AvailableResource) + sum(xit*Exist_Upper*AvailableResource) <= Renew_Cap
            LHS = 0;
            % Existing Renewables
            counter = 1;        
            for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
                LHS = LHS + Exist_xit(i,t)*sum((Exist_Upper(i)*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_exist_data.RenewableScenario{counter}))).*(GEP_Load.LD));
                counter = counter + 1;
            end
            % Candidate Renewables
            counter = 1;        
            for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
                LHS = LHS + xit(i,t)*sum((Cand_Upper(i)*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_cand_data.RenewableScenario{counter}))).*(GEP_Load.LD));
                counter = counter + 1;
            end
            % Actual Constraint
            Constraints = Constraints + [LHS <= Renew_Cap(t)];
        end
    end
    
    %% Renewables generation resource: This constraint is read that the generation from renewables needs to be at most the available resource
    % Existing renewables
    for s = 1:scenarios
        for t = 1:width(PLoad)    
                counter = 1; % fictitious variable to start the count for the length of Solar
                for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist   %% i starts where the thermal generators ends, and for an amount of solar generators as defined by the file Exist_Solar.
                    Constraints = Constraints + [(Exist_git(i,:,t,s)' <= Exist_xit(i,t)*Exist_Upper(i).*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_exist_data.RenewableScenario{counter}))): 'Existing Renewable resource'];
                    counter = counter + 1;      % This constraint reads: Git <= InstalledCapacity*resource for that time and block.
                end                             % This constraints compare a vector with another vector directly.
        end
    end    
    % Candidate renewables
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1; % fictitious variable to start the count for the length of Solar
            for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand     %% i starts where thermal generators end, and for an amount of solar as defined in Cand_Solar.
                Constraints = Constraints + [(Cand_git(i,:,t,s)' <= xit(i,t)*Cand_Upper(i).*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_cand_data.RenewableScenario{counter}))): 'Existing Renewable resource'];
                counter = counter + 1;      % This constraint reads: Git <= InstalledCapacity*resource for that time and block.
            end
        end
    end

    % Add constraint to make solar and wind to not grow in an uneven way
    % Find where the data for wind starts.
    first_wind_exist = find(Renewables_exist_data.Tech == "Wind",1);
    first_wind_cand = find(Renewables_cand_data.Tech == "Wind",1);  
    for s = 1:scenarios
        for t = 5:width(PLoad)
            Constraints = Constraints + [0 <= sum(Exist_git(Flags.flag_EcoThermal_exist+1:Flags.flag_EcoThermal_exist+first_wind_exist-1,:,t,s)*(GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+1:Flags.flag_EcoThermal_cand+first_wind_cand-1,:,t,s)*(GEP_Load.LD)) <= 3*(sum(Exist_git(Flags.flag_EcoThermal_exist+first_wind_exist:Flags.flag_renewables_exist,:,t,s)*(GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+first_wind_cand:Flags.flag_renewables_cand,:,t,s)*(GEP_Load.LD)))];
            Constraints = Constraints + [0 <= sum(Exist_git(Flags.flag_EcoThermal_exist+first_wind_exist:Flags.flag_renewables_exist,:,t,s)*(GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+first_wind_cand:Flags.flag_renewables_cand,:,t,s)*(GEP_Load.LD)) <= 3*(sum(Exist_git(Flags.flag_EcoThermal_exist+1:Flags.flag_EcoThermal_exist+first_wind_exist-1,:,t,s)*(GEP_Load.LD)) + sum(Cand_git(Flags.flag_EcoThermal_cand+1:Flags.flag_EcoThermal_cand+first_wind_cand-1,:,t,s)*(GEP_Load.LD)))];
        end
    end
    % Renewable quota contraint. It reads: Sum(Cand_Git(StartAtEcoThermal:end,:,final year)*                          *(vector of LoadDuration)) + Sum(Exist_Git(StartAtEcoThermak:end,:,final year)*                                 *(vector of LoadDuration)) >= Renewable commitment*LoadDuration*Clusters;
    for s = 1:scenarios
        Constraints = Constraints + [[sum(Cand_git(Flags.flag_thermal_cand+1:Flags.flag_renewables_cand,:,str2double(opt_year)-fst_study_horizon+1,s)*(GEP_Load.LD)) + sum(Exist_git(Flags.flag_thermal_exist+1:Flags.flag_renewables_exist,:,str2double(opt_year)-fst_study_horizon+1,s)*(GEP_Load.LD)) >= Renew_commitment*GEP_Load.LD'*GEP_Load.(opt_year)]: 'Renewable Quota']; 
    end

    % Measure Renewable Curtailment
    % Non negativity of curtailment
    Constraints = Constraints + [Curt_exist_git >= 0];
    Constraints = Constraints + [Curt_cand_git >= 0];
    % Existing renewables
    for s = 1:scenarios
        for t = 1:width(PLoad)    
            counter = 1;                % We create this counter because the index i is different when pointing in the Renewable Scenario list, than when pointing towards the position inside the Exist_Upper and Exist_git matrices.
            for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
                Constraints = Constraints + [[Curt_exist_git(counter,:,t,s)' == Exist_xit(i,t)*Exist_Upper(i).*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_exist_data.RenewableScenario{counter})) - Exist_git(i,:,t,s)']: 'Existing Renewable Curtailment'];
                counter = counter + 1;      % Constraint reads: Curtailment == Renewable Energy Available - Actual Renewable Energy Used
            end
        end
    end
    %Candidates
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;                % We create this counter because the index i is different when pointing in the Renewable Scenario list, than when pointing towards the position inside the Cand_Upper and Cand_git matrices.
            for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
                Constraints = Constraints + [[Curt_cand_git(counter,:,t,s)' == xit(i,t)*Cand_Upper(i).*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_cand_data.RenewableScenario{counter})) - Cand_git(i,:,t,s)']: 'Candidate Renewable Curtailment'];
                counter = counter + 1;      % Constraint reads: Curtailment == Renewable Energy Available - Actual Renewable Energy Used
            end
        end
    end

    % Maximum Curtailment per year
    % Existing Generators
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
                % This constraint reads:    sum(Curtailment of Renewable i at year t) <= 0.05*Installed Capacity*Actual resource. So, we are saying that we can only curtail 5% of the resource at Planning stage.
                Constraints = Constraints + [sum(Curt_exist_git(counter,:,t,s)) <= 0.03*sum(Exist_Upper(i)*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_exist_data.RenewableScenario{counter})))];
                counter = counter + 1;
            end
        end
    end
    % Candidate Generators
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
                Constraints = Constraints + [sum(Curt_cand_git(counter,:,t,s)) <= 0.03*sum(xit(i,t)*Cand_git(i)*renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(string(Renewables_cand_data.RenewableScenario{counter})))]; 
            end
        end
    end

    %% Minimum year for investments: All Candidates
    % Thermal
    for t = 1:width(PLoad)
        for i = 1:Flags.flag_thermal_cand
            if Thermal_cand_data.MinYear(i) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'Thermal minimum year constraint'];
            end
        end
    end
    % EcoThermal
    for t = 1:width(PLoad)
        counter = 1;
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            if EcoThermal_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'EcoThermal minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    % Renewables
    for t = 1:width(PLoad)
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            if Renewables_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'Renewables minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    % Energy Storage Systems
    for t = 1:width(PLoad)
        counter = 1;
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
            if ESS_cand_data.MinYear(counter) > fst_study_horizon +t -1
                Constraints = Constraints + [(xit(i,t) == 0): 'ESS minimum year constraint'];
            end
            counter = counter +1;
        end
    end
    
    %% Battery Additional Constraints
    % Battery Energy Conservation Constraint: Annually
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand % the charged and discharged energy decision variables do not need a counter flag because it only is defined by the amount of batteries in the system
                Constraints = Constraints + [(q_Discharge_itb(counter,:,t,s)*GEP_Load.LD == ESS_cand_data.Efficiency(counter)*(q_Charge_itb(counter,:,t,s)*GEP_Load.LD)): 'Batteries Energy conservation constraint annually'];
                counter = counter + 1;
            end
        end
    end
    % Maximum energy discharged and charged per year
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
                Constraints = Constraints + [(q_Discharge_itb(counter,:,t,s)*GEP_Load.LD <= 365*xit(i,t)*ESS_cand_data.Duration(counter)*ESS_cand_data.UpperLimit(counter)): 'Maximum energy discharged per year'];
                Constraints = Constraints + [(q_Charge_itb(counter,:,t,s)*GEP_Load.LD <= 365*xit(i,t)*ESS_cand_data.Duration(counter)*ESS_cand_data.UpperLimit(counter)/ESS_cand_data.Efficiency(counter)): 'Maximum energy Charged per year'];
                counter = counter + 1;      % Discharge*LoadDuration ([MW]*[h])   <=  365[days]*InstalledQuantity[n/a]*Duration[h/day]*UpperLimit[MW]  -> [MWh]
            end
        end
    end

    % Battery Energy Conservation Constraint: Monthly
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand % the charged and discharged energy decision variables do not need a counter flag because it only is defined by the amount of batteries in the system
                for b = 1:blocks:stages*blocks
                    Constraints = Constraints + [(q_Discharge_itb(counter,b:b+blocks-1,t,s)*GEP_Load.LD(b:b+blocks-1) == ESS_cand_data.Efficiency(counter)*(q_Charge_itb(counter,b:b+blocks-1,t,s)*GEP_Load.LD(b:b+blocks-1))): 'Batteries Energy conservation constraint monthly'];
                end
                counter = counter + 1;
            end
        end
    end
    % Maximum energy discharged and charged per month
    for s = 1:scenarios
        for t = 1:width(PLoad)
            counter = 1;
            for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
                for b = 1:blocks:stages*blocks
                    Constraints = Constraints + [(q_Discharge_itb(counter,b:b+blocks-1,t,s)*GEP_Load.LD(b:b+blocks-1) <= sum(GEP_Load.LD(b:b+blocks-1,:))/24*xit(i,t)*ESS_cand_data.Duration(counter)*ESS_cand_data.UpperLimit(counter)): 'Maximum energy discharged per month']; 
                    Constraints = Constraints + [(q_Charge_itb(counter,b:b+blocks-1,t,s)*GEP_Load.LD(b:b+blocks-1) <= sum(GEP_Load.LD(b:b+blocks-1,:))/24*xit(i,t)*ESS_cand_data.Duration(counter)*ESS_cand_data.UpperLimit(counter)/ESS_cand_data.Efficiency(counter)): 'Maximum energy Charged per month'];
                end
                counter = counter + 1;
            end
        end
    end

    %% Investment related constraints
    % Decision related to existing technologies:
    Constraints = Constraints + [Exist_xit(:,1) == 1];      % Existing technologies start their life as connected and decided.
    % Measure Exist_yit as the difference between consecutive years
    for t = 1:width(Exist_xit)
        if t == 1
            Constraints = Constraints + [Exist_yit(:,1) == Exist_xit(:,1)];
        else
            Constraints = Constraints + [Exist_yit(:,t) == Exist_xit(:,t) - Exist_xit(:,t-1)];
        end
    end
    
    % Constraint Exist_yit to retire the plant only once.
    for i = 1:height(Exist_xit)
        Constraints = Constraints + [0 <= sum(Exist_yit(i,:)) <= 1];
    end

    % Difference between investments needed: yit
    for t = 1:width(xit)
        if t == 1
            Constraints = Constraints + [(yit(:,1) == xit(:,1)): 'Accumulated and New Investment equality for first year'];   % For the first year, Investment and math_indicator are the same.
        else
            Constraints = Constraints + [(yit(:,t) == xit(:,t)-xit(:,t-1)): 'Investment difference between consecutive years']; % However, for the rest, is the difference between the capacity needed in the previous year
        end                                                                                         % And this year. In case there is no new capacity added, then the incurred cost should be zero.
    end

    % Measure uninstalling infrastructure: alpha
    for t = 1:width(PLoad)
        % Constraints = Constraints + [alpha_it(i,t) == min(0,yit(i,t))];
        Constraints = Constraints + [alpha_it(:,t) <= 0; alpha_it(:,t) <= yit(:,t)];        % The original equation is the one above. We linearize it in this constraint
    end                             % Alpha needs to be equal or less than the min value between 0 or yit.                     

    % Annuities costs: wit
    % Candidate Thermal Generators
    for t = 1:width(PLoad)
         for i = 1:Flags.flag_thermal_cand
             rolling_window = max(1,t-Thermal_cand_data.Amortization(i)+1):t;    % This define the window of time for when the new investment is valid
             Constraints = Constraints + [wit(i,t) == sum(yit(i,rolling_window) - alpha_it(i,rolling_window))];   % The investment is equal to the active payments on technology i. This include decision variables that were made even when retired.
         end                             
    end
    % Candidate EcoThermal Generators
    for t = 1:width(GEP_Load)-1
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
             rolling_window = max(1,t-EcoThermal_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             Constraints = Constraints + [wit(i,t) == sum(yit(i,rolling_window) - alpha_it(i,rolling_window))];   % The investment is equal to the active payments on technology i. This include decision variables that were made even when retired.
             counter = counter +1;      % wit must be more than yit and more than the rolling window sum of yit.    
         end
    end
    % Candidate Renewable Generators
    for t = 1:width(PLoad)
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
             rolling_window = max(1,t-Renewables_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             Constraints = Constraints + [wit(i,t) == sum(yit(i,rolling_window) - alpha_it(i,rolling_window))];   % The investment is equal to the active payments on technology i. This include decision variables that were made even when retired.
             counter = counter +1;      
         end
    end
    % Batteries
    for t = 1:width(PLoad)
        counter = 1;       % We keep this counter because the technology i which we are pointing in the vector Amortization is probably different than the in wit or yit
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
             rolling_window = max(1,t- ESS_cand_data.Amortization(counter)+1):t;    % This define the window of time for when the new investment is valid
             Constraints = Constraints + [wit(i,t) == sum(yit(i,rolling_window) - alpha_it(i,rolling_window))];   % The investment is equal to the active payments on technology i. This include decision variables that were made even when retired.
             counter = counter +1;      
         end
    end

    %% Objective Function  
    % Investment Cost
    % Information is drawn from: Capital Cost and Performance Characteristic Estimates for Utility Scale Electric Power Generating Technologies:
    % https://www.eia.gov/analysis/studies/powerplants/capitalcost/pdf/capital_cost_AEO2020.pdf
    Cand_invest_Cost = [Thermal_cand_data.InvestmentCost;EcoThermal_cand_data.InvestmentCost;Renewables_cand_data.InvestmentCost;ESS_cand_data.InvestmentCost];
    
    % Initialize Cost
    Cost = 0;
    % Fuel Costs
    % Add Operative Cost - Existing Thermal
    pi_s = 1/scenarios;         % We declare the probability of each scenario
    for s = 1:scenarios
        OC_Cost = 0;            % Initialize the OC Cost at zero.
        for t = 1:str2double(opt_year)-fst_study_horizon+1
            for i = 1:Flags.flag_thermal_exist
                OC_Cost = OC_Cost + Fuel_price_data.(string(Thermal_exist_data.Fuel{i}))(num2str(fst_study_horizon+t-1))*Thermal_exist_data.VarCost(i).*(Exist_git(i,:,t)*GEP_Load.LD)/((1+rate)^(t-1));
            end               % This formula reads: [AnnualVariation]*                                            *[ExistVariableCost(USD/MWh)]* *[Production(MW)]*[Duration(h)]/(1+r)^t;
        end
        % Add variable cost - Existing EcoThermal           
        for t = 1:str2double(opt_year)-fst_study_horizon+1
            counter = 1;
            for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
                if not(string(EcoThermal_exist_data.Fuel{counter}) == "Geothermal")   %% I skip the cost of Geothermal given that the renewable production cost of Geothermal is zero.
                    OC_Cost = OC_Cost + Fuel_price_data.(string(EcoThermal_exist_data.Fuel{counter}))(num2str(fst_study_horizon+t-1))*EcoThermal_exist_data.VarCost(counter).*(Exist_git(i,:,t)*GEP_Load.LD)/((1+rate)^(t-1)); 
                end               % This formula reads: [AnnualVariation]*                                                     *[ExistVariableCost(USD/MWh)]*          *[Production(MW)]*[Duration(h)]/(1+r)^t;
                counter = counter + 1;
            end
        end
        % Add variable cost - Candidate Thermal
        for t = 1:str2double(opt_year)-fst_study_horizon+1
            for i = 1:Flags.flag_thermal_cand
                OC_Cost = OC_Cost + Fuel_price_data.(string(Thermal_cand_data.Fuel{i}))(num2str(fst_study_horizon+t-1))*Thermal_cand_data.VarCost(i).*(Cand_git(i,:,t)*GEP_Load.LD)/((1+rate)^(t-1));
            end               % This formula reads: [AnnualVariation]*                                           *[ExistVariableCost(USD/MWh)]**[Production(MW)]*[Duration(h)]/(1+r)^t;
        end
        % Add variable cost - Candidate EcoThermal
        for t = 1:str2double(opt_year)-fst_study_horizon+1
            counter = 1;
            for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
                if not(string(EcoThermal_cand_data.Fuel{counter}) == "Geothermal")        %% I skip the cost of Geothermal given that the renewable production cost of Geothermal is zero.
                    OC_Cost = OC_Cost + Fuel_price_data.(string(EcoThermal_cand_data.Fuel{counter}))(num2str(fst_study_horizon+t-1))*EcoThermal_cand_data.VarCost(counter).*(Cand_git(i,:,t)*GEP_Load.LD)/((1+rate)^(t-1)); 
                end               % This formula reads: [AnnualVariation]*                                                    *[ExistVariableCost(USD/MWh)]*         *[Production(MW)]*[Duration(h)]/(1+r)^t;
                counter = counter + 1;
            end
        end
        
        % Curtailment Penalties: Both done at the same time.
        for t = 1:str2double(opt_year)-fst_study_horizon+1
            OC_Cost = OC_Cost + Curt_penalty*sum(Curt_exist_git(:,:,t,s)*GEP_Load.LD)/((1+rate)^(t-1));  % FIX: The curtailment cost is fixed to 100. But, I should change it for it to be a variable. This comment applies all over the script. Especially at the UCP.
            OC_Cost = OC_Cost + Curt_penalty*sum(Curt_cand_git(:,:,t,s)*GEP_Load.LD)/((1+rate)^(t-1));      
        end        % Penalty = sum(Curtailment at any given year)*PenaltyCost/(1+rate)^year

        Cost = Cost + pi_s*OC_Cost;            % We add the weighted Operative Cost of each scenario.
    end
    
    % O&M for Existing Plants
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        % Thermal
        for i = 1:Flags.flag_thermal_exist
            Cost = Cost + Exist_Upper(i)*Exist_xit(i,t)*Thermal_exist_data.FixedOM(i)/((1+rate)^(t-1));
        end
        % EcoThermal
        counter = 1;
        for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
            Cost = Cost + Exist_Upper(i)*Exist_xit(i,t)*EcoThermal_exist_data.FixedOM(counter)/((1+rate)^(t-1));
            counter = counter +1;
        end
        % Renewables
        counter = 1;
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            Cost = Cost + Exist_Upper(i)*Exist_xit(i,t)*Renewables_exist_data.FixedOM(counter)/((1+rate)^(t-1));
            counter = counter +1;
        end
    end
    
    % Investment Costs and O&M for Candidate Plants
    % Add investment cost: Candidate Thermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        for i = 1:Flags.flag_thermal_cand
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((rate*(1+rate)^Thermal_cand_data.Amortization(i))/((1+rate)^Thermal_cand_data.Amortization(i)-1))/((1+rate)^(t-1));
            % Interpretration: Cost = vCap*xit*O&M/(1+rate)^year
            Cost = Cost + Cand_Upper(i)*xit(i,t)*Thermal_cand_data.FixedOM(i)/((1+rate)^(t-1));
        end 
    end
    % Add investment cost: Candidate EcoThermal
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((rate*(1+rate)^EcoThermal_cand_data.Amortization(counter))/((1+rate)^EcoThermal_cand_data.Amortization(counter)-1))/((1+rate)^(t-1));
            % Interpretration: Cost = vCap*xit*O&M/(1+rate)^year
            Cost = Cost + Cand_Upper(i)*xit(i,t)*EcoThermal_cand_data.FixedOM(counter)/((1+rate)^(t-1));
            counter = counter + 1;
        end 
    end
    % Add investment cost: Candidate Renewables
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((rate*(1+rate)^Renewables_cand_data.Amortization(counter))/((1+rate)^Renewables_cand_data.Amortization(counter)-1))/((1+rate)^(t-1));
            % Interpretration: Cost = vCap*xit*O&M/(1+rate)^year
            Cost = Cost + Cand_Upper(i)*xit(i,t)*Renewables_cand_data.FixedOM(counter)/((1+rate)^(t-1));
            counter = counter + 1;
        end 
    end
    % Add investment cost: ESS
    for t = 1:str2double(opt_year)-fst_study_horizon+1
        counter = 1;
        for i = Flags.flag_renewables_cand+1:Flags.flag_ESS_cand
            % Interpretation:    Cost = sum(Installed*Wit*UnitInvestmentCost*(rate(1+rate)^n)/((1+rate)^n-1)))/(1+rate)^year,       where n: is the amortization years of that machine and it might differ from t.
            Cost = Cost + (Cand_Upper(i)*wit(i,t))*Cand_invest_Cost(i)*((rate*(1+rate)^ESS_cand_data.Amortization(counter))/((1+rate)^ESS_cand_data.Amortization(counter)-1))/((1+rate)^(t-1));
            % Interpretration: Cost = vCap*xit*O&M/(1+rate)^year
            Cost = Cost + Cand_Upper(i)*xit(i,t)*ESS_cand_data.FixedOM(counter)/((1+rate)^(t-1));
            counter = counter + 1;
        end 
    end
    
    % Penalties
    % Early Retirement penalty factor: 10% of the investment costs
    for t = 1:width(PLoad)
        for i=1:Flags.quant_cand
            Cost = Cost - 0.1*Cand_Upper(i)*alpha_it(i,t)*Cand_invest_Cost(i)/((1+rate)^(t-1));
        end             % We add the negative of alpha because in case of a uninstallation, alpha is negative. But the cost should be positive.
    end
    

    %% Set options for YALMIP
    %                    verbose: amount of output information; % Our OF has big parameters so we increase Numeric Focus; % We set a time limit of 20 minutes and a Gap of 1.5%.   
    options = sdpsettings('verbose',2,'solver','gurobi','debug',1,'gurobi.NumericFocus',2,'gurobi.TimeLimit',1800,'gurobi.MIPgap',0.01,'gurobi.Seed',rand()*1000);    % We try to add a random seed variable to check if we find alternative solutions.
    options = sdpsettings(options, 'usex0',1, 'gurobi.Heuristics',0.15,'gurobi.ImproveStartGap',0.03,'gurobi.MIPFocus',3);
    %                     % ,'gurobi.ImproveStartGap',0.01,'gurobi.MIPFocus',3                    
    % Let gurobi try to warm start the solution; %% Double the amount of time invested in Heuristics;% After we reach 1% optimality gap, we can look
    %                     % for better feasible solutions. % MIPFocus = 3 because the bound is not moving fast enough.
    

    % We don't solve the problem anymore, but instead create the optimizer object:
    GEP_optimizer = optimizer(Constraints,Cost,options,{sOR,ramp_capacity,Renew_Cap},{xit,Cand_git,q_Discharge_itb,q_Charge_itb,Exist_git,Cost});
end

function [Investment,Output_G_exist,Output_G_cand,Cost,GEP_optimizer] = GEP_evaluate(GEP_optimizer,sOR,ramp_capacity,DBB_info,Renew_Cap,iterations)
    import java.text.*;     v = DecimalFormat;                         % To use when printing results for easier readibility.

    % First: Evaluate the GEP Optimizer:
    [sol,errorcode,~,~,GEP_optimizer] = GEP_optimizer({sOR,ramp_capacity,Renew_Cap});
    
    % Second: call_BDD info to display better the information in Command Window and create a vector with the years
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;

    Name_cand  = [Thermal_cand_data.Name;EcoThermal_cand_data.Name;Renewables_cand_data.Name;ESS_cand_data.Name];       % Candidate Names

    ColNaams = ["2021"];
    for i = 2:30
        ColNaams(i) = [convertCharsToStrings(num2str(2021-1+i))];
    end


    % Third: we will test the errorcode. In case there is no error we will convert the information to be ready to graph
    if errorcode == 0 || errorcode == 3         % Yalmip sees Gurobi running out of time as an error. We make a quickfix by including this errorcode as an admissible number.
        disp("The accumulated investment decision variable Xit: " + newline);
        
        disp(array2table(sol{1},'RowNames',Name_cand,'VariableNames',ColNaams));
        disp("The total cost of this GEP is: " + char(v.format(sol{6})) + " USD" + newline)   % 6: Cost 
        Cost = sol{6};
        % The variable Output_G have four dimensions [tech,cluster,year,scenario], so we make a mean over the scenarios for the graphs.
        Output_G_exist = mean(sol{5},4);                                   % 5: Exist_git
        Output_G_cand = mean([sol{2};sol{3};-sol{4}],4);                   % 2: Cand_git; 3: q_Discharge_ih; 4: q_Charge_ih
        % The variable Investment only has two dimensions, which are technology and year.
        Investment      = sol{1};                                          % 1: xit
        % Print Investment to Excel File
        Excel_Investment = array2table(Investment,'RowNames',Name_cand,'VariableNames',ColNaams);
        filename = pwd + "\Reports\" +"GEP_Investment_results.xlsx";
        sheet = "Iteration_" + num2str(iterations);
        writetable(Excel_Investment,filename,'Sheet',sheet);        
    else
        disp("The GEP case is not feasible.")
        return   % In case the mathematical formulation is not feasible, we stop the code.
    end
end

function Renew_commitment = calc_renew_commitment(opt_year)
    % Renewable quota
    % We assume that by 2020 the quota is already 50%, and by 2050 it must be
    % at least 90%, so we will do a linear function to represent this:
    % y = m*x + b
    final = 0.9; % final value commitment at the end of the study horizon
    m = (final - 0.5)/(31-1); % slope of the line
    b = 0.5 - m;              % y-intercept
    year = str2double(opt_year) - 2019;    % normalize so 2019 will be year 0
    
    Renew_commitment = m*year+ b;       % output
end

function test_days = ident_days_to_test(Global_Chronological_demand,opt_year,renewables_gen,Investment,basedemand_length,DBB_info)  
    global qt_days
    % Find the day with the least net load: Demand - Renewable generation
    % Chronological demand of the optimized year
    Load = Global_Chronological_demand.(opt_year);

    % Retrieve renewable generation data
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    Flags = DBB_info.Flags;
    
    %% Hourly Energy Generation
    renew_gen = zeros(height(Load),1); % vector that will be used to add all of the generation of the renewables
    % Renewables generation - existing:
    for i = 1:height(Renewables_exist_data)
        renew_gen = renew_gen + renewables_gen.(Renewables_exist_data.RenewableScenario{i})*Renewables_exist_data.UpperLimit(i);
    end
    % Renewable generation - candidate:
    for i = 1:height(Renewables_cand_data)
        renew_gen = renew_gen + renewables_gen.(Renewables_cand_data.RenewableScenario{i})*Renewables_cand_data.UpperLimit(i)*Investment(Flags.flag_EcoThermal_cand+i);
    end
    
    % NetLoad
    NetLoad = Load - renew_gen;
    
    % Table for ordered information: timestamps - NetLoad
    % Use blocks flag in ordered demand base
    year = str2double(opt_year);
    % First make a first date
    t = datetime([year,1,1]);
    % Second, convert it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(basedemand_length,1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','day','current');
    
    % Concatenate both time_stamps and NetLoad
    NetLoad_table = table(time_stamps,NetLoad);
    
    % First, we sort NetLoad_table, in ascending order, so at the begginning we will have the least NetLoad Look for the minimum 5 NetLoads
    NetLoad_table = sortrows(NetLoad_table,2);
    % Take out the first column, while at the same time apply a unique
    sort_time = unique(NetLoad_table.time_stamps,'stable');
    % Finally, take out the first qt_days unique instances:
    least_NetLoad = sort_time(1:qt_days);
    
    % Second, we sort NetLoad_table, in descending order, so at the begginng we will have the biggest NetLoad, so instances in which Load is big and renewable generation is low
    NetLoad_table = sortrows(NetLoad_table,2,'descend');
    % Take out the first column, while at the same time apply a unique, in case that we have repeating days at the begginning
    sort_time = unique(NetLoad_table.time_stamps,'stable');
    % Finally, we extract the data that we want:
    most_NetLoad = sort_time(1:qt_days);

    % Third, we will look for the instances with the most and least renewable energy production
    Renew_gen_table = table(time_stamps,renew_gen);
    % Sort the table based on renew_gen, which is the second variable:
    Renew_gen_table = sortrows(Renew_gen_table,2);
    % Take out the first column, while at the same time apply a unique:
    sort_time = unique(Renew_gen_table.time_stamps,'stable');
    % Extract the data we want:
    least_Renew_gen = sort_time(1:qt_days);
    
    % At the same time, we just reshuffle the table to order it in descending order:
    Renew_gen_table = sortrows(Renew_gen_table,2,'descend');
    % Take out the first column, while applying a unique:
    sort_time = unique(Renew_gen_table.time_stamps,'stable');
    % Extract the top qt_days data:
    most_Renew_gen = sort_time(1:qt_days);

    % Fourth, we decided we want to use Most Load. We make this change based on the logic of the four previous 
    Load_table = table(time_stamps,Load);   % So we first convert Load into a table partnered with their respective time stamps
    % Now we sort the Load table, according to Load, so according to our second variable
    Load_table = sortrows(Load_table,2,'descend');   % the descend attribute means that looking down the vector, the values descend in value. So the biggest values sit atop.
    % Now we take out the dates we desire from the first column
    sort_time = unique(Load_table.time_stamps,'stable');
    % We only extract the amount of unique days we desire according to the quantity of days we will be testing
    most_Demand = sort_time(1:qt_days);

    % Fifth, we want a set of random days.
    Random = datetime([repmat(year,qt_days,1), round((12-1).*rand(qt_days,1)) + 1,round((12-1).*rand(qt_days,1)) + 1]);


    %% Order all of the output in a single Table
    test_days = table(least_NetLoad,most_NetLoad,least_Renew_gen,most_Renew_gen,most_Demand,Random);
end

function [Week_Chrono_LNL,Week_Chrono_MNL,Week_Chrono_LRG,Week_Chrono_MRG,Week_Chrono_MDM,Week_Chrono_RAN] = Load_Chrono_to_test(opt_year,Global_Chronological_demand,test_days)
    % Extract Chronological Demand for optimized year
    Load = Global_Chronological_demand.(opt_year);
    % Table with time and load
    % First make a first date
    t = datetime([str2double(opt_year),1,1]);
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(length(Load),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    
    % Finally, we convert back to date format, while at the same time round
    % to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');
    
    % 1. Create a matrix for the Loads of all tested days:
    % demand for Least Net Load
    Load_Chrono_LNL = zeros(24,length(test_days.least_NetLoad)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.least_NetLoad)
        idx = find(time_stamps == string(test_days.least_NetLoad(i)));        % the command find verifies in which position does the condition happen. In this case the condition is where does the date matches in the vector of dates
        Load_Chrono_LNL(:,i) = Load(idx:idx+23);                              % we then use that flag and extract the Load for the same position + the following 23 cases.
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_LNL = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_LNL = [Week_Chrono_LNL;Load_Chrono_LNL(:,i)];
    end

    % 1. demand for Most Net Load
    Load_Chrono_MNL = zeros(24,length(test_days.most_NetLoad)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.most_NetLoad)
        idx = find(time_stamps == string(test_days.most_NetLoad(i)));
        Load_Chrono_MNL(:,i) = Load(idx:idx+23);
    end
     % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MNL = [];
    for i = 1:width(Load_Chrono_MNL)
        Week_Chrono_MNL = [Week_Chrono_MNL;Load_Chrono_MNL(:,i)];
    end

    % 1. demand for least Renewable Generation
    Load_Chrono_LRG = zeros(24,length(test_days.least_Renew_gen)); 
    % 2. Find the position for the respective days:
    for i = 1:length(test_days.least_Renew_gen)
        idx = find(time_stamps == string(test_days.least_Renew_gen(i)));
        Load_Chrono_LRG(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_LRG = [];
    for i = 1:width(Load_Chrono_LRG)
        Week_Chrono_LRG = [Week_Chrono_LRG;Load_Chrono_LRG(:,i)];
    end

    % 1. demand for Most Renewable Generation
    Load_Chrono_MRG = zeros(24,length(test_days.most_Renew_gen)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.most_Renew_gen)
        idx = find(time_stamps == string(test_days.most_Renew_gen(i)));
        Load_Chrono_MRG(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MRG = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_MRG = [Week_Chrono_MRG;Load_Chrono_MRG(:,i)];
    end

    % 1. demand for Most Demand
    Load_Chrono_MDM = zeros(24,length(test_days.most_Demand)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.most_Demand)
        idx = find(time_stamps == string(test_days.most_Demand(i)));
        Load_Chrono_MDM(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_MDM = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_MDM = [Week_Chrono_MDM;Load_Chrono_MDM(:,i)];
    end

    % 1. demand for Random
    Load_Chrono_RAN = zeros(24,length(test_days.Random)); 
    % 2. Find the position for the respective days:    
    for i = 1:length(test_days.Random)
        idx = find(time_stamps == string(test_days.Random(i)));
        Load_Chrono_RAN(:,i) = Load(idx:idx+23);
    end
    % 3. Interlink all of the days in a single vector to be used as Load
    Week_Chrono_RAN = [];
    for i = 1:width(Load_Chrono_LNL)
        Week_Chrono_RAN = [Week_Chrono_RAN;Load_Chrono_RAN(:,i)];
    end
end

function [Renewables_chrono_exist_scenario_concat,Renewables_chrono_cand_scenario_concat,Fuel_exist_index,Fuel_cand_index] = UCP_renew_gen_and_fuel(renewables_gen,dates,opt_year,DBB_info)
    % Call information to know the 
    %% Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    Fuel_price_data = DBB_info.Fuel_price_data;
    Flags = DBB_info.Flags;

    %%% RENEWABLES
    %% 1. Extract the information for the generators
    % First make a first date
    t = datetime([str2double(opt_year),1,1]);
    % Second, converti it to a number
    t = datenum(t);
    % Third, assign to time_stamps(1) as initial date
    time_stamps = NaN(height(renewables_gen),1);
    time_stamps(1) = t;
    % fourth, fill up the correct numbers for all the long vector time
    for ii = 2:length(time_stamps)
        time_stamps(ii) = time_stamps(ii-1) + 1/24; % for the length in time_stamps we add one hour
    end
    % Finally, we convert back to date format, while at the same time round
    % to the nearest hour
    time_stamps = dateshift(datetime(time_stamps,'ConvertFrom','datenum'),'start','hour','nearest');

    % We now delete the variable time stamps from renewables_gen, making a
    % new variable with a slightly different name
    Renewables_gen  = removevars(renewables_gen,{'time_stamp'});
    % Immediately after, insert the new 'time_stamp' column:
    Renewables_gen  = addvars(Renewables_gen,time_stamps,'Before','block_flag');
    
    % Extract the ids of the dates that we need to extract
    id   = zeros(length(dates),1);
    for i = 1:length(dates)                                     % Warning: no need to fix when stumbling into trouble.          
        id(i) = find(Renewables_gen.time_stamps == dates(i));   % Warning: it stumbles into trouble when used outside the main code because the test daays change between loop,
    end                                                         % but when tested, the variable test_days is not updated while Renewables_gen is.

    %% 2. Extract the Chronological Generation for the Existing generators
    % Solar Chronological resource
    Renewables_chrono_exist_scenario = nan(24,height(Renewables_exist_data),length(dates));  % Make a NAN 3-D array to be replaced with information extracted from the chronological renewable scenario
    for j = 1:length(id) 
        counter = 1;        % Given that we need to start the counter at the top line of the RenewablesExist.xlsb file
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist           % for depending on the amount of existing solar generators
            Renewables_exist = Renewables_gen.(Renewables_exist_data.RenewableScenario{counter});     % We extract the information for the specific renewable scenario (all year)
            Renewables_exist = Renewables_exist(id(j):id(j)+23);                                      % We cut the information just for the dates we are interested in
            Renewables_chrono_exist_scenario(:,counter,j) = Renewables_exist;                         % We add to the master array
            counter = counter + 1;
        end
    end
    % Make the 3-D array into a 2-D array by concatenating the matrices into a single matrix.
    Renewables_chrono_exist_scenario_concat = []; % Empty matrix variable that in which the other information will be pasted to.
    for i = 1:size(Renewables_chrono_exist_scenario,3)
        Renewables_chrono_exist_scenario_concat = [Renewables_chrono_exist_scenario_concat;Renewables_chrono_exist_scenario(:,:,i)];  %% We now just append below the information of all existing scenarios in a single 2-D Matrix
    end

    %% 3. Extract the Chronological Generation for the Candidate generators
    % Solar Chronological resource
    Renewables_chrono_cand_scenario = nan(24,height(Renewables_cand_data),length(dates));  % Make a NAN 3-D array to be replaced with information extracted from the chronological renewable scenario
    for j = 1:length(id) 
        counter = 1;        % Given that we need to start the counter at the top line of the SolarExit.xlsb file
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand          % for depending on the amount of existing solar generators
            Renewables_cand = Renewables_gen.(Renewables_cand_data.RenewableScenario{counter});     % We extract the information for the specific renewable scenario (all year)
            Renewables_cand = Renewables_cand(id(j):id(j)+23);                                      % We cut the information just for the dates we are interested in
            Renewables_chrono_cand_scenario(:,counter,j) = Renewables_cand;                         % We add to the master array
            counter = counter + 1;
        end
    end
    % Make the 3-D array into a 2-D array by concatenating the matrices into a single matrix.
    Renewables_chrono_cand_scenario_concat = []; % Empty matrix variable that in which the other information will be pasted to.
    for i = 1:size(Renewables_chrono_exist_scenario,3)
        Renewables_chrono_cand_scenario_concat = [Renewables_chrono_cand_scenario_concat;Renewables_chrono_cand_scenario(:,:,i)];  %% We now just append below the information of all existing scenarios in a single 2-D Matrix
    end

    %%% FUEL PRICES
    % Existing Fuel Forecast: Thermal
    Fuel_exist_index = [];
    for i = 1:Flags.flag_thermal_exist             % First: finds the fuel of technology i. Then, extracts the Fuel_price_data for that fuel.
        Fuel_exist_index = [Fuel_exist_index,Fuel_price_data.(string(Thermal_exist_data.Fuel{i}))(opt_year)]; % Finally, it obtains just a number indexed for year opt_year
    end
    % Existing Fuel Variation: EcoThermal
    for i = 1: Flags.flag_EcoThermal_exist - Flags.flag_thermal_exist
        Fuel_exist_index = [Fuel_exist_index,Fuel_price_data.(string(EcoThermal_exist_data.Fuel{i}))(opt_year)];
    end
    
    % Candidate Fuel Variation: Thermal
    Fuel_cand_index = [];
    for i = 1:Flags.flag_thermal_cand
        Fuel_cand_index = [Fuel_cand_index,Fuel_price_data.(string(Thermal_cand_data.Fuel{i}))(opt_year)];
    end
    % Candidate Fuel Variatin: EcoThermal
    for i = 1: Flags.flag_EcoThermal_cand - Flags.flag_thermal_cand
        Fuel_cand_index = [Fuel_cand_index,Fuel_price_data.(string(EcoThermal_cand_data.Fuel{i}))(opt_year)];
    end
end

function UCP_optimizer = UCP_optimizer_creator(DBB_info)
    global Curt_penalty qt_days
    size_Load = 24*qt_days;
    % Nomenclature:
    % i: Generators. Will be used both for aggregated generators in Existing or Candidate decision variable object.
    % h: hours.
    % t: years. -> However, there is not instance in here where we iterate over years. We only use the optimized year variable to keep track of the year in which we are right now.
    
    % Call BDD information/data
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;

    % General Information
    % Existing Generator Bounds
    Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
    Exist_Lower = [Thermal_exist_data.LowerLimit; EcoThermal_exist_data.LowerLimit; zeros(height(Renewables_exist_data),1)]; % Lower Generation Bounds Existing
    
    % Candidate Generator Bounds
    Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
    Cand_Lower = [Thermal_cand_data.LowerLimit; EcoThermal_cand_data.LowerLimit; zeros(height(Renewables_cand_data),1)]; % Lower Generation Bounds Candidate
      
    %% Mathematical Formulation
    % Decision variables
    % Operation
    Exist_git       = sdpvar(Flags.quant_exist,size_Load,'full');            % continuous variable to define generation output from exist generation
    Cand_git        = sdpvar(Flags.flag_renewables_cand,size_Load,'full');   % continuous variable to define generation output from future generation, excluding ESS
    % Commitment
    CommitExist_zit = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to determine if a Existing unit is committed or not at a given moment
    CommitCand_zit  = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to determine if a Candidate unit is committed or not at a given moment, excluding ESS
    % Startup- and shutdown
    StartExist_yit  = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to flag if an Existing Thermal Technology starts-up
    StartCand_yit   = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to flag if a Candidate Thermal Technology starts-up, excluding ESS
    ShutD_exist_xit = binvar(Flags.quant_exist,size_Load,'full');            % binary variable to flag if an Existing Thermal Technology Shuts-down
    ShutD_cand_xit  = intvar(Flags.flag_renewables_cand,size_Load,'full');   % integer variable to flag if a Candidate Thermal Technology Shuts-down, excluding ESS
    % Battery operation
    q_SOC_ih        = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable decision of amount of energy stored in each type of battery
    q_Discharge_ih  = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable of Discharge by ESS i at hour h
    q_Charge_ih     = sdpvar(height(ESS_cand_data),size_Load,'full');        % continuous variable of Charge by ESS i at hour h
    d_ih            = binvar(height(ESS_cand_data),size_Load,'full');        % binary variable to determine if a battery should be charging or discharging itself
    % Curtailment decision variables
    Curt_Exist_git  = sdpvar(Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,size_Load,'full');    % continuous variable to keep track of curtailment of existing renewables
    Curt_Cand_git  = sdpvar(Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,size_Load,'full');       % continuous variable to keep track of curtailment of existing renewables

    % Parameters
    % Generators
    Investment = sdpvar(Flags.quant_cand,1,'full'); % Parameter for the vector of Investment passed by the GEP
    % Load
    Load       = sdpvar(size_Load,1,'full');        % Parameter for load changing per UCP tested. This changes by year, and by iteration of the GEP
    % Renewables
    Renewables_chrono_exist_scenario = sdpvar(size_Load,Flags.flag_renewables_exist-Flags.flag_EcoThermal_exist,'full');    % Parameter that will dictated the resource available for that specific hour for existing renewables
    Renewables_chrono_cand_scenario  = sdpvar(size_Load,Flags.flag_renewables_cand-Flags.flag_EcoThermal_cand,'full');      % Parameter that will dictated the resource available for that specific hour for candidate renewables
    % Fuel Price
    Fuel_exist_index = sdpvar(1,Flags.flag_EcoThermal_exist,'full');     % Parameter that will be multiplied against the energy generated by existing thermal power, including ecothermal.
    Fuel_cand_index  = sdpvar(1,Flags.flag_EcoThermal_cand,'full');      % Parameter that will be multiplied against the energy generated by candidate thermal power, including ecothermal.
    
    %% Constraints
    % General Constraints
    Constraints = [];
    % Power Balance Constraint
    for h = 1:size_Load
        Constraints = Constraints + [(sum(Exist_git(:,h)) + sum(Cand_git(:,h)) + sum(q_Discharge_ih(:,h)) == Load(h) + sum(q_Charge_ih(:,h))) : 'Power Balance Constraint']; % 
    end

    % Generator Bounds: Min and Max Constraints
    % Existing Generator Bounds
    for h = 1:size_Load          % This constraint reads: [Binvar*LowerBound <= Generation <= Binvar*UpperBound ]; In other words, if a technology is committed its generation must be within bounds. If not, it must be zero.
        Constraints = Constraints + [(CommitExist_zit(:,h).*Exist_Lower <= Exist_git(:,h) <= CommitExist_zit(:,h).*Exist_Upper): 'Existing Generator Bounds'];
    end
    % Candidate Generator Bounds
    for h = 1:size_Load          % This constraint reads: [Binvar*LowerBound <= Generation <= Binvar*UpperBound ]; In other words, if a technology is committed its generation must be within bounds. If not, it must be zero.
        Constraints = Constraints + [(CommitCand_zit(:,h).*Cand_Lower <= Cand_git(:,h) <= CommitCand_zit(:,h).*Cand_Upper): 'Candidate Generator Bounds']; 
    end

    % Commitment Integer Upper bound
    for h = 1:size_Load          % This constraints reads: zit <= Ni, wher zit is the commitment and Ni is the installed number of technology i
        Constraints = Constraints + [(CommitCand_zit(:,h) <= Investment(1:Flags.flag_renewables_cand)) : 'Commitment at most to Installed Candidate'];
    end


    %% Thermal Constraints
    % Ramp up and Ramp down Constraints
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        for h = 2:size_Load
            Constraints = Constraints + [(Exist_git(i,h-1) - Thermal_exist_data.RampDown(i)*CommitExist_zit(i,h) <= Exist_git(i,h) <= Exist_git(i,h-1) + Thermal_exist_data.RampUp(i)*CommitExist_zit(i,h)): 'Existing Thermal Ramp up and Down Constraint'];   % Ramp-up and ramp-down constraint
        end
    end
    % Existing EcoThermal Generators
    counter = 1;        % counter to distinguish when counting for the Ecothermal data and the Git matrix.
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        for h = 2:size_Load
            Constraints = Constraints + [(Exist_git(i,h-1) - EcoThermal_exist_data.RampDown(counter)*CommitExist_zit(i,h) <= Exist_git(i,h) <= Exist_git(i,h-1) + EcoThermal_exist_data.RampUp(counter)*CommitExist_zit(i,h)): 'Existing EcoThermal Ramp up and down Constraint'];   % Ramp-up and ramp-down constraint
        end
        counter = counter + 1;
    end
    % Cand Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        for h = 2:size_Load          % Ramp-up and ramp-down constraint
            Constraints = Constraints + [(Cand_git(i,h-1) - Thermal_cand_data.RampDown(i)*CommitCand_zit(i,h) <= Cand_git(i,h) <= Cand_git(i,h-1) + Thermal_cand_data.RampUp(i)*CommitCand_zit(i,h)) : 'Candidate Thermal Ramp up and down Constraint'];
        end
    end
    % Candidate EcoThermal Generators
    counter = 1;        % counter to distinguish when counting for the Ecothermal data and the Git matrix.
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        for h = 2:size_Load          % Ramp-up and ramp-down constraint
            Constraints = Constraints + [(Cand_git(i,h-1) - EcoThermal_cand_data.RampDown(counter)*CommitCand_zit(i,h) <= Cand_git(i,h) <= Cand_git(i,h-1) + EcoThermal_cand_data.RampUp(counter)*CommitCand_zit(i,h)) :'Candidate EcoThermal Ramp up and down Constraint'];
        end
        counter = counter + 1;
    end

    % Start Indicator: All Technologies at the same time.          % Observation: We are defining a flag for start for Renewables, but we will defined a Cost Zero, so the solver might take whatever value it decides to be okay.
    for h = 2:size_Load
        Constraints = Constraints + [(CommitExist_zit(:,h) <= CommitExist_zit(:,h-1) + StartExist_yit(:,h)) : 'Startup Indicator Existing Generators'];
        Constraints = Constraints + [(CommitCand_zit(:,h) <=  CommitCand_zit(:,h-1)  + StartCand_yit(:,h)): 'Startup Indicator Candidate Generators'];
    end
    % Shutdown Indicator: All Techonologies done at the same time. % Observation: We are defining a flag for start for Renewables, but we will defined a Cost Zero, so the solver might take whatever value it decides to be okay.                             
    for h = 2:size_Load                                
        Constraints = Constraints + [(CommitExist_zit(:,h) >= CommitExist_zit(:,h-1) - ShutD_exist_xit(:,h)): 'Shutdown Indicator Existing Generators'];
        Constraints = Constraints + [(CommitCand_zit(:,h) >= CommitCand_zit(:,h-1) - ShutD_cand_xit(:,h)): 'Shutdown Indicator Candidate Generators'];
    end    

    % MinUp and MinDown time
    % MinUp and MinDown Existing Thermal Generators
    % Existing Thermal Generators
    for h = 2:size_Load
        for i = 1:Flags.flag_thermal_exist
            rangeUp   = h:min(size_Load,h+Thermal_exist_data.MinUp(i)-1);
            rangeDown = h:min(size_Load,h+Thermal_exist_data.MinDown(i)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitExist_zit(i,rangeUp) >= StartExist_yit(i,h)): 'MinUp Time Existing Thermal'];        % MinUp time
            Constraints = Constraints + [(CommitExist_zit(i,rangeDown) <= 1 - ShutD_exist_xit(i,h)): 'MinDown Time Existing Thermal']; % MinDown time
        end
    end
    % Existing EcoThermal Generators
    for h = 2:size_Load
        counter = 1;        % counter to distinguish EcoThermal and Git.
        for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
            rangeUp =   h:min(size_Load,h+EcoThermal_exist_data.MinUp(counter)-1);
            rangeDown = h:min(size_Load,h+EcoThermal_exist_data.MinDown(counter)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitExist_zit(i,rangeUp) >= StartExist_yit(i,h)): 'MinUp Time Existing EcoThermal'];        % MinUp time
            Constraints = Constraints + [(CommitExist_zit(i,rangeDown) <= 1 - ShutD_exist_xit(i,h)): 'MinDown Time Existing EcoThermal']; % MinDown time
            counter = counter + 1;
        end
    end
    % Candidate Thermal Generators
    for h = 2:size_Load
        for i = 1:Flags.flag_thermal_cand
            rangeUp =   h:min(size_Load,h+Thermal_cand_data.MinUp(i)-1);
            rangeDown = h:min(size_Load,h+Thermal_cand_data.MinDown(i)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitCand_zit(i,rangeUp) >= StartCand_yit(i,h)): 'MinUp Time Candidate Thermal'];         % MinUp time
            Constraints = Constraints + [(CommitCand_zit(i,rangeDown) <= Investment(i) - ShutD_cand_xit(i,h)): 'MinDown Time Candidate EcoThermal'];  % MinDown time
        end
    end
    % Candidate EcoThermal Generators
    for h = 2:size_Load
        counter = 1;        % counter to distinguish EcoThermal and Git.
        for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
            rangeUp =   h:min(size_Load,h+EcoThermal_cand_data.MinUp(counter)-1);
            rangeDown = h:min(size_Load,h+EcoThermal_cand_data.MinDown(counter)-1);
            % Add actual Constraints
            Constraints = Constraints + [(CommitCand_zit(i,rangeUp) >= StartCand_yit(i,h)): 'MinUp Time Candidate EcoThermal'];        % MinUp time
            Constraints = Constraints + [(CommitCand_zit(i,rangeDown) <= Investment(i) - ShutD_cand_xit(i,h)): 'MinDown Time Candidate EcoThermal']; % MinDown time
            counter = counter + 1;
        end
    end
    
    %% Renewables generation resource: This constraint is read that the generation from renewables needs to be at most the available resource
    % Existing renewables
    counter = 1; % fictitious variable to start the count for the length of renewables
    for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
        Constraints = Constraints + [(Exist_git(i,:)' <= Exist_Upper(i).*Renewables_chrono_exist_scenario(:,counter)): 'Existing Renewable Available Resource'];
        counter = counter +1;
    end
    % Candidate renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
        Constraints = Constraints + [(Cand_git(i,:)' <= (Investment(i)*Cand_Upper(i)).*Renewables_chrono_cand_scenario(:,counter)): 'Candidate Renewable Available Resource'];
        counter = counter +1;
    end

    % Constraints to keep track of the curtailment
    % Existing renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
        Constraints = Constraints + [(Curt_Exist_git(counter,:)' == Exist_Upper(i).*Renewables_chrono_exist_scenario(:,counter) - Exist_git(i,:)'): 'Existing Renewable Curtailment'];
        counter = counter +1;       % Given that Curt_Exist_git dimensions are only for renewables, we keep track with 'counter' instead of 'i'
    end
    % Candidate renewables
    counter = 1;
    for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
        Constraints = Constraints + [(Curt_Cand_git(counter,:)' == (Investment(i)*Cand_Upper(i)).*Renewables_chrono_cand_scenario(:,counter) - Cand_git(i,:)'): 'Candidate Renewable Curtailment'];
        counter = counter +1;       % Given that Curt_Exist_git dimensions are only for renewables, we keep track with 'counter' instead of 'i'
    end                

    %% ESS Additional constraints
    % State of charge constraints
    for i = 1:height(ESS_cand_data)
        for h = 2:size_Load
            Constraints = Constraints + [(q_SOC_ih(i,h) == q_SOC_ih(i,h-1) - q_Discharge_ih(i,h) + ESS_cand_data.Efficiency(i)*q_Charge_ih(i,h)): 'SOC energy conservation'];
        end
    end
    % Initial and final state of charge per day
    hours = [0:24:size_Load];
    hours(1)= 1;
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:length(hours)
            Constraints = Constraints + [(q_SOC_ih(i,hours(h)) == 1/2*ESS_cand_data.Duration(i)*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'SOC same at the end of everyday'];
        end
        counter = counter +1;
    end
    % Maximum energy stored at the batteries
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:size_Load
            Constraints = Constraints + [(0 <= q_SOC_ih(i,h) <= ESS_cand_data.Duration(i)*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Energy Bounds in Batteries'];
        end
        counter = counter +1;
    end
    % Charge and Discharge bounds
    counter = Flags.flag_renewables_cand+1;
    for i = 1:height(ESS_cand_data)
        for h = 1:size_Load
            Constraints = Constraints + [(0 <= q_Discharge_ih(i,h) <= d_ih(i,h)*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Batteries Discharge Power Bound'];
            Constraints = Constraints + [(0 <= q_Charge_ih(i,h) <= (1-d_ih(i,h))*ESS_cand_data.UpperLimit(i)*Investment(counter)): 'Batteries Charge Power Bound' ];
        end
        counter = counter +1;
    end

    %% Security Constraints
    % Maximum Curtailment: um(Curt_Exist_git) + sum(Curt_Cand_exist) <= 0.6(Exist_ResourceAvailability*RenewablesInstalledCapacity + Cand_RenewablesInstalledCapacity*Investment*ResourceAvailability)
    % Constraints = Constraints + [sum(sum(Curt_Exist_git)) + sum(sum(Curt_Cand_git)) <= 0.5*(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand))))];

    %% Cap for UCP Payment
    % We will define a maximum payment for UCP within a week. The cap is 1billion USD
    % This Constraint was built after the Objective Function. However, it will be based on the structure of the Objective Function.
    % We define the LHS as the sum of costs, while the RHS is just the cap = 1billion USD
    Cap_LHS = 0;
    % Add Thermal generation cost - Existing
    for i = 1:Flags.flag_thermal_exist
        Cap_LHS = Cap_LHS + sum(Fuel_exist_index(i)*Thermal_exist_data.VarCost(i)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Existing
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cap_LHS = Cap_LHS + sum(Fuel_exist_index(i)*EcoThermal_exist_data.VarCost(counter)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % Add Thermal generation cost - Candidate
    for i = 1:Flags.flag_thermal_cand
        Cap_LHS = Cap_LHS + sum(Fuel_cand_index(i)*Thermal_cand_data.VarCost(i)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Candidate
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cap_LHS = Cap_LHS + sum(Fuel_cand_index(i)*EcoThermal_cand_data.VarCost(counter)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % On the section "Add generation cost" the second term is a vector, so we will insert a sum function to add the whole vector and minimize a single Objective Function
    Cap_LHS = sum(Cap_LHS);
    
    % Add StartUp Cost and Shutdown Cost                % We only sum from the second hour onwards, because for the first hour it might not be defined.
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        Cap_LHS = Cap_LHS + sum(StartExist_yit(i,2:end))*Thermal_exist_data.StartUpCost(i);
        Cap_LHS = Cap_LHS + sum(ShutD_exist_xit(i,2:end))*Thermal_exist_data.StartUpCost(i);
    end
    % Existing EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cap_LHS = Cap_LHS + sum(StartExist_yit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        Cap_LHS = Cap_LHS + sum(ShutD_exist_xit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Candidate Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        Cap_LHS = Cap_LHS + sum(StartCand_yit(i,2:end))*Thermal_cand_data.StartUpCost(i);
        Cap_LHS = Cap_LHS + sum(ShutD_cand_xit(i,2:end))*Thermal_cand_data.StartUpCost(i);
    end
    % Candidate EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cap_LHS = Cap_LHS + sum(StartCand_yit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        Cap_LHS = Cap_LHS + sum(ShutD_cand_xit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Add Curtailment penalty cost: 1. When we first sum(Curtailment) we obtain a vector. 2. When we sum(sum(Curtailment)) we get a number. This number is added to the cost multiplied by 100 USD/MWh.
    % Existing Curtailment
    Cap_LHS = Cap_LHS + 300*sum(sum(Curt_Exist_git));
    % Candidate Curtailment
    Cap_LHS = Cap_LHS + 300*sum(sum(Curt_Cand_git));

    % Finally add the constraint
    Constraints = Constraints + [(Cap_LHS <= 1000000000): 'Cap UCP Financial Resource'];

    %% Objective Function    % FIX: THE UCP IS NOT INDEXING THE COST OF FUEL ACROSS TIME. Refer to v8. UCP_optimizer(sOR,ramp,opt_year). This also affects the 'Cap UCP Financial Resource' contraint.
    % Initialize Cost
    Cost = 0;
    % Add Thermal generation cost - Existing
    for i = 1:Flags.flag_thermal_exist
        Cost = Cost + sum(Fuel_exist_index(i)*Thermal_exist_data.VarCost(i)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Existing
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cost = Cost + sum(Fuel_exist_index(i)*EcoThermal_exist_data.VarCost(counter)*Exist_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % Add Thermal generation cost - Candidate
    for i = 1:Flags.flag_thermal_cand
        Cost = Cost + sum(Fuel_cand_index(i)*Thermal_cand_data.VarCost(i)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
    end
    % Add EcoThermal generation cost - Candidate
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cost = Cost + sum(Fuel_cand_index(i)*EcoThermal_cand_data.VarCost(counter)*Cand_git(i,:)); % This formula reads: sum([VariableCost(USD/MWh)]*[Production(MW)]);
        counter = counter + 1;
    end
    % On the section "Add generation cost" the second term is a vector, so we will insert a sum function to add the whole vector and minimize a single Objective Function
    Cost = sum(Cost);
    
    % Add StartUp Cost and Shutdown Cost                % We only sum from the second hour onwards, because for the first hour it might not be defined.
    % Existing Thermal Generators
    for i = 1:Flags.flag_thermal_exist
        Cost = Cost + sum(StartExist_yit(i,2:end))*Thermal_exist_data.StartUpCost(i);
        Cost = Cost + sum(ShutD_exist_xit(i,2:end))*Thermal_exist_data.StartUpCost(i);
    end
    % Existing EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_exist+1:Flags.flag_EcoThermal_exist
        Cost = Cost + sum(StartExist_yit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        Cost = Cost + sum(ShutD_exist_xit(i,2:end))*EcoThermal_exist_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Candidate Thermal Generators
    for i = 1:Flags.flag_thermal_cand
        Cost = Cost + sum(StartCand_yit(i,2:end))*Thermal_cand_data.StartUpCost(i);
        Cost = Cost + sum(ShutD_cand_xit(i,2:end))*Thermal_cand_data.StartUpCost(i);
    end
    % Candidate EcoThermal Generators
    counter = 1;
    for i = Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand
        Cost = Cost + sum(StartCand_yit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        Cost = Cost + sum(ShutD_cand_xit(i,2:end))*EcoThermal_cand_data.StartUpCost(counter);
        counter = counter +1;
    end
    % Add Curtailment penalty cost: 1. When we first sum(Curtailment) we obtain a vector. 2. When we sum(sum(Curtailment)) we get a number. This number is added to the cost multiplied by 100 USD/MWh.
    % Existing Curtailment
    Cost = Cost + Curt_penalty*sum(sum(Curt_Exist_git));
    % Candidate Curtailment
    Cost = Cost + Curt_penalty*sum(sum(Curt_Cand_git));

    %% Solve the problem
    % Set options for YALMIP
    options = sdpsettings('verbose',1,'solver','gurobi','debug',1,'gurobi.NumericFocus',1,'gurobi.MIPGap',0.01,'gurobi.TimeLimit',1200,'gurobi.Seed',rand()*1000,'gurobi.Heuristics',0.1,'usex0',1);

    % Create optimizer object
    UCP_optimizer = optimizer(Constraints,Cost,options,{Investment,Load,Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_exist_index,Fuel_cand_index},{Exist_git,Cand_git,q_Discharge_ih,q_Charge_ih,Cost,Curt_Exist_git,Curt_Cand_git});
end

function [Output_G_cand,Output_G_exist,errorcode,UCP_optimizer,Case_Curtailment] = UCP_evaluate(renewables_gen,dates,opt_year,Investment,Load,UCP_optimizer,UCP_type,DBB_info)
    import java.text.*;     v = DecimalFormat;                         % To use when printing results for easier readibility.
    
    % First: we calculate the renewables and fuel variation to use with the UCP optimizer
    % For that we will call another function that does that
    [Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_Forecast_exist,Fuel_Forecast_cand] = UCP_renew_gen_and_fuel(renewables_gen,dates,opt_year,DBB_info);
    
    % Second: With all of this information and the Load and investment that are inputs from the original function, we can use the optimizer.
    [sol,errorcode,~,~,UCP_optimizer] = UCP_optimizer({Investment,Load,Renewables_chrono_exist_scenario,Renewables_chrono_cand_scenario,Fuel_Forecast_exist,Fuel_Forecast_cand});

    % Third: we will test the errorcode. In case there is no error we will convert the information to be ready to graph
    if errorcode == 3   % In case the errorcode is 3 it means it runs out of time before reaching the desired Optimality Gap.
        errorcode = 0;  % We don't consider this to be an actual error, so we change it back to zero (0).
    end
    if errorcode == 0
        disp("The cost associated with this UCP (" + UCP_type +") is: " + char(v.format(sol{5})) + " USD")   % 5: Cost 
        Output_G_exist = [sol{1}];                                         % 1: Exist_git
        Output_G_cand = [sol{2};sol{3};-sol{4}];                           % 2: Cand_git; 3: q_Discharge_ih; 4: q_Charge_ih
        % Curtailment output to reach a decision about how to give back information to the GEP:
        Thermal_exist_data = DBB_info.Thermal_exist_data;
        EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
        Renewables_exist_data = DBB_info.Renewables_exist_data;    
        Thermal_cand_data = DBB_info.Thermal_cand_data;
        EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
        Renewables_cand_data = DBB_info.Renewables_cand_data;
        Flags = DBB_info.Flags;
    
        % General Information
        % Existing Generator Bounds
        Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
        
        % Candidate Generator Bounds
        Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
        Case_Curtailment = (sum(sum(sol{6})) + sum(sum(sol{7})))/(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand))));
        if Case_Curtailment >= 0.2
            disp("This case curtailment is: " + string(v.format((sum(sum(sol{6})) + sum(sum(sol{7})))/(sum(Renewables_chrono_exist_scenario*Exist_Upper(Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist)) + sum(Renewables_chrono_cand_scenario*(Investment(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand).*Cand_Upper(Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand)))))))
        end
        
    else
        disp("This case is not feasible.")
        Output_G_cand = nan;
        Output_G_exist = nan;
        Case_Curtailment = 0;
    end
end

function fig = graph_EnergyProduction(Output_G_exist,Output_G_cand,opt_year,graph_name_flag,sOR,eval_criteria,visibility,GEP_Load,DBB_info,ramp_capacity,Scenario)
    import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
    %% Extract information of generators Name
    % Call Input Information about existing and candidate generators
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
    Renewables_exist_data = DBB_info.Renewables_exist_data;    
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    Renewables_cand_data = DBB_info.Renewables_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;

    % Name of Generators: Existing
    Name_exist = [Thermal_exist_data.Name;EcoThermal_exist_data.Name;Renewables_exist_data.Name];
    % Name of Generators: Candidate
    Name_cand  = [Thermal_cand_data.Name;EcoThermal_cand_data.Name;Renewables_cand_data.Name;ESS_cand_data.Name + " discharges"];       % We add to the batteries the caveat that those are injections or discharges.
    % Append names
    Name_vector = [Name_exist;Name_cand];

    % Table: Name vs Technology
    fuel_tech = [Thermal_exist_data.Fuel;EcoThermal_exist_data.Fuel;Renewables_exist_data.Tech;Thermal_cand_data.Fuel;EcoThermal_cand_data.Fuel;Renewables_cand_data.Tech;ESS_cand_data.Tech];

    %% Graph output
    switch graph_name_flag
        case "UCP"
            Fig_Naam = "UCP-"+ opt_year +" Hourly Dispatch - "+ eval_criteria;;
        case "GEP"
            Fig_Naam = "GEP Generation by block-"+opt_year + "-sOR: " + num2str(sOR*100) + "% -ramp: " + string(v.format(ramp_capacity)) + "[MW]";
            % For GEP we will make some markers.
            GEP_Load = GEP_Load/1000;
            x_axis_GEP_marker = [];
            for i = 1:height(GEP_Load)
                x_axis_GEP_marker = [x_axis_GEP_marker;string(v.format(GEP_Load(i)))];
            end
    end

    fig = figure("Name",Fig_Naam,'NumberTitle','off',"HandleVisibility",'on',"Visible",visibility);

    % For vertical axis readibility
    flag = 0;
    if max(sum([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)])) - min(sum(Output_G_cand(Flags.quant_cand+1:end,:))) > 20000
        Output_G_cand = Output_G_cand/1000;
        Output_G_exist = Output_G_exist/1000;
        flag = flag + 1;
    end
    
    hold on;
    switch graph_name_flag
        case "UCP"
            title(Fig_Naam,'FontSize',9);
            subtitle("sOR: " + num2str(sOR*100) + "% -ramp: " + string(v.format(ramp_capacity)) + "[MW]-Scenario: " + Scenario,'FontSize',8);
        case "GEP"
            title(Fig_Naam,'FontSize',10);
    end
    % First we will graph the charging of the batteries           
    b = bar(Output_G_cand(Flags.quant_cand+1:end,:)','stacked',"FaceColor","flat");          % When GEP finished, it already gives the charging information negative, so no need to add a negative sign while making the bar graph.    
    % Color code of the bar graph that wll be below the y-axis
    for ii = 1:size(b,2)
        result = char(fuel_tech(ii));
        switch result
            case 'Coal'
                b(ii).CData = [0 0 0];                  % black
            case 'NaturalGas'
                b(ii).CData = [.5 .5 .5];               % dark gray
            case 'Bunker'
                b(ii).CData = [.25 .25 .25];            % half black half gray
            case 'Solar'
                b(ii).CData = [0.9290 0.6940 0.1250];   % opaque yellow
            case 'Wind'
                b(ii).CData = [0.3010 0.7450 0.9330];   % turquoise
            case 'Nuclear'
                b(ii).CData = [0.6350 0.0780 0.1840];   % violet
            case 'Geothermal'
                b(ii).CData = [0.8500 0.3250 0.0980];   % orange
            case 'Li'
                b(ii).CData = [.8 .8 .8];               % silver
        end
    end
    % Now we create the usual information
    b = bar([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)]','stacked',"FaceColor","flat");
    % Color code of the bar graph that will be above the y-axis
    for ii = 1:size(b,2)
        result = char(fuel_tech(ii));
        switch result
            case 'Coal'
                b(ii).CData = [0 0 0];                  % black
            case 'NaturalGas'
                b(ii).CData = [.5 .5 .5];               % dark gray
            case 'Bunker'
                b(ii).CData = [.25 .25 .25];            % half black half gray
            case 'Solar'
                b(ii).CData = [0.9290 0.6940 0.1250];   % opaque yellow
            case 'Wind'
                b(ii).CData = [0.3010 0.7450 0.9330];   % turquoise
            case 'Nuclear'
                b(ii).CData = [0.6350 0.0780 0.1840];   % violet
            case 'Geothermal'
                b(ii).CData = [0.8500 0.3250 0.0980];   % orange
            case 'Li'
                b(ii).CData = [.8 .8 .8];               % silver
        end
    end
    % Figure readibility
    legend([ESS_cand_data.Name + " charging";Name_vector],"Location","bestoutside","FontSize",4);   % We will try to locate it at "best", but only because "bestoutside" cuts out the Title of the graphs
    % For vertical axis readibility
    ytickformat('%,4.0f');
    if flag == 1
        ylabel("Power [GW]");
    else
        ylabel("Power [MW]");
    end
    switch graph_name_flag
        case "UCP"
            xlabel('Time (hours)');
        case "GEP"
            xlabel('Clusters - Representative Demands [GW]');
            xticklabels(x_axis_GEP_marker)
    end
    
    ylim([min(sum(Output_G_cand(Flags.quant_cand+1:end,:))) max(sum([Output_G_exist;Output_G_cand(1:Flags.quant_cand,:)]))*1.05]);

    hold off;
    
end    

function DBB_info = call_BDD_info(Database_name)
    %% Draw all of the information from the BDD files         
    % Define existing generation data
    % Thermal
    Thermal_exist_data = readtable(Database_name,'Sheet',"ThermalExist",'ReadVariableNames',true,"UseExcel",false);
    % Variable cost According to https://www.ods.org.hn/index.php/informes/costes-marginales/2022-costomarginales/mayo22-costosmarginales:
    flag_thermal_exist = height(Thermal_exist_data);
    % EcoThermal: describing a technology that in theory is thermal, but it is not polluting (Nuclear and Geothermal)
    EcoThermal_exist_data = readtable(Database_name,'Sheet',"EcoThermalExist",'ReadVariableNames',true,"UseExcel",false);
    flag_EcoThermal_exist = flag_thermal_exist + height(EcoThermal_exist_data);
    % Renewables
    Renewables_exist_data = readtable(Database_name,'Sheet',"RenewablesExist",'ReadVariableNames',true,"UseExcel",false);
    Renewables_exist_data = sortrows(Renewables_exist_data,"Tech","ascend");
    % Suggestion: make the function to overwrite the Excel file to have it already sorted.
    flag_renewables_exist =  flag_EcoThermal_exist + height(Renewables_exist_data);
    % If this were python, we could use dictionaries and do something similar as use keys to call up specific generators
    quant_exist  = height(Thermal_exist_data) + height(EcoThermal_exist_data) + height(Renewables_exist_data);
       
    % Define candidate generation data
    % Thermal
    Thermal_cand_data = readtable(Database_name,'Sheet',"ThermalCand",'ReadVariableNames',true,"UseExcel",false);
    flag_thermal_cand = height(Thermal_cand_data);
    % EcoThermal: describing a technology that in theory is thermal, but it is not polluting (Nuclear and Geothermal)
    EcoThermal_cand_data = readtable(Database_name,'Sheet',"EcoThermalCand",'ReadVariableNames',true,"UseExcel",false);
    flag_EcoThermal_cand = flag_thermal_cand + height(EcoThermal_cand_data);
    % Renewables
    Renewables_cand_data = readtable(Database_name,'Sheet',"RenewablesCand",'ReadVariableNames',true,"UseExcel",false);
    Renewables_cand_data = sortrows(Renewables_cand_data,"Tech","ascend");
    flag_renewables_cand = flag_EcoThermal_cand + height(Renewables_cand_data);
    % Energy Storage Systems
    % Batteries
    ESS_cand_data = readtable(Database_name,'Sheet',"ESSCand",'ReadVariableNames',true,"UseExcel",false);
    flag_ESS_cand = flag_renewables_cand + height(ESS_cand_data);
    % If this were python, we could use dictionaries and do something similar as use keys to call up specific generators
    quant_cand  = height(Thermal_cand_data) + height(Renewables_cand_data) + height(EcoThermal_cand_data) + height(ESS_cand_data);
    
    % Define the fuel price variation data
    Fuel_price_data = readtable(Database_name,'Sheet',"FuelPrice",'ReadVariableNames',true,"UseExcel",false);
    RowNaams = ["2021"];
    for i = 2:30
        RowNaams(i) = [convertCharsToStrings(num2str(2021-1+i))];
    end
    Fuel_price_data.Properties.RowNames = RowNaams;

    %% Organize all flags in a single table
    Flags = table(flag_thermal_exist,flag_EcoThermal_exist,flag_renewables_exist,quant_exist,flag_thermal_cand,flag_EcoThermal_cand,flag_renewables_cand,flag_ESS_cand,quant_cand);

    %% Close all Excel files possibly still open
    fclose('all');

    DBB_info = struct('Thermal_exist_data',Thermal_exist_data,'EcoThermal_exist_data',EcoThermal_exist_data,'Renewables_exist_data',Renewables_exist_data,'Thermal_cand_data',Thermal_cand_data,'EcoThermal_cand_data',EcoThermal_cand_data,'Renewables_cand_data',Renewables_cand_data,'ESS_cand_data',ESS_cand_data,'Fuel_price_data',Fuel_price_data,'Flags',Flags);
end

function GEP_Load = calc_GEP_Load(Global_LDC_Demand,opt_LDC)
    global fst_study_horizon last_study_horizon
    % This functions takes in as argument: Global Demand with clusters and the clustered table for the last year of the study horizon
    % It will give out the Load needed for the GEP problem. So it needs to be a matrix or table with dimension:
    % (blocks*stages,last_year - fst_year + 1)


    % Extract the headers of Global_LDC_Demand variable
    header = Global_LDC_Demand.Properties.VariableNames';
    % Store in a vector the header information in double form
    horizon = str2double(header);

    % Create directly a table with the variables that we will use
    GEP_Load = unique(Global_LDC_Demand,'stable');
    
    % Now we will constrain the horizon. In case 
    % Create the matrix storing the load duration and the centroids
    for j = 1:length(horizon)
        if fst_study_horizon > horizon(j) || last_study_horizon < horizon(j)    % Skip years that are below or above the study horizon
            GEP_Load = removevars(GEP_Load,string(header(j)));
        end
    end
    
    % We add the Load Duration
    GEP_Load = addvars(GEP_Load,opt_LDC.Hours,'NewVariableNames',"LD",'Before',1);    % We add the load or cluster duration to the GEP Load
end

function graph_GEP(Output_G_exist,Output_G_cand,iteration,visibility,GEP_Load,DBB_info,sOR,Investment)        
    global fst_study_horizon last_study_horizon
    % Documentation: https://undocumentedmatlab.com/articles/export_fig
    %                https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig?s_tid=mwa_osa_a
    % Read ME:       https://github.com/altmany/export_fig/blob/master/README.md
    % Ghostscript:   https://www.mathworks.com/matlabcentral/answers/252411-how-to-actually-use-export_fig
    
    % Basic syntax: export_fig(filename, [handle], options...)
    Route = pwd + "\\Reports";     % Directory name

    % Update Ramp_Capacity measurement
    for t = 1:last_study_horizon-fst_study_horizon+1
        ramp_capacity(t) = ramp_calc(Investment(:,t),nan,DBB_info);     % In this case the ramp capacity referes to the observed ramp capacity. When fed to the GEP, actually it means the minimum ramp capacity.
    end

    for t = 1:size(Output_G_exist,3)    % We then create all of the figures
        opt_year = string(fst_study_horizon+t-1);
        fig = graph_EnergyProduction(Output_G_exist(:,:,t),Output_G_cand(:,:,t),opt_year,"GEP",sOR(t),"",visibility,GEP_Load.(opt_year),DBB_info,ramp_capacity(t),[]);
        % [imageData, alpha] = export_fig:          we disregard imageData and alpha because there is no need for it.
        [~,~] = export_fig(Route + "\GEP_iteration-"+ num2str(iteration) +"_sOR-"+ num2str(max(sOR)*100)+"_ramp "+ num2str(round(max(ramp_capacity))) +"MW.pdf",[fig],'-append');           % In this step we export the figures to pdf. ,'-depsc'
    end

    %% Close all Excel files possibly still open
    fclose('all');      close all hidden;       close all force;
end

function graph_UCP(Output_G_exist,Output_G_cand,sOR,opt_year,eval_criteria,visibility,DBB_info,Investment,Scenario)          %%FIX: This module should use export_fig(). However, for some reason due to Ghoscript internal errors, the Script is not working. I am using right now an ugly print setup.
    global fst_study_horizon last_study_horizon
    Route = pwd + "\\Reports";     % Directory name
    
    % Update Ramp_Capacity measurement
    for t = 1:last_study_horizon-fst_study_horizon+1
        ramp_capacity(t) = ramp_calc(Investment(:,t),nan,DBB_info);     % In this case the ramp capacity referes to the observed ramp capacity. When fed to the GEP, actually it means the minimum ramp capacity.
    end

    fig = graph_EnergyProduction(Output_G_exist,Output_G_cand,opt_year,"UCP",sOR,eval_criteria,visibility,[],DBB_info,ramp_capacity(t),Scenario);     % We create a figure which displays the UCP for a given week. 
    exportgraphics(fig,Route + "\\UCP-year-"+ opt_year +".pdf","Append",true,'Resolution',900);      % https://www.mathworks.com/help/matlab/ref/exportgraphics.html#namevaluepairarguments

    %% Close all Excel files possibly still open
    fclose('all');      close all hidden;       close all force;
end

function [yearly_Installed_Capacity,legend_description] = graph_Peak_Invest(Investment,GEP_Load,flag_peakLoad_InstalledCapacity,yearly_Installed_Capacity,legend_description,renew_gen_block,DBB_info,sOR)
    global blocks stages fst_study_horizon last_study_horizon scenarios
    % For readibility
    import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
    
    % We start by outright trying to graph
    if flag_peakLoad_InstalledCapacity == 1
        % Also, we need to extract the peak demand for each year
        PLoad = max(GEP_Load{:,2:end}); % We only extract peak demand
        
        % We extract the years being graphed
        years = GEP_Load(1,2:end);                  % we extract the information as table
        years = years.Properties.VariableNames;     % we now have the years stored in a cell char array
        label = years;                              % we use this helper variable called label
        years = [];
        for ii = 1:length(label)
            years = [years; str2double(cell2mat(label(ii)))];
        end
    
        fig = figure("Name","Yearly Peak Demand and Firm Power per iteration",'NumberTitle','off',"HandleVisibility",'on',"Visible",'on');     
        
        % First graph is the bar plot demand
        hold on;
        title("Peak Demand and Firm Power per iteration");          
        axis([-inf inf 0 inf]);                                 % We define the floor of the graph as zero. The rest are automatic.
        if max(max(yearly_Installed_Capacity)) >= 30000 || max(PLoad) >= 30000  % We do this to properly scale the graph in case there are big peak demands.
            yearly_Installed_Capacity = yearly_Installed_Capacity/1000;
            PLoad = PLoad / 1000;
            ylabel('Power [GW]')
        else
            ylabel('Power [MW]')
        end
        bar(years,PLoad);                         % We graph peak demand in bar plot, while we will plot the firm power in lines.
        
        % Rest of the graphs are the line plots for each Firm Power demand
        % Plot the iterations Installed Capacity
        for jj = 1:height(yearly_Installed_Capacity)
            plot(years,yearly_Installed_Capacity(jj,:),'-*')
        end
        
        % Figure readibility
        % Change legend description to fit in less information
        divisor = ceil(length(legend_description)/10);              % We will attempt to only have 10 entries in the legend, although there might be more information
        if divisor >= 2
            new_legend_index = length(legend_description):-divisor:1;
            if ~isempty(find(new_legend_index,1))
                new_legend_index = [new_legend_index,1];
            end
            new_legend_index = sort(new_legend_index,'ascend');
    
            % Correct legend description to only allocate 
            for kk = 1:length(legend_description)
                if isempty(find(new_legend_index == kk))
                    legend_description(kk) = '';
                end
            end
        end
    
        % Add Legend
        legend(["Peak Demand";legend_description],"Location","southeast","FontSize",8);
        xlabel('Years');
        ytickformat('%,4.0f');
        hold off
        
        % Export the fig file
        Route = pwd + "\\Reports";     % Directory name
        % [imageData, alpha] = export_fig:          we disregard imageData and alpha because there is no need for it.
        [~,~] = export_fig(Route + "\Yearly Peak Demand and Firm Power per iteration.pdf",[fig],'-append');           % In this step we export the figures to pdf.
    
    else    % In case the previous is not successful, we will add parameters to the output matrices that are keeping track of the firm power and legend
        % 1. First, we need to call the BDD info to make sure we are using the correct placement of all info
        Thermal_exist_data = DBB_info.Thermal_exist_data;
        EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
        Renewables_exist_data = DBB_info.Renewables_exist_data;    
        Thermal_cand_data = DBB_info.Thermal_cand_data;
        EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
        Renewables_cand_data = DBB_info.Renewables_cand_data;
        ESS_cand_data = DBB_info.ESS_cand_data;
        % 2. We now create a vector with the installed capacity of each of the technologies
        % Existing Generator Bounds
        Exist_Installed = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
        % Candidate Generator Bounds
        Cand_Installed = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit;ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
        % 3. Multiply the Installed capacities by their respective peak demand contribution factors Peak Factor coefficient factor contribution for existing generators
        % Renewables Peak Factor coefficient factor
        Exist_renewables_capacity_factor = ones(stages*blocks,height(Renewables_exist_data));
        for s = 1:scenarios
            for t = 1: last_study_horizon - fst_study_horizon +1
                for i = 1:height(Renewables_exist_data)
                    if t == 1
                        Exist_renewables_capacity_factor(:,i) = [renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_exist_data.RenewableScenario{i})];
                    else
                        Exist_renewables_capacity_factor(:,i) = mean([Exist_renewables_capacity_factor(:,i),renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_exist_data.RenewableScenario{i})],2);
                    end
                end
            end
        end
        % Take out the capacity factor of renewables for peak demand
        Exist_renewables_pPF = zeros(1,height(Renewables_exist_data));
        for i = 1:height(Renewables_exist_data)
            for mm = 1:stages
                if mm == 1
                    Exist_renewables_pPF(1,i) = Exist_renewables_capacity_factor((mm-1)*blocks+1,i);
                else
                    Exist_renewables_pPF(1,i) = mean([Exist_renewables_pPF(1,i),Exist_renewables_capacity_factor((mm-1)*blocks+1,i)]);
                end
            end
        end

        % Take out the capacity factor of renewables for peak demand
        pPC_exist = [Thermal_exist_data.HistoricalAvailability;EcoThermal_exist_data.HistoricalAvailability;Exist_renewables_pPF']; % Peak factor of all technologies in peak demand.
        
        % Peak Factor coefficient factor contribution for candidate generators
        Cand_renewables_capacity_factor = ones(stages*blocks,height(Renewables_cand_data));
        for s = 1:scenarios
            for t = 1: last_study_horizon - fst_study_horizon +1
                for i = 1:height(Renewables_cand_data)
                    if t == 1
                        Cand_renewables_capacity_factor(:,i) = [renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_cand_data.RenewableScenario{i})];
                    else
                        Cand_renewables_capacity_factor(:,i) = mean([Cand_renewables_capacity_factor(:,i),renew_gen_block.("S"+s).("Y"+num2str(fst_study_horizon+t-1)).(Renewables_cand_data.RenewableScenario{i})],2);
                    end
                end
            end
        end
        % Take out the capacity factor of renewables for peak demand
        Cand_renewables_pPF = zeros(1,height(Renewables_cand_data));
        for i = 1:height(Renewables_cand_data)
            for mm = 1:stages
                if mm == 1
                    Cand_renewables_pPF(1,i) = Cand_renewables_capacity_factor((mm-1)*blocks+1,i);
                else
                    Cand_renewables_pPF(1,i) = mean([Cand_renewables_pPF(1,i),Cand_renewables_capacity_factor((mm-1)*blocks+1,i)]);
                end
            end
        end
        % Take out the capacity factor of renewables for peak demand
        pPC_cand = [Thermal_cand_data.HistoricalAvailability;EcoThermal_cand_data.HistoricalAvailability;Cand_renewables_pPF';ones(height(ESS_cand_data),1)]; % Peak factor of all technologies in peak demand.        We define the factor of Batteries as able to deliver all of their installed capacity.
    
        % Multiply Exist_Installed.*pPC_exist and also Can_Installed.*pPC_cand
        FirmPower_Exist = Exist_Installed.*pPC_exist;       % Existing
        FirmPower_Cand  = Cand_Installed.*pPC_cand;         % Candidate
    
        % Multiply vector FirmPower_Cand for each of the years in Investment
        yearly_FirmPower_Cand = FirmPower_Cand'*Investment;     % this results in a vector in which each entry is the firm power for each year
        
        % 4. We will make some tags for the legend that differentiates each iteration with sOR and ramp at year it failed.
        % Update Ramp_Capacity measurement
        for t = 1:last_study_horizon-fst_study_horizon+1
            ramp_capacity(t) = ramp_calc(Investment(:,t),nan,DBB_info);     % In this case the ramp capacity referes to the observed ramp capacity. When fed to the GEP, actually it means the minimum ramp capacity.
        end
        % legend_description = [legend_description; "sOR: " + string(max(sOR)*100) + "% ramp: "+ string(v.format(max(ramp_capacity))) + "MW/h"];
        legend_description = [legend_description; "sOR: " + string(max(sOR)*100)];
        yearly_Installed_Capacity = [yearly_Installed_Capacity; sum(FirmPower_Exist) + yearly_FirmPower_Cand];          % This vector stores the yearly Firm Power per iteration    
    end
end

function ramp = ramp_calc(Investment,iterations,DBB_info)
    % Call in information to be used to calculate information
    Thermal_exist_data = DBB_info.Thermal_exist_data;
    EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;   
    Thermal_cand_data = DBB_info.Thermal_cand_data;
    EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
    ESS_cand_data = DBB_info.ESS_cand_data;
    Flags = DBB_info.Flags;
    
    % Add the ramping capacities for existing thermal and EcoThermal
    exist_ramp = sum(Thermal_exist_data.RampUp) + sum(EcoThermal_exist_data.RampUp);
    
    % Calculate the cand_ramp for year evaluated
    cand_ramp = 0;  % Initialize the ramp at zero
    if iterations ~= 0  % Cand_ramp = Thermal Ramp + EcoThermal Ramp + Battery Installed Capacity. We consider that batteries have instantaneous ramps, so it is only limited by their Installed Capacity  
        cand_ramp = cand_ramp + Thermal_cand_data.RampUp'*Investment(1:Flags.flag_thermal_cand) + EcoThermal_cand_data.RampUp'*Investment(Flags.flag_thermal_cand+1:Flags.flag_EcoThermal_cand) + ESS_cand_data.UpperLimit'*Investment(Flags.flag_renewables_cand+1:Flags.flag_ESS_cand);
    end
    
    % Total ramp
    ramp = exist_ramp + cand_ramp;    
end

function Renew_Cap = renew_cap_calc(Investment,renew_gen_block,iterations,DBB_info,GEP_Load,Renew_Cap,t,constrain_factor)
    import java.text.*;     v = DecimalFormat;  % To use when printing results for easier readibility.
    if iterations == 0
        Renew_Cap = GEP_Load.LD'*table2array(GEP_Load(:,2:end));
    else
        % Take out informaiton from DBB_info
        Thermal_exist_data = DBB_info.Thermal_exist_data;
        EcoThermal_exist_data = DBB_info.EcoThermal_exist_data;
        Renewables_exist_data = DBB_info.Renewables_exist_data;    
        Thermal_cand_data = DBB_info.Thermal_cand_data;
        EcoThermal_cand_data = DBB_info.EcoThermal_cand_data;
        Renewables_cand_data = DBB_info.Renewables_cand_data;
        ESS_cand_data = DBB_info.ESS_cand_data;
        Flags = DBB_info.Flags;
        % Existing Generator Bounds
        Exist_Upper = [Thermal_exist_data.UpperLimit; EcoThermal_exist_data.UpperLimit; Renewables_exist_data.UpperLimit]; % Upper Generation Bounds Existing
        % Candidate Generator Bounds
        Cand_Upper = [Thermal_cand_data.UpperLimit; EcoThermal_cand_data.UpperLimit; Renewables_cand_data.UpperLimit; ESS_cand_data.UpperLimit]; % Upper Generation Bounds Candidate
        
        % Measure the desired upper cap
        Measure_Cap = 0;
        % Existing Renewables
        counter = 1;
        for i = Flags.flag_EcoThermal_exist+1:Flags.flag_renewables_exist
            Measure_Cap = Measure_Cap + sum(Exist_Upper(i).*renew_gen_block.(string(Renewables_exist_data.RenewableScenario{counter})).*GEP_Load.LD);
            counter = counter +1;
        end
        % Candidate Renewables
        counter = 1;
        for i = Flags.flag_EcoThermal_cand+1:Flags.flag_renewables_cand
            Measure_Cap = Measure_Cap + sum(Investment(i)*Cand_Upper(i).*renew_gen_block.(string(Renewables_cand_data.RenewableScenario{counter})).*GEP_Load.LD);
            counter = counter + 1;
        end
        % Measurement of Renew Cap
        disp("The amount of renewable energy has been constrained from " + string(v.format(Renew_Cap(t))) + " to " + string(v.format(constrain_factor*min(Measure_Cap,Renew_Cap(t)))))
        Renew_Cap(t) = constrain_factor*min(Measure_Cap,Renew_Cap(t));
    end
end

%% Content