clc
clear
tic
%% Initialization parameters and variable used to save the results

mc_iter = 1000;        % Monte Carlo iterations for statistical analysis (put more than 1 if you want to carry out some statistical analysis via Monte Carlo simulation)
time_steps = 1;  % number of considered time steps (in this code, 1 day with 1 minute resolution)
num_nodes = 99;     % number of nodes of the grid used below (TO BE CHANGED if a different grid is used)

Vtime_true = zeros(time_steps,num_nodes,mc_iter);   %c vector where the "true" voltage results are saved
Vtime_est = zeros(time_steps,num_nodes,mc_iter);    %c vector where the "actual" (coming from the SE + voltage control applications) voltage results are saved


%% Start computations over the decided time steps

for t = 1:1    %%% --> change this if you want just to focus on a specific time step
    %% Part associated to network data;
    %%% --> CHANGE the part in this subsection with the data of the grid to be used
    % ATL-Industrial 100 nodes
    grid_data = Data_reader_industrial();
    Base.V=15e3/sqrt(3); %In V
    Base.S=1e6; % In VA
    time.year=2030;
    time.month=6;
    time.week=7;
    time.day=47;

    [branch, node, ~] = Network_industrial_Atlantide(Base, time, grid_data);
    toc
    tic
    branch.B = zeros(1,branch.num);
   

    [Y, Bs] = Ymatrix(branch, node);                          % computation admittance matrix of the network

    % Power flow calculation for calculating the state assumed as reference and as true operating condition of the network;
    [V, I, S, ~] = Power_flow_NV(branch, node, Y, Bs);

    %% Start Monte Carlo procedure in case of statistical analysis
    branch2=branch;
    node2=node;
    %% Meter placement
    V.mag.meas.index = [1:node.num]';               % indication of the nodes where voltage magnitude meas are placed.
    S.inj.meas.index = [2:node.num]';          % indication of the nodes for the power injection meas
    S.br1.meas.index = [1:1:branch.num]';                % indication of the branches where active and reactive power meas are placed (measurement on the starting side of the branch).
    S.br2.meas.index = []';                % indication of the branches where active and reactive power meas are placed (measurement on the ending side of the branch).
    I.br1.mag.meas.index = []';            % indication of the branches where current magnitude meas are placed (there are no shunt admittances, thus the side of the branch is not important).
    I.br1.mag.measPmu.index = []';                                             % indication of the branches where PMU current meas are placed.
    I.br1.phase.measPmu.index = I.br1.mag.measPmu.index;                       % indication of the branches where PMU current meas are placed.
    V.mag.measPmu.index = []';                                                 % indication of the nodes where PMU voltage meas are placed.
    V.phase.measPmu.index = V.mag.measPmu.index;                               % indication of the nodes where PMU voltage meas are placed.

    %% Uncertainty of the measurements
    Vindex= [1,32,82];                      %Vector containing the indexes of existing voltage meters
    Sbrindex= [];                           %Vector containing the indexes of existing powerflow meters
    V.mag.meas.unc = 3*ones(99,1);           % uncertainty (in percent) of voltage magnitude measurements
    V.mag.meas.unc(Vindex) = 1;              % change the uncertainty of existing meters to 1%
    S.inj.meas.unc = 50;           % uncertainty (in percent) of active and reactive power injection measurements (pseudo-measurements)
    S.br1.meas.unc = 50*ones(98,1);           % uncertainty (in percent) of active and reactive branch power measurements
    S.br1.meas.unc(Sbrindex) = 2;             % change the uncertainty of existing meters to 2%
    S.br2.meas.unc = 0;           % uncertainty (in percent) of active and reactive branch power measurements
    I.br1.mag.meas.unc = 0;       % uncertainty (in percent) of current magnitude measurements
    I.br1.mag.measPmu.unc = 0;                                               % uncertainty (in percent) of current magnitude measurements in PMUs
    I.br1.phase.measPmu.unc = 0;                                             % uncertainty (in centiradians) of current phase-angle measurements in PMUs
    V.mag.measPmu.unc = 0;                                                   % uncertainty (in percent) of voltage magnitude measurements in PMUs
    V.phase.measPmu.unc = 0;                                                 % uncertainty (in centiradians) of voltage phase-angle measurements in PMUs

    %% Creation of the measurements
    [zdatatrue] = Meas_data_true(branch, S, V, I);
        %% Sample creation
    [~,~,zdata_samples]=sampling_call(branch,node,zdatatrue,mc_iter);    %Sample the measurements/branch parameters in their respective PDFs

    %% Uncertainty of the measurements
    %         V.mag.meas.unc = [1];           % uncertainty (in percent) of voltage magnitude measurements
    %         S.inj.meas.unc = 1;           % uncertainty (in percent) of active and reactive power injection measurements (pseudo-measurements)
    %         S.br1.meas.unc = [1];           % uncertainty (in percent) of active and reactive branch power measurements
    %         S.br2.meas.unc = 0;           % uncertainty (in percent) of active and reactive branch power measurements
    %         I.br1.mag.meas.unc = 0;       % uncertainty (in percent) of current magnitude measurements
    %         I.br1.mag.measPmu.unc = 0;                                               % uncertainty (in percent) of current magnitude measurements in PMUs
    %         I.br1.phase.measPmu.unc = 0;                                             % uncertainty (in centiradians) of current phase-angle measurements in PMUs
    %         V.mag.measPmu.unc = 0;                                                   % uncertainty (in percent) of voltage magnitude measurements in PMUs
    %         V.phase.measPmu.unc = 0;                                                 % uncertainty (in centiradians) of voltage phase-angle measurements in PMUs

    %% Creation of the measurements
        zdata=zdatatrue;

    for MC=1:mc_iter
        if rem(MC,100) == 0
            MC
        end
       
        %% Update of branch parameters
%                 branch2.R = branch_samples.R(:,MC); %pick respective set of sampled R
%                 branch2.X = branch_samples.X(:,MC); %pick respective set of sampled X
%                 branch2.B = branch_samples.B(MC,:); %pick respective set of sampled B
        [Y2,Bs2] = Ymatrix(branch2, node2);  % computation admittance matrix of the network with corrupted network parameters


        %% update of measured valued
        %%% --> REPLACE HERE IN THE SECOND COLUMN OF ZDATA THE SAMPLED VALUES
        %%% OF VOLTAGE MAGNITUDES, POWER INJECTIONS AND BRANCH POWERS;
        %%% NB: in the first column: 1=voltage meas; 2=active pwr inj;
        %%% 3=reactive pwr inj; 4=active branch pwr; 5=reactive branch pwr
        zdata(:,2)=zdata_samples(:,MC); %pick respective set of sampled measurements
        

        %% State Estimation

        [V, I, S, num_iter, uncVm_percent] = NV_DSSE_call(branch2, node2, zdata, Y2, Bs2, V, I, S);    % function where the SE algorithm is called

        %% Saving of the results
        Vtime_true(t,:,MC) = V.mag.true_val'; 
        Vtime_est(t,:,MC) = V.mag.est'; 
        Vreal_est(t,:,MC) = V.real.est';
        Vimag_est(t,:,MC) = V.imag.est';
        Vreal_true(t,:,MC) = V.real.true_val';
        Vimag_true(t,:,MC) = V.imag.true_val';


    end
    toc
    XX = [zdata_samples(1:num_nodes,:)' zdata_samples(3*num_nodes-1:end,:)']; %Create X matrix for PCE
    Ypce=squeeze(Vtime_est)'; %Create Y matrix for PCE
    %% UQLAB
    tic
    uqlab %Activate uqlab
    input_params=[zdata([1:num_nodes 3*num_nodes-1:end],2) zdata([1:num_nodes 3*num_nodes-1:end],6)]; %Select mean and std of each of the inputs considered in the PCE
    for ii=1:size(zdata,1)-2*(num_nodes-1)
        IOpts.Marginals(ii).Type = 'Gaussian' ; %select type of PDF for each input
        IOpts.Marginals(ii).Parameters = [input_params(ii,1), input_params(ii,2)] ; %Assign mean and std parameters for each of the inputs
    end
    myInput = uq_createInput(IOpts); %Create input object for uqlab

    %Select PCE as the metamodeling tool
    MetaOpts.Type = 'Metamodel';
    MetaOpts.MetaType = 'PCE';
    %     MetaOpts2.Method = 'LARS';
    % Use input parameters and build PCE model
    MetaOpts.ExpDesign.X = XX;
    MetaOpts.ExpDesign.Y = Ypce;
    MetaOpts.Degree = 2;
    MetaOpts.TruncOptions.qNorm = 0.5;
    MetaOpts.TruncOptions.MaxInteraction = 2;
    myPCE = uq_createModel(MetaOpts); %Run PCE routine
    uq_print(myPCE);
    toc
    tic
    % Sobol indices
    SobolOpts.Type = 'Sensitivity'; %Select sensitivity indices type
    SobolOpts.Method = 'Sobol'; %Select sensitivity indices type
    SobolOpts.Sobol.Order = 1; %Select order of sobol inideces to calculate
    mySobolAnalysisPCE = uq_createAnalysis(SobolOpts); %Compute Sobol indices
    first_order=mySobolAnalysisPCE.Results.FirstOrder; %Save first order indices
    total_order=mySobolAnalysisPCE.Results.Total; %Save total order indices
    toc
    tic


    %% Heatmaps creation
    figure(),h1=heatmap(first_order); %First order heatmap
    xlabel('Voltage magnitude error')
    ylabel('Measured value')
    for ii=1:num_nodes
        h1.XDisplayLabels{ii,1}=['V',num2str(ii)];
    end
    for ii=1:size(zdata,1)
        switch zdata(ii,1)
            case 1
                h1.YDisplayLabels{ii,1}=['V',num2str(zdata(ii,4))];
            case 2
                %                 h1.YDisplayLabels{ii,1}=['P',num2str(zdata(ii,4))];
            case 3
                %                 h1.YDisplayLabels{ii,1}=['Q',num2str(zdata(ii,4))];
            case 4
                h1.YDisplayLabels{ii-2*(num_nodes-1),1}=['P',num2str(zdata(ii,4)),'-',num2str(zdata(ii,5))];
            case 5
                h1.YDisplayLabels{ii-2*(num_nodes-1),1}=['Q',num2str(zdata(ii,4)),'-',num2str(zdata(ii,5))];
        end
    end
    h1.Colormap=parula;
    h1.Title='First order indexes';
    
    figure(),ht=heatmap(total_order); %Total order heatmap
    xlabel('Voltage magnitude error')
    ylabel('Measured value')
    ht.XDisplayLabels=h1.XDisplayLabels;
    ht.YDisplayLabels=h1.YDisplayLabels;
    ht.Colormap=parula;
    ht.Title='Total order indexes';
end
toc
%% Rankings
tic
stdV_mc=UAcode(Vindex,Sbrindex,time); %Compute extended uncertainty profile
%First metric
rank1=total_order*(stdV_mc.^2); %compute SI*unc^2
rank1=[ht.YDisplayLabels num2cell(rank1)];
rank1=sortrows(rank1,2,'descend');
first_order_cell=num2cell(total_order);

%Second metric
rank2=total_order*(stdV_mc.*(stdV_mc>0.85)).^2; %compute SI*unc^2 for unc>threshold (change 0.85 accordingly)
rank2=[ht.YDisplayLabels num2cell(rank2)];
rank2=sortrows(rank2,2,'descend');

%Third metric
feeders=cell2mat(grid_data.ff); 
feeders=[(2:99)' str2num(feeders(:,3))];
[~,locs] = findpeaks([min(stdV_mc);stdV_mc;min(stdV_mc)],'SortStr','descend','NPeaks',1); %Find feeder with highest uncertainty
locs=locs-1;
feeder_max=feeders(feeders(:,1)==locs,2);
idx=feeders(feeders(:,2)==feeder_max,1);
rank3=total_order(:,idx)*(stdV_mc(idx).^2);
rank3=[ht.YDisplayLabels num2cell(rank3)];
rank3=sortrows(rank3,2,'descend');

%Modified third metric
unc_sum=zeros(feeders(end,2),1);
for i=1:feeders(end,2)
    unc_sum(i)=sum(stdV_mc([false;feeders(:,2)==i]));
end
[~,feeder_max]=max(unc_sum);
idx=feeders(feeders(:,2)==feeder_max,1);
rank4=total_order(:,idx)*(stdV_mc(idx).^2);
rank4=[ht.YDisplayLabels num2cell(rank4)];
rank4=sortrows(rank4,2,'descend');

jj=[rank1 rank2 rank3 rank4]; %join rankings to compare
toc
