clc
clear

%% Initialization parameters and variable used to save the results
tic
mc_iter = 1000;        % Monte Carlo iterations for statistical analysis (put more than 1 if you want to carry out some statistical analysis via Monte Carlo simulation)
time_steps = 1440;  % number of considered time steps (in this code, 1 day with 1 minute resolution)
num_nodes = 99;     % number of nodes of the grid used below (TO BE CHANGED if a different grid is used)

Vtime_true = zeros(num_nodes,mc_iter);   %vector where the "true" voltage results are saved
Vtime_est = zeros(num_nodes,mc_iter);    %vector where the "actual" (coming from the SE + voltage control applications) voltage results are saved

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
    branch.B = zeros(1,branch.num);
    [Y, Bs] = Ymatrix(branch, node);                          % computation admittance matrix of the network

    % Power flow calculation for calculating the state assumed as reference and as true operating condition of the network;
    [V, I, S, num_iter] = Power_flow_NV(branch, node, Y, Bs);

    %% Start Monte Carlo procedure in case of statistical analysis
    branch2=branch;
    for MC=1:mc_iter
        if rem(MC,100) == 0
            MC
        end

        %% Meter placement
        V.mag.meas.index = [1]';               %indication of the nodes where voltage magnitude meas are placed.
        S.inj.meas.index = [2:1:node.num]';          % indication of the nodes for the power injection meas
        S.br1.meas.index = []';              % indication of the branches where active and reactive power meas are placed (measurement on the starting side of the branch).
        S.br2.meas.index = []';                % indication of the branches where active and reactive power meas are placed (measurement on the ending side of the branch).
        I.br1.mag.meas.index = []';            % indication of the branches where current magnitude meas are placed (there are no shunt admittances, thus the side of the branch is not important).
        I.br1.mag.measPmu.index = []';                                             % indication of the branches where PMU current meas are placed.
        I.br1.phase.measPmu.index = I.br1.mag.measPmu.index;                       % indication of the branches where PMU current meas are placed.
        V.mag.measPmu.index = []';                                                 % indication of the nodes where PMU voltage meas are placed.
        V.phase.measPmu.index = V.mag.measPmu.index;                                   % indication of the nodes where PMU voltage meas are placed.

        %% Uncertanty of the measurements
        V.mag.meas.unc = [1];           % uncertainty (in percent) of voltage magnitude measurements
        S.inj.meas.unc = [50];           % uncertainty (in percent) of active and reactive power injection measurements (pseudo-measurements)
        S.br1.meas.unc = [2];           % uncertainty (in percent) of active and reactive branch power measurements
        S.br2.meas.unc = 0;           % uncertainty (in percent) of active and reactive branch power measurements
        I.br1.mag.meas.unc = 0;       % uncertainty (in percent) of current magnitude measurements
        I.br1.mag.measPmu.unc = 0;                                               % uncertainty (in percent) of current magnitude measurements in PMUs
        I.br1.phase.measPmu.unc = 0;                                             % uncertainty (in centiradians) of current phase-angle measurements in PMUs
        V.mag.measPmu.unc = 0;                                                   % uncertainty (in percent) of voltage magnitude measurements in PMUs
        V.phase.measPmu.unc = 0;                                                 % uncertainty (in centiradians) of voltage phase-angle measurements in PMUs

        %% Creation of the measurements
        [zdatatrue] = Meas_data_true(branch, S, V, I);
        errn = randn(size(zdatatrue,1),1);
        [zdata] = Meas_data_generation(zdatatrue, errn);     % random generation of the measurements, starting from the true values, according to the distribution uncertainty of each measurement

        %% State Estimation
        [Y2,Bs2] = Ymatrix(branch2, node); % computation admittance matrix of the network with corrupted network parameters
        [V, I, S, num_iter, uncVm_percent] = NV_DSSE_call(branch2, node, zdata, Y2, Bs2, V, I, S); % function where the SE algorithm is called

        %% Saving of the results

        Vtime_true(:,MC) = V.mag.true_val';
        Vtime_est(:,MC) = V.mag.est';
        Vreal_est(:,MC) = V.real.est';
        Vimag_est(:,MC) = V.imag.est';
        Vreal_true(:,MC) = V.real.true_val';
        Vimag_true(:,MC) = V.imag.true_val';
        Vtime_true_ang(:,MC)=V.phase.true_val;
        Vtime_est_ang(:,MC)=V.phase.est;
    end
    stdV_mc = std(abs(Vtime_est) - V.mag.true_val*ones(1,mc_iter),0,2); %Computation of std of the errors of the magnitude estimates
    stdV_mc=300*stdV_mc./V.mag.true_val; %Calculation of extended uncertainty of the magnitude estimations
end
%%Plot uncertainty profile
figure(4)
if ~ishold()
    plot([1:99],ones(1,99),'DisplayName','1% limit') % 1% uncertianty threshold
end
hold on
legend()
plot(stdV_mc,'DisplayName','last plot'); %Plot current uncertainty profile
xlabel('Node')
ylabel('Voltage magnitude uncertainty (%)')
ylim([0,4.5])
toc