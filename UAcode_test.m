function stdV_mc = UAcode_test(Vindex,Sbrindex)

%% Initialization parameters and variable used to save the results

s = RandStream.create('mt19937ar','seed',152489);                          % random extraction starting from the specified seed
RandStream.setGlobalStream(s);

mc_iter = 1000;        % Monte Carlo iterations for statistical analysis (put more than 1 if you want to carry out some statistical analysis via Monte Carlo simulation)
ctrl_steps = 1;    % number of iterations considered for the voltage control at each time step (to allow it to converge to the final result)
% time_steps = 1440;  % number of considered time steps (in this code, 1 day with 1 minute resolution)
num_nodes = 10;     % number of nodes of the grid used below (TO BE CHANGED if a different grid is used)

Vtime_true = zeros(num_nodes,mc_iter);   %c vector where the "true" voltage results are saved
Vtime_est = zeros(num_nodes,mc_iter);    %c vector where the "actual" (coming from the SE + voltage control applications) voltage results are saved

%% Start computations over the decided time steps
for t = 1:1    %%% --> change this if you want just to focus on a specific time step

    [branch, node] = Network_toy_grid();
    [Y, Bs] = Ymatrix(branch, node);                          % computation admittance matrix of the network
    [V, I, S, ~] = Power_flow_NV(branch, node, Y, Bs);

    %% Start Monte Carlo procedure in case of statistical analysis
    branch2=branch;
    for MC=1:mc_iter
        if rem(MC,100) == 0
            disp(MC)
        end

        %% Iterations over the voltage control algorithm
        for k=1:ctrl_steps
            %% Meter placement
            V.mag.meas.index = Vindex';               % indication of the nodes where voltage magnitude meas are placed.
%             V.mag.meas.index = [1:1:node.num]';               % indication of the nodes where voltage magnitude meas are placed.
            S.inj.meas.index = (2:1:node.num)';          % indication of the nodes for the power injection meas
            S.br1.meas.index = Sbrindex';               % indication of the branches where active and reactive power meas are placed (measurement on the starting side of the branch).
            S.br2.meas.index = []';                % indication of the branches where active and reactive power meas are placed (measurement on the ending side of the branch).
            I.br1.mag.meas.index = []';            % indication of the branches where current magnitude meas are placed (there are no shunt admittances, thus the side of the branch is not important).
            I.br1.mag.measPmu.index = []';                                             % indication of the branches where PMU current meas are placed.
            I.br1.phase.measPmu.index = I.br1.mag.measPmu.index;                       % indication of the branches where PMU current meas are placed.
            V.mag.measPmu.index = []';                                                 % indication of the nodes where PMU voltage meas are placed.
            V.phase.measPmu.index = V.mag.measPmu.index;                                   % indication of the nodes where PMU voltage meas are placed.

            %% Uncertanty of the measurements
            V.mag.meas.unc = 1;           % uncertainty (in percent) of voltage magnitude measurements
%             V.mag.meas.unc = [1;50;50;50;50;50;1;50;50;1];          % uncertainty (in percent) of voltage magnitude measurements
            S.inj.meas.unc = [50;50;50;50;50;50;50;50;50];           % uncertainty (in percent) of active and reactive power injection measurements (pseudo-measurements)
            S.br1.meas.unc = 1;%[500;500;500;500;500;500;500;500;500];           % uncertainty (in percent) of active and reactive branch power measurements
            S.br2.meas.unc = 0;           % uncertainty (in percent) of active and reactive branch power measurements
            I.br1.mag.meas.unc = 0;       % uncertainty (in percent) of current magnitude measurements
            I.br1.mag.measPmu.unc = 0;                                               % uncertainty (in percent) of current magnitude measurements in PMUs
            I.br1.phase.measPmu.unc = 0;                                             % uncertainty (in centiradians) of current phase-angle measurements in PMUs
            V.mag.measPmu.unc = 0;                                                   % uncertainty (in percent) of voltage magnitude measurements in PMUs
            V.phase.measPmu.unc = 0;                                                 % uncertainty (in centiradians) of voltage phase-angle measurements in PMUs

            %% Creation of the measurements
            [zdatatrue] = Meas_data_true(branch, S, V, I);
            if k == 1
                errn = randn(size(zdatatrue,1),1);
            end
            [zdata] = Meas_data_generation(zdatatrue, errn);     % random generation of the measurements, starting from the true values, according to the distribution uncertainty of each measurement
            
            %% State Estimation
            r=0.01*[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
            x=0.01*[0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
            for jj=1:9
                branch2.R(jj)=unifrnd((1-r(jj))*branch.R(jj),(1+r(jj))*branch.R(jj));
                branch2.X(jj)=unifrnd((1-x(jj))*branch.X(jj),(1+x(jj))*branch.X(jj));
            end
            [Y2,Bs2] = Ymatrix(branch2, node);
            [V, I, S, ~, ~] = NV_DSSE_call(branch2, node, zdata, Y2, Bs2, V, I, S); %Bs2 Branch2    % function where the SE algorithm is called

            %% Saving of the results
            if k == ctrl_steps
                Vtime_true(:,MC) = V.mag.true_val'; %changed
                Vtime_est(:,MC) = V.mag.est'; %changed
            end
        end
    end
    stdV_mc = std(abs(Vtime_est) - V.mag.true_val*ones(1,mc_iter),0,2);
    stdV_mc=300*stdV_mc./V.mag.true_val;
end
