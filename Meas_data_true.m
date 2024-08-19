function [zdatatrue] = Meas_data_true(branch, S, V, I)
             
%%% Computation number of measurement for each different category
numV = length(V.mag.meas.index);                                           % number voltage magnitude meas
numPi = length(S.inj.meas.index);                                          % number active pwr inj meas
numQi = length(S.inj.meas.index);                                          % number reactive pwr inj meas
numPf1 = length(S.br1.meas.index);                                         % number branch active pwr meas
numQf1 = length(S.br1.meas.index);                                         % number branch reactive pwr meas
numPf2 = length(S.br2.meas.index);                                         % number branch active pwr meas
numQf2 = length(S.br2.meas.index);                                         % number branch reactive pwr meas
numI = length(I.br1.mag.meas.index);                                           % number current magnitude meas
numVsync = length(V.mag.measPmu.index);                                    % number voltage PMU meas
numIsync = length(I.br1.mag.measPmu.index);                                    % number current PMU meas
num_meas = numPi + numQi + numPf1 + numQf1 + numPf2 + numQf2 + numV + numI + 2*numVsync + 2*numIsync;    % total number meas

%%% Definition codes associated to each type of measurement
t1 = ones(numV,1);                                                         % code = 1 -> voltage magnitude meas
t2 = 2*ones(numPi,1);                                                      % code = 2 -> active pwr inj meas
t3 = 3*ones(numQi,1);                                                      % code = 3 -> reactive pwr inj meas
t4 = 4*ones(numPf1+numPf2,1);                                              % code = 4 -> active branch pwr meas
t5 = 5*ones(numQf1+numQf2,1);                                              % code = 5 -> reactive branch pwr meas
t6 = 6*ones(numI,1);                                                       % code = 6 -> current magnitude meas
t7 = 7*ones(numVsync,1);                                                   % code = 7 -> PMU voltage magnitude meas
t8 = 8*ones(numVsync,1);                                                   % code = 8 -> PMU voltage phase angle meas
t9 = 9*ones(numIsync,1);                                                   % code = 9 -> PMU current magnitude meas
t10 = 10*ones(numIsync,1);                                                 % code = 10 -> PMU current phase angle meas

type = [t1; t2; t3; t4; t5; t6; t7; t8; t9; t10];                          %%% 1st column of zdata matrix

%%% Definition of the true values at the measurement points
V.mag.meas.true_val = V.mag.true_val(V.mag.meas.index);
S.inj.meas.P.true_val = S.inj.real.true_val(S.inj.meas.index);
S.inj.meas.Q.true_val = S.inj.imag.true_val(S.inj.meas.index);
S.br1.meas.P.true_val = S.br1.real.true_val(S.br1.meas.index);
S.br2.meas.P.true_val = S.br2.real.true_val(S.br2.meas.index);
S.br1.meas.Q.true_val = S.br1.imag.true_val(S.br1.meas.index);
S.br2.meas.Q.true_val = S.br2.imag.true_val(S.br2.meas.index);
I.br1.mag.meas.true_val = I.br1.mag.true_val(I.br1.mag.meas.index);
V.mag.measPmu.true_val = V.mag.true_val(V.mag.measPmu.index);
V.phase.measPmu.true_val = V.phase.true_val(V.phase.measPmu.index);
I.br1.mag.measPmu.true_val = I.br1.mag.true_val(I.br1.mag.measPmu.index);
I.br1.phase.measPmu.true_val = I.br1.phase.true_val(I.br1.phase.measPmu.index);

z = [V.mag.meas.true_val; S.inj.meas.P.true_val; S.inj.meas.Q.true_val; 
    S.br1.meas.P.true_val; S.br2.meas.P.true_val; S.br1.meas.Q.true_val; 
    S.br2.meas.Q.true_val; I.br1.mag.meas.true_val; V.mag.measPmu.true_val; 
    V.phase.measPmu.true_val; I.br1.mag.measPmu.true_val; I.br1.phase.measPmu.true_val];   %%% 2nd columm zdata matrix

%%% from = definition start_node of branch meas or node for either voltage or pwr inj meas
%%% to = definition end_node of branch meas (0 for voltage or pwr inj meas)
%%% br = definition branch code of the branch meas (0 for voltage or pwr inj meas)
V.mag.meas.from = V.mag.meas.index;
V.mag.meas.to = zeros(numV,1);
V.mag.meas.br = zeros(numV,1);
S.inj.meas.from = S.inj.meas.index;          
S.inj.meas.to = zeros(numPi,1);
S.inj.meas.br = zeros(numPi,1);
S.br1.meas.from = branch.start(S.br1.meas.index);
S.br1.meas.to = branch.end(S.br1.meas.index);
S.br1.meas.br = branch.cod(S.br1.meas.index);
S.br2.meas.from = branch.end(S.br2.meas.index);
S.br2.meas.to = branch.start(S.br2.meas.index);
S.br2.meas.br = branch.cod(S.br2.meas.index);
I.br1.mag.meas.from = branch.start(I.br1.mag.meas.index);
I.br1.mag.meas.to = branch.end(I.br1.mag.meas.index);
I.br1.mag.meas.br = branch.cod(I.br1.mag.meas.index);
V.mag.measPmu.from = V.mag.measPmu.index;
V.mag.measPmu.to = zeros(numVsync,1);
V.mag.measPmu.br = zeros(numVsync,1);
V.phase.measPmu.from = V.phase.measPmu.index;
V.phase.measPmu.to = zeros(numVsync,1);
V.phase.measPmu.br = zeros(numVsync,1);
I.br1.mag.measPmu.from = branch.start(I.br1.mag.measPmu.index);
I.br1.mag.measPmu.to = branch.end(I.br1.mag.measPmu.index);
I.br1.mag.measPmu.br = branch.cod(I.br1.mag.measPmu.index);
I.br1.phase.measPmu.from = branch.start(I.br1.phase.measPmu.index);
I.br1.phase.measPmu.to = branch.end(I.br1.phase.measPmu.index);
I.br1.phase.measPmu.br = branch.cod(I.br1.phase.measPmu.index);

br = [V.mag.meas.br; S.inj.meas.br; S.inj.meas.br; 
    S.br1.meas.br; -S.br2.meas.br; S.br1.meas.br; 
    -S.br2.meas.br; I.br1.mag.meas.br; V.mag.measPmu.br; 
    V.phase.measPmu.br; I.br1.mag.measPmu.br; I.br1.phase.measPmu.br];                            % 3rd column zdata matrix

from = [V.mag.meas.from; S.inj.meas.from; S.inj.meas.from; 
    S.br1.meas.from; S.br2.meas.from; S.br1.meas.from; 
    S.br2.meas.from; I.br1.mag.meas.from; V.mag.measPmu.from; 
    V.phase.measPmu.from; I.br1.mag.measPmu.from; I.br1.phase.measPmu.from];                      % 4th column zdata matrix

to = [V.mag.meas.to; S.inj.meas.to; S.inj.meas.to; 
    S.br1.meas.to; S.br2.meas.to; S.br1.meas.to; 
    S.br2.meas.to; I.br1.mag.meas.to; V.mag.measPmu.to; 
    V.phase.measPmu.to; I.br1.mag.measPmu.to; I.br1.phase.measPmu.to];                            % 5th column zdata matrix

%%% definition standard deviation for the different types of measurements    

if numV > 0
    V.mag.meas.std_dev = ((V.mag.meas.unc/300).*ones(numV,1)).*V.mag.meas.true_val;
else
    V.mag.meas.std_dev = [];
end
if numPi > 0
    S.inj.meas.P.std_dev = abs(((S.inj.meas.unc/300).*ones(numPi,1)).*S.inj.meas.P.true_val);
else
    S.inj.meas.P.std_dev = [];
end
if numQi > 0
    S.inj.meas.Q.std_dev = abs(((S.inj.meas.unc/300).*ones(numQi,1)).*S.inj.meas.Q.true_val);
else
    S.inj.meas.Q.std_dev = [];
end
if numPf1 > 0
    S.br1.meas.P.std_dev = abs(((S.br1.meas.unc/300).*ones(numPf1,1)).*S.br1.meas.P.true_val);
else
    S.br1.meas.P.std_dev = [];
end
if numQf1 > 0
    S.br1.meas.Q.std_dev = abs(((S.br1.meas.unc/300).*ones(numQf1,1)).*S.br1.meas.Q.true_val);
else
    S.br1.meas.Q.std_dev = [];
end
if numPf2 > 0
    S.br2.meas.P.std_dev = abs(((S.br2.meas.unc/300).*ones(numPf2,1)).*S.br2.meas.P.true_val);
else
    S.br2.meas.P.std_dev = [];
end
if numQf2 > 0
    S.br2.meas.Q.std_dev = abs(((S.br2.meas.unc/300).*ones(numQf2,1)).*S.br2.meas.Q.true_val);
else
    S.br2.meas.Q.std_dev = [];
end
if numI > 0
    I.br1.mag.meas.std_dev = ((I.br1.mag.meas.unc/300).*ones(numI,1)).*I.br1.mag.meas.true_val;
else
    I.br1.mag.meas.std_dev = [];
end
if numVsync > 0
    V.mag.measPmu.std_dev = ((V.mag.measPmu.unc/300).*ones(numVsync,1)).*V.mag.measPmu.true_val;
    V.phase.measPmu.std_dev = ((V.phase.measPmu.unc/300).*ones(numVsync,1));
else
    V.mag.measPmu.std_dev = [];
    V.phase.measPmu.std_dev = [];
end
if numIsync > 0
    I.br1.mag.measPmu.std_dev = ((I.br1.mag.measPmu.unc/300).*ones(numIsync,1)).*I.br1.mag.measPmu.true_val;
    I.br1.phase.measPmu.std_dev = ((I.br1.phase.measPmu.unc/300).*ones(numIsync,1));
else
    I.br1.mag.measPmu.std_dev = [];
    I.br1.phase.measPmu.std_dev = [];
end

if ~isempty(S.inj.meas.from)
    if S.inj.meas.from(1) == 1
        S.inj.meas.P.std_dev(1) = S.inj.meas.P.std_dev(1)*(S.br1.meas.unc/S.inj.meas.unc);
        S.inj.meas.Q.std_dev(1) = S.inj.meas.Q.std_dev(1)*(S.br1.meas.unc/S.inj.meas.unc);
    end
end
    

std_dev = [V.mag.meas.std_dev; S.inj.meas.P.std_dev; S.inj.meas.Q.std_dev; 
    S.br1.meas.P.std_dev; S.br2.meas.P.std_dev; S.br1.meas.Q.std_dev; 
    S.br2.meas.Q.std_dev; I.br1.mag.meas.std_dev; V.mag.measPmu.std_dev; 
    V.phase.measPmu.std_dev; I.br1.mag.measPmu.std_dev; I.br1.phase.measPmu.std_dev];


%%% Definition zdata matrix (with true measurement values)
zdatatrue = zeros(num_meas,6);
zdatatrue2 = zeros(num_meas,6);
zdatatrue2(:,1) = type;
zdatatrue2(:,2) = z;
zdatatrue2(:,3) = br;
zdatatrue2(:,4) = from;
zdatatrue2(:,5) = to;
zdatatrue2(:,6) = std_dev;


%%% Part created to alternate the presence of magnitude and phase angle meas in case of PMUs
Vsyncstart = 1 + numPi + numQi + numPf1 + numQf1 + numPf2 + numQf2 +numV + numI;
Isyncstart = Vsyncstart + 2*numVsync;
zdatatrue(1:Vsyncstart-1 , :) = zdatatrue2(1:Vsyncstart-1 , :);
for i=1:numVsync
    zdatatrue(Vsyncstart + 2*(i-1) , :) = zdatatrue2(Vsyncstart + i-1 , :);
    zdatatrue(Vsyncstart + 2*i - 1 , :) = zdatatrue2(Vsyncstart + i-1 + numVsync , :);
end
for i=1:numIsync
    zdatatrue(Isyncstart + 2*(i-1) , :) = zdatatrue2(Isyncstart + i-1 , :);
    zdatatrue(Isyncstart + 2*i - 1 , :) = zdatatrue2(Isyncstart + i-1 + numIsync , :);
end
end
