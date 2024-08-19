%% Overall main code
clc
clear all

% Read grid data for whole year
grid_data = Data_reader_industrial();
Base.V=15e3/sqrt(3); %In V
Base.S=1e6; % In VA

%% Pre data handling

%% Loop over one year (2030) with initial meter config
time.year=2030;
c=0; %Variable for retrieving the results of a specific timestep that can be reused in the current timestep
num_nodes=99;
Vyear.complex=zeros(num_nodes,96*365); % Matrix for storing the complex voltages troughout the year
Vyear.mag=zeros(num_nodes,96*365); % Matrix for storing the magnitude of voltages troughout the year
Vyear.phase=zeros(num_nodes,96*365); % Matrix for storing the phase of voltages troughout the year
Vyear.real=zeros(num_nodes,96*365); % Matrix for storing the real part of voltages troughout the year
Vyear.imag=zeros(num_nodes,96*365); % Matrix for storing the imaginary part of voltages troughout the year
% Main loop for the year
check_vector=[]; %Vector containing the months for which the first 7 days have already been simulated
counter=0; %Variable to check if a week has been calculated (All weeks in the month are the same)
for m=datenum(datetime(2030,01,01)):1/96:datenum(datetime(2030,12,31))+1-1/96
    c=c+1;
    if rem(c,100) == 0
        c
    end
    t=datevec(m);
    time.month=t(2);
    time.week=weekday(datetime(t));
    time.day=round(mod(m,floor(m))*96+1);
    if ismember(time.month,check_vector) %Check if the first week of the current month has already been calculated
        Vyear.complex(:,c)=Vyear.complex(:,c-7*96); %Reuse complex voltages if available
        Vyear.mag(:,c)=Vyear.mag(:,c-7*96); %Reuse magnitude of voltages if available
        Vyear.phase(:,c)=Vyear.phase(:,c-7*96); %Reuse phase of voltages if available
        Vyear.real(:,c)=Vyear.real(:,c-7*96); %Reuse real part of voltages if available
        Vyear.imag(:,c)=Vyear.imag(:,c-7*96); %Reuse imaginary part of voltages if available

    else
        %% Update grid data for current date
        [branch, node] = Network_industrial_Atlantide(Base, time, grid_data);
        branch.B = zeros(1,branch.num);      %B for the branches is set to zero
        [Y, Bs] = Ymatrix(branch, node);                          % computation admittance matrix of the network

        %% Power flow for current date
        % Run power flow
        [V, I, S, ~] = Power_flow_NV(branch, node, Y, Bs);

        %% Save results
        Vyear.complex(:,c)=V.complex.true_val; % Save complex voltages
        Vyear.mag(:,c)=V.mag.true_val; %Save magnitude of voltages
        Vyear.phase(:,c)=V.phase.true_val; %Save phase of voltages
        Vyear.real(:,c)=V.real.true_val; %Save real part of voltages
        Vyear.imag(:,c)=V.imag.true_val; %Save imaginary part of voltages
        counter=counter+1;
        if counter==7*96 % Check if the first seven days of the month have been completely calculated
            counter=0;
            check_vector=[check_vector time.month];
        end
    end
end
