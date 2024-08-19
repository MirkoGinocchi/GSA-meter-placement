function [branch, P, Q] = Data_reader_rural(Base, time)

filename = 'C:\Users\Usuario\OneDrive - Universidad del Norte\Universidad+\Maestría\RWTH\3. Semester\Thesis\Sciebo\David\Code\Atlantide networks\ATL_Grid_Rural.xls';
[~,from] = xlsread(filename,'Linee','B2:B115');
[~,to] = xlsread(filename,'Linee','C2:C115');
[~,type] = xlsread(filename,'Linee','D2:D115');
len = xlsread(filename,'Linee','F2:F115');
[~,state_from] = xlsread(filename,'Linee','H2:H115');
[~,state_to] = xlsread(filename,'Linee','I2:I115');

filename2 = 'C:\Users\Usuario\OneDrive - Universidad del Norte\Universidad+\Maestría\RWTH\3. Semester\Thesis\Sciebo\David\Code\Atlantide networks\ATL_Components_Rural_Grid.xls';
[~,cable] = xlsread(filename2,'Linee','A2:A15');
R = xlsread(filename2,'Linee','C2:C15');
L = xlsread(filename2,'Linee','D2:D15');
counter = 0;

for i = 1:length(from)
    if strcmp(from{i,1},'')
        continue
    end
    if strcmp(state_from{i,1},'A') || strcmp(state_to{i,1},'A')
        continue
    end
    counter = counter +1;
    num = from{i,1}(end-2:end);
    branch.start(counter,1) = str2double(num);
    num = to{i,1}(end-2:end);
    branch.end(counter,1) = str2double(num);
    branch.length(counter,1) = len(i,1);
    type_now = type{i,1};
    check = 0;
    iter = 0;
    while check == 0
        iter = iter+1;
        if strcmp(type_now, cable{iter,1})
            check = 1;
        else 
            continue
        end
        branch.R(counter,1) = R(iter)*len(i,1)/Base.Z; %per unit
        branch.X(counter,1) = 2*pi*50*L(iter)/(Base.Z*1000); %per unit
    end
end

[~,ff] = xlsread(filename,'Nodi','F2:F104');
branch.feeder = zeros(length(branch.start),1);
for i = 3:103
    idx = find(branch.end == i);
    branch.feeder(idx,1) = str2double(ff{i-2,1}(end));
end

branch.start = branch.start - 1;
branch.end = branch.end - 1;

[branch.end, idx] = sort(branch.end);
branch.start = branch.start(idx);
branch.R = branch.R(idx);
branch.X = branch.X(idx);
branch.length = branch.length(idx);
branch.feeder = branch.feeder(idx);

if time.year < 2010
    time.year = 2010;
elseif time.year > 2030
    time.year = 2030;
else
    time.year = time.year - 2009;
end

filename3 = 'C:\Users\Usuario\OneDrive - Universidad del Norte\Universidad+\Maestría\RWTH\3. Semester\Thesis\Sciebo\David\Code\Atlantide networks\ATL_Profiles_Rural_Grid.xls';
coeff.agr.day = xlsread(filename3, 'Profili Carico Giornalieri', 'F2:F97');
coeff.rlv.day = xlsread(filename3, 'Profili Carico Giornalieri', 'G2:G97');
coeff.rmv1.day = xlsread(filename3, 'Profili Carico Giornalieri', 'H2:H97');
coeff.rmv2.day = xlsread(filename3, 'Profili Carico Giornalieri', 'I2:I97');
coeff.rmv3.day = xlsread(filename3, 'Profili Carico Giornalieri', 'J2:J97');

coeff.agr.week = xlsread(filename3, 'Coeff Sett Carico', 'F2:F8');
coeff.rlv.week = xlsread(filename3, 'Coeff Sett Carico', 'G2:G8');
coeff.rmv1.week = xlsread(filename3, 'Coeff Sett Carico', 'H2:H8');
coeff.rmv2.week = xlsread(filename3, 'Coeff Sett Carico', 'I2:I8');
coeff.rmv3.week = xlsread(filename3, 'Coeff Sett Carico', 'J2:J8');

coeff.agr.month = xlsread(filename3, 'Coeff Mens Carico', 'F2:F13');
coeff.rlv.month = xlsread(filename3, 'Coeff Mens Carico', 'G2:G13');
coeff.rmv1.month = xlsread(filename3, 'Coeff Mens Carico', 'H2:H13');
coeff.rmv2.month = xlsread(filename3, 'Coeff Mens Carico', 'I2:I13');
coeff.rmv3.month = xlsread(filename3, 'Coeff Mens Carico', 'J2:J13');

coeff.agr.year = xlsread(filename3, 'Coeff Ann Carico', 'F2:F22');
coeff.rlv.year = xlsread(filename3, 'Coeff Ann Carico', 'G2:G22');
coeff.rmv1.year = xlsread(filename3, 'Coeff Ann Carico', 'H2:H22');
coeff.rmv2.year = xlsread(filename3, 'Coeff Ann Carico', 'I2:I22');
coeff.rmv3.year = xlsread(filename3, 'Coeff Ann Carico', 'J2:J22');

coeff.pv.day = xlsread(filename3, 'Profili Gen Giornalieri', 'D2:D97');

coeff.pv.week = xlsread(filename3, 'Coeff Sett Gen', 'D2:D8');

coeff.pv.month = xlsread(filename3, 'Coeff Mens Gen', 'D2:D13');

coeff.pv.year = xlsread(filename3, 'Coeff Ann Gen', 'D2:D22');

coeff.agr.total = coeff.agr.day(time.day)*coeff.agr.week(time.week)*coeff.agr.month(time.month)*coeff.agr.year(time.year);
coeff.rlv.total = coeff.rlv.day(time.day)*coeff.rlv.week(time.week)*coeff.rlv.month(time.month)*coeff.rlv.year(time.year);
coeff.rmv1.total = coeff.rmv1.day(time.day)*coeff.rmv1.week(time.week)*coeff.rmv1.month(time.month)*coeff.rmv1.year(time.year);
coeff.rmv2.total = coeff.rmv2.day(time.day)*coeff.rmv2.week(time.week)*coeff.rmv2.month(time.month)*coeff.rmv2.year(time.year);
coeff.rmv3.total = coeff.rmv3.day(time.day)*coeff.rmv3.week(time.week)*coeff.rmv3.month(time.month)*coeff.rmv3.year(time.year);
coeff.pv.total = coeff.pv.day(time.day)*coeff.pv.week(time.week)*coeff.pv.month(time.month)*coeff.pv.year(time.year);

P = zeros(103,1);
Q = zeros(103,1);

[~,load_bus] = xlsread(filename, 'Carichi_A', 'B2:B191');
load_P = xlsread(filename, 'Carichi_A', 'C2:C191');
load_Q = xlsread(filename, 'Carichi_A', 'D2:D191');
[~,load_type] = xlsread(filename, 'Carichi_A', 'F2:F191');
for i = 1:length(load_P)
    num = str2double(load_bus{i,1}(end-2:end));
    switch load_type{i,1}
        case 'AGR'
            load_Ppu = load_P(i,1)*coeff.agr.total/Base.S;
            load_Qpu = load_Q(i,1)*coeff.agr.total/Base.S;
        case 'RLV'
            load_Ppu = load_P(i,1)*coeff.rlv.total/Base.S;
            load_Qpu = load_Q(i,1)*coeff.rlv.total/Base.S;
        case 'RMV_CUST1'
            load_Ppu = load_P(i,1)*coeff.rmv1.total/Base.S;
            load_Qpu = load_Q(i,1)*coeff.rmv1.total/Base.S;
        case 'RMV_CUST2'
            load_Ppu = load_P(i,1)*coeff.rmv2.total/Base.S;
            load_Qpu = load_Q(i,1)*coeff.rmv2.total/Base.S;
        case 'RMV_CUST3'
            load_Ppu = load_P(i,1)*coeff.rmv3.total/Base.S;
            load_Qpu = load_Q(i,1)*coeff.rmv3.total/Base.S;
    end
    P(num,1) = P(num,1) + load_Ppu;
    Q(num,1) = Q(num,1) + load_Qpu;
end

[~,gen_bus] = xlsread(filename, 'Generatori Statici_A', 'B2:B6');
gen_P = xlsread(filename, 'Generatori Statici_A', 'C2:C6');
for i = 1:length(gen_P)
    num = str2double(gen_bus{i,1}(end-2:end));
    gen_Ppu = gen_P(i,1)*coeff.pv.total/Base.S;
    P(num,1) = P(num,1) - gen_Ppu;
end

P(1:2,:) = [];
Q(1:2,:) = [];

P = (P/3)*1e6;
Q = (Q/3)*1e6;