function [branch, P, Q, dg_lim] = Data_updater_industrial(Base, time, grid)

counter = 0;

for i = 1:length(grid.from)
    if strcmp(grid.from{i,1},'')
        continue
    end
    if strcmp(grid.state_from{i,1},'A') || strcmp(grid.state_to{i,1},'A')
        continue
    end
    counter = counter +1;
    num = grid.from{i,1}(end-2:end)-'0';
    %     branch.start(counter,1) = sscanf(num,'%d'); %sscanf
    branch.start(counter,1) = num(end)+num(end-1)*10+num(end-2)*100; %sscanf
    num = grid.to{i,1}(end-2:end)-'0';
    %     branch.end(counter,1) = sscanf(num,'%d'); %sscanf
    branch.end(counter,1) = num(end)+num(end-1)*10+num(end-2)*100; %sscanf
    branch.length(counter,1) = grid.len(i,1);
    type_now = grid.type{i,1};
    iter=find(strcmp(grid.cable,type_now),1);
%     check = 0;
%     iter = 0;
%     while check == 0
%         iter = iter+1;
%         if strcmp(type_now, grid.cable{iter,1}) %sscanf
%             check = 1;
%         else 
%             continue
%         end
        branch.R(counter,1) = grid.R(iter)*grid.len(i,1)/Base.Z; %per unit
        branch.X(counter,1) = 2*pi*50*grid.L(iter)/(Base.Z*1000); %per unit
%     end
end
branch.feeder = zeros(length(branch.start),1);
for i = 3:100
    idx = find(branch.end == i);
%     branch.feeder(idx,1) = sscanf(grid.ff{i-2,1}(end),'%d'); %sscanf
    branch.feeder(idx,1) = grid.ff{i-2,1}(end)-'0'; %sscanf
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

coeff.res.total = 0.8*grid.coeff.res.day(time.day)*grid.coeff.res.week(time.week)*grid.coeff.res.month(time.month)*grid.coeff.res.year(time.year);
coeff.ind.total = grid.coeff.ind.day(time.day)*grid.coeff.ind.week(time.week)*grid.coeff.ind.month(time.month)*grid.coeff.ind.year(time.year);
coeff.com.total = grid.coeff.com.day(time.day)*grid.coeff.com.week(time.week)*grid.coeff.com.month(time.month)*grid.coeff.com.year(time.year);
coeff.pv.total = grid.coeff.pv.day(time.day)*grid.coeff.pv.week(time.week)*grid.coeff.pv.month(time.month)*grid.coeff.pv.year(time.year);
coeff.chp.total = grid.coeff.chp.day(time.day)*grid.coeff.chp.week(time.week)*grid.coeff.chp.month(time.month)*grid.coeff.chp.year(time.year);
coeff.wind.total = grid.coeff.wind.day(time.day)*grid.coeff.wind.week(time.week)*grid.coeff.wind.month(time.month)*grid.coeff.wind.year(time.year);

P = zeros(100,1);
Q = zeros(100,1);
dg_lim=P;

for i = 1:length(grid.load_P)
    num = grid.load_bus{i,1}(end-2:end)-'0'; %sscanf
    num=num(end)+num(end-1)*10+num(end-2)*100;
    switch grid.load_type{i,1}
        case 'RES'
            load_Ppu = grid.load_P(i,1)*coeff.res.total/Base.S;
            load_Qpu = grid.load_Q(i,1)*coeff.res.total/Base.S;
        case 'IND'
            load_Ppu = grid.load_P(i,1)*coeff.ind.total/Base.S;
            load_Qpu = grid.load_Q(i,1)*coeff.ind.total/Base.S;
        case 'COM'
            load_Ppu = grid.load_P(i,1)*coeff.com.total/Base.S;
            load_Qpu = grid.load_Q(i,1)*coeff.com.total/Base.S;
    end
    P(num,1) = P(num,1) + load_Ppu;
    Q(num,1) = Q(num,1) + load_Qpu;
end


for i = 1:length(grid.gen_P_rotanti)
    num = sscanf(grid.gen_bus_rotanti{i,1}(end-2:end),'%d');
    switch grid.gen_type_rotanti{i,1}
        case 'CHP'
            gen_Ppu = grid.gen_P_rotanti(i,1)*coeff.chp.total/Base.S;
            gen_Qpu = grid.gen_Q_rotanti(i,1)*coeff.chp.total/Base.S;
            dg_lim(num)=dg_lim(num)+gen_Ppu;
        case 'WIND'
            gen_Ppu = grid.gen_P_rotanti(i,1)*coeff.wind.total/Base.S;
            gen_Qpu = grid.gen_Q_rotanti(i,1)*coeff.wind.total/Base.S;
            dg_lim(num)=dg_lim(num)+gen_Ppu;
    end
    P(num,1) = P(num,1) - gen_Ppu;
    Q(num,1) = Q(num,1) - gen_Qpu;
end

for i = 1:length(grid.gen_P_statici)
    num = sscanf(grid.gen_bus_statici{i,1}(end-2:end),'%d');
    gen_Ppu = grid.gen_P_statici(i,1)*coeff.pv.total/Base.S;
    dg_lim(num)=dg_lim(num)+gen_Ppu;
    P(num,1) = P(num,1) - gen_Ppu;
end

P(1:2,:) = [];
Q(1:2,:) = [];
dg_lim(1:2,:) = [];

P = (P/3)*1e6;
Q = (Q/3)*1e6;
dg_lim = (dg_lim/3)*1e6;