%% Save SI from current sim
load('Sensitivity_Indeces_all.mat')
first_order_all.wst.OV.V1_32_50_84_93_S62=first_order;
total_order_all.wst.OV.V1_32_50_84_93_S62=total_order;
save("Sensitivity_Indeces_all","first_order_all","total_order_all")

%% Load SI from OV and UV and create the combined ranking
tic
%load SIs
load('Sensitivity_Indeces_all.mat')
load('inputs.mat')
load('feeders.mat')
%define timestep for under and overvoltage
time.year=2030;
time.month=6;%3-6clc
time.week=7;%4-7
time.day=47;%33-47
time_u.year=2030;
time_u.month=3;%3-6
time_u.week=4;%4-7
time_u.day=33;%33-47
%define SIs for under and overvoltage
total_order=total_order_all.wst.OV.V1_32_50_84_S62;%--------------------------------------------
total_order_u=total_order_all.wst.OV.V1_32_50_84_S62;%--------------------------------------------
%set measuremnt configuration
Vindex=[1,32,50,84]; %--------------------------------------------
Sbrindex=[62]; %--------------------------------------------
%Compute uncertainty profile for under and over voltage
stdV_mc=UAcode(Vindex,Sbrindex,time);
stdV_mc_u=UAcode(Vindex,Sbrindex,time_u);
%Calculate ws rankings for over and under voltage
rank=total_order*(stdV_mc).^2;
rank_u=total_order_u*(stdV_mc_u).^2;
% Join rankings and sort
wo=1;%2474/35040; %2409(2025)
wu=1;%1-wo;
rank=wo*rank+wu*rank_u;
rank=[s num2cell(rank)];
rank=sortrows(rank,2,'descend');

%Calculate ws (threshold) rankings for over and under voltage
rank2=total_order*(stdV_mc.*(stdV_mc>1)).^2;
rank2_u=total_order_u*(stdV_mc_u.*(stdV_mc_u>1)).^2;
% Join rankings and sort
rank2=wo*rank2+wu*rank2_u;
rank2=[s num2cell(rank2)];
rank2=sortrows(rank2,2,'descend');

% %Calculate fp rankings for over and under voltage
[~,locs] = findpeaks([min(stdV_mc);stdV_mc;min(stdV_mc)],'SortStr','descend','NPeaks',1);
locs=locs-1;
feeder_max=feeders(feeders(:,1)==locs,2);
idx=feeders(feeders(:,2)==feeder_max,1);
rank3=total_order(:,idx)*(stdV_mc(idx).^2);

[~,locs_u] = findpeaks([min(stdV_mc_u);stdV_mc_u;min(stdV_mc_u)],'SortStr','descend','NPeaks',1);
locs_u=locs_u-1;
feeder_max_u=feeders(feeders(:,1)==locs_u,2);
idx_u=feeders(feeders(:,2)==feeder_max_u,1);
rank3_u=total_order_u(:,idx_u)*(stdV_mc_u(idx_u).^2);

% % Join rankings and sort
rank3=wo*rank3+wu*rank3_u;
rank3=[s num2cell(rank3)];
rank3=sortrows(rank3,2,'descend');

jj=[rank rank2 rank3];



ranko=[s num2cell(rank)];
ranko=sortrows(ranko,2,'descend');

rank2o=[s num2cell(rank2)];
rank2o=sortrows(rank2o,2,'descend');

rank3o=[s num2cell(rank3)];
rank3o=sortrows(rank3o,2,'descend');

jjo=[ranko rank2o rank3o];

ranku=[s num2cell(rank_u)];
ranku=sortrows(ranku,2,'descend');

rank2u=[s num2cell(rank2_u)];
rank2u=sortrows(rank2u,2,'descend');

rank3u=[s num2cell(rank3_u)];
rank3u=sortrows(rank3u,2,'descend');

jju=[ranku rank2u rank3u];
figure()
hold on
plot(stdV_mc_u)
plot(stdV_mc)
ylim([0,4.5])
toc

%% heatmaps for document
% load('inputs.mat')
% load('Sensitivity_Indeces_all.mat')
% total_order=total_order_all.ws.y2018.V1;
% figure(),ht=heatmap(total_order);
%     xlabel('Voltage magnitude')
%     ylabel('Measured value')
%     for ii=1:size(total_order,2)
%         ht.XDisplayLabels{ii,1}=['V',num2str(ii)];
%     end
%     ht.YDisplayLabels=s;
%     ht.Colormap=parula;
%     ht.Title='Total order indexes';
%% get idx for meter_configs and set them
% load('meter_configs_ws_3.mat');
% asd=fields(first_order_all.wst.U1);
% for i=1:length(asd)
%     asd{i}(asd{i}=='_')=';';
% end
% asd
% 
% meter_configs.V.index{1}=[];
% meter_configs.V.index{2}=[];
% meter_configs.V.index{3}=[];
% meter_configs.V.index{4}=[];
% meter_configs.V.index{5}=[];
% meter_configs.V.index{6}=[];
% meter_configs.V.index{7}=[];
% meter_configs.V.index{8}=[];
% meter_configs.V.index{9}=[];
% meter_configs.V.index{10}=[];
% 
% meter_configs.Sbr.index{1}=[];
% meter_configs.Sbr.index{2}=[];
% meter_configs.Sbr.index{3}=[];
% meter_configs.Sbr.index{4}=[];
% meter_configs.Sbr.index{5}=[];
% meter_configs.Sbr.index{6}=[];
% meter_configs.Sbr.index{7}=[];
% meter_configs.Sbr.index{8}=[];
% meter_configs.Sbr.index{9}=[];
% meter_configs.Sbr.index{10}=[];

%% add feeder lines
% y=[0 4.5];
% x=[2,2;34,34;56,56;59,59;79,79;85,85;94,94];
% for i=1:2
%     figure(i)
%     hold on
%     plot(x(1,:),y,'k',x(2,:),y,'k',x(3,:),y,'k',x(4,:),y,'k',x(5,:),y,'k',x(6,:),y,'k',x(7,:),y,'k')
% end

%% save year results
% save("year wst","Vyear_all","Pess_all","Pdg_all","Qdg_all","CF_vector")

%% Animation for ten steps
% open('..\..\Presentation\10 steps ws 5-50%.fig')
% h=gca;
% h.Children(19).Visible='on';
% for i=9:18
%     h.Children(i).Visible='off';
% end
% pause()
% for i=18:-1:9
%     h.Children(i).Visible='on';
%     pause()
%     h.Children(i+1).Visible='off';
% end
% for i=9:19
%     h.Children(i).Visible='on';
% end

%% ws with normalized SI over inputs
% load('Sensitivity_Indeces_all.mat')
% asd=fields(first_order_all.ws.U3_50);
% f=openfig('..\..\Presentation\10 steps ws 3-50%.fig','invisible');
% load('inputs.mat')
% 
% for i=1:10
%     total_order=total_order_all.ws.U3_50.(asd{i});
%     total_norm=normalize(total_order,2,'Range');
%     total_norm(isnan(total_norm))=0;
%     stdV_mc2=f.Children(2).Children(round(-10*(i)/9+20.111111)).YData';
%     ranks=total_norm*(stdV_mc2.^2);
%     rank{i}=[s num2cell(ranks)];
%     rank{i}=sortrows(rank{i},2,'descend');
% end
% 
% [rank{1} rank{2} rank{3} rank{4} rank{5} rank{6} rank{7} rank{8} rank{9} rank{10}];

%% Recalculate CF_vectors
% load('year fp 3_50.mat')
% fp=[];
% for i=1:size(Pdg_all,3)
%     Pdg=Pdg_all(:,:,i);
%     Pdg(Pdg<1e-4)=0;
%     Pess=Pess_all(Pess_all(:,:,i)>1e-4,i);
%     fp=[fp sum(sum(abs(Pess)+abs(Pdg),2))];
% end

% load('year fp 3_50.mat')
CF_vectors.wsl=CF_vector;
Pdg_alls.wsl=Pdg_all;
Qdg_alls.wsl=Qdg_all;
Pess_alls.wsl=Pess_all;
Vyear_alls.wsl=Vyear_all;
% 
a=[CF_vectors.wss;CF_vectors.wst;CF_vectors.fps;CF_vectors.fss;CF_vectors.wsc;CF_vectors.fpc;CF_vectors.wsl];
a(1:6,1)=253747;
figure
bar([4e4/15*ones(5,1),-diff(a,1,2)'])
ylabel('Savings (EUR/yr)')
xticklabels({'1st meter','2nd meter','3rd meter','4th meter','5th meter'})
legend({'Cost of meter','1^{st} metric - OV','2^{nd} metric - OV','3^{rd} metric - OV','Modified 3^{rd} metric - OV','1^{st} and 2^{nd} metric - OV+UV','3^{rd} metric - OV+UV','1^{st} metric - OV 2018'})
title('Cost of meter vs savings')
ylim([0,12e3])
%%
close(figure(3))
figure(3)
titles={'1^{st} metric - OV','2^{nd} metric - OV','3^{rd} metric - OV','Modified 3^{rd} metric - OV','1^{st} and 2^{nd} metric - OV+UV','3^{rd} metric - OV+UV','1^{st} metric - OV 2018'};
a2=[a(1:5,:); a(7,:)];
colors=[0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;0.6350 0.0780 0.1840];
fs=16;
for i=1:7
    if i==7
    subplot(4,2,7.5)
    else
    subplot(4,2,i)
    end
    bar(-diff(a(i,:),1,2)',0.4,'FaceColor',colors(i,:))
    grid on
    hold on
    plot([0.5 5.5],4e4/15*[1 1],'--k','LineWidth',2)
    ylim([0,12e3])
    ylabel('Savings (EUR/yr)','FontSize', fs)
    xticklabels({'1st meter','2nd meter','3rd meter','4th meter','5th meter'})
    title(titles(i),'FontSize', fs)
    ax = gca;
    ax.FontSize = fs;
end
%%
% 
% FP=a*4/100;
% figure()
% bar(a'*4/100)
% ylabel('Flexible power (MWh/yr)')
% xticklabels({'Initial configuration','1st meter','2nd meter','3rd meter','4th meter','5th meter'})
% legend({'Weighted sum single','Feeder peak single','Feeder sum single','Weighted sum combined','Feeder peak combined'})
% title('Requested flexible power per year for different meter configurations')

TC=[CF_vectors.wss(end-1) 4*4e4/15;
    CF_vectors.wst(end-1) 3*4e4/15;
    CF_vectors.fps(end-1) 4*4e4/15;
    CF_vectors.fss(end-3) 2*4e4/15;
    CF_vectors.wsc(end-2) 3*4e4/15;
    CF_vectors.fpc(end-1) 4*4e4/15;
    CF_vectors.wsl(end-3) 2*4e4/15];
figure()
b=bar(TC,'stacked');
ylabel('Cost (EUR/yr)')
xticklabels({'1^{st} metric - OV','2^{nd} metric - OV','3^{rd} metric - OV','Modified 3^{rd} metric - OV','1^{st} and 2^{nd} metric - OV+UV','3^{rd} metric - OV+UV','1^{st} metric - OV 2018'})
legend({'Cost of flexible power','Cost of meters'})
title('Final total cost for different methods')

text(1:7,TC(:,1)/2,num2str(round(TC(:,1)),7),'HorizontalAlignment','center','Color','w','FontSize',18)
text(1:7,(TC(:,1)+sum(TC,2))/2,num2str(round(TC(:,2)),7),'HorizontalAlignment','center','Color','w','FontSize',18)
text(1:7,b(2).YEndPoints,num2str(round(sum(TC,2)),7),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',18)
grid on

mc=(0:5)*4e4/15;
CC1=[106000 CF_vectors.wss(2:6); mc(1:6)]';
CC2=[106000 CF_vectors.wst(2:5); mc(1:5)]';
CC3=[106000 CF_vectors.fps(2:6); mc(1:6)]';
CC4=[106000 CF_vectors.fss(2:4); mc(1:4)]';
CC5=[106000 CF_vectors.wsc(2:5); mc(1:5)]';
CC6=[106000 CF_vectors.fpc(2:6); mc(1:6)]';
CC7=[CF_vectors.wsl(1:4); mc(1:4)]';
figure()
hold on
b1=bar([-0.6,-0.4,-0.2,0,0.2,0.4]+1,CC1,'stacked');
b2=bar([-0.6,-0.4,-0.2,0,0.2]+2.4,CC2,'stacked');
b3=bar([-0.6,-0.4,-0.2,0,0.2,0.4]+3.6,CC3,'stacked');
b4=bar([-0.4,-0.2,0,0.2]+4.8,CC4,'stacked');
b5=bar([-0.6,-0.4,-0.2,0,0.2]+6,CC5,'stacked');
b6=bar([-0.6,-0.4,-0.2,0,0.2,0.4]+7.2,CC6,'stacked');
b7=bar([-0.4,-0.2,0,0.2]+8.4,CC7,'stacked');
ylabel('Cost (EUR/yr)')
xticks([1 2.2 3.4 4.6 5.8 7.2 8.4])
xticklabels({'1^{st} metric - OV', ...
    '2^{nd} metric - OV', ...
    '3^{rd} metric - OV', ...
    'Modified 3^{rd} metric - OV', ...
    '1^{st} and 2^{nd} metric - OV+UV', ...
    '3^{rd} metric - OV+UV', ...
    '1^{st} metric - OV 2018'})
% figure(3), legend({'Cost of meter','1^{st} and 2^{nd} metric – OV','3^{rd} metric – OV','Modified 3^{rd} metric – OV','1^{st} and 2^{nd} metric – OV+UV','3^{rd} metric – OV+UV'})
legend({'Cost of flexible power','Cost of meters'})
title('Total cost at each step for different methods')
grid on

%% CDF and clusters for year profile
% clc, clearvars,close all
openfig('\Voltage_profile_year.fig');
h=gca;
for i=1:99
    Vyear(i,:)=h.Children(102-i).YData;
end
% [f,x] = ecdf(Vyear(1,:));
% opts = statset('Display','off','UseParallel',true,'MaxIter',100);
% [~,b1]=min(min(Vyear,[],1));
% [~,b2]=max(max(Vyear,[],1));
% init=Vyear(:,[b1,b2]);
% d1=pdist2(Vyear',init','squaredeuclidean');
% d2=pdist2(Vyear',init');
% init3d=repmat(init',[1 1 10]);
% [idx,C,sumd,d,midx,info] = kmedoids(Vyear',4,'Replicate',10,'Options',opts);
% figure()
% figure(2) %6177 14543
% plot(Vyear)
% figure(3)
% colororder([0 0 1
%     1 0 0
%     0 1 0
%     0 0 0]);
% plot(Vyear(:,midx))
% figure(4)
% hold on
% for i=1:size(Vyear,2)
%     switch idx(i)
%         case 1
%             color='b';
%         case 2
%             color='r';
%         case 3 
%             color='g';
%         case 4
%             color='k';
%     end
%     plot(Vyear(:,i),color)
% end
% 
% figure(5)
% hold on
% for i=1:size(Vyear,2)
%     switch idx(i)
%         case 1
%             color='b';
%         case 2
%             color='r';
%         case 3 
%             color='g';
%         case 4
%             color='k';
%     end
%     plot(i*ones(99,1),Vyear(:,i),color)
% end
% 
mm=[];
for m=datenum(datetime(2015,01,01)):1/96:datenum(datetime(2015,12,31))+1-1/96
    mm=[mm m];
end
% leg={};
% for i=1:length(midx)
    asd=mm(midx(i));
    t=datevec(asd);
    time.month=t(2);
    time.week=weekday(datetime(t));
    time.day=round(mod(asd,floor(asd))*96+1);
%     leg{i}=['month:',num2str(time.month),' weekday:',num2str(time.week),' time:',num2str(t(end-2)),':',num2str(t(end-1))];
% end
% figure(4)
% legend(leg)
% figure(4)
% hold on
% plot(Vyear(:,b2),'r','DisplayName','Worst OV')
% plot(Vyear(:,b1),'k','DisplayName','Worst UV')
% title('silhouette for 8 clusters')

%% Combining P and Q
clc
%load SIs
load('Sensitivity_Indeces_all.mat')
load('inputs.mat')
load('feeders.mat')
%define timestep 
time.year=2030;
time.month=6;
time.week=7;
time.day=47;
%define SIs and check if first_order=total order
first_order=first_order_all.fs.U3_50.V1_7_14_19_30_32_49_S58_62_84;%--------------------------------------------
total_order=total_order_all.fs.U3_50.V1_7_14_19_30_32_49_S58_62_84;%--------------------------------------------
sum(sum(abs(first_order-total_order)))
%Add P and Q SIs
total_order(100:197,:)=total_order(100:197,:)+total_order(198:end,:);
%Remove Q SIs and Q labels
total_order(198:end,:)=[];
s(198:end,:)=[];
%set measuremnt configuration
Vindex=[1,7,14,19,30,32,49]; %--------------------------------------------
Sbrindex=[58,62,84]; %--------------------------------------------
%Compute uncertainty profile 
stdV_mc=UAcode(Vindex,Sbrindex,time);
%Calculate ws ranking
rank=total_order*(stdV_mc).^2;
rank=[s num2cell(rank)];
rank=sortrows(rank,2,'descend');

%Calculate ws (threshold) ranking
rank2=total_order*(stdV_mc.*(stdV_mc>1)).^2;
rank2=[s num2cell(rank2)];
rank2=sortrows(rank2,2,'descend');

% %Calculate fp ranking
[~,locs] = findpeaks([min(stdV_mc);stdV_mc;min(stdV_mc)],'SortStr','descend','NPeaks',1);
locs=locs-1;
feeder_max=feeders(feeders(:,1)==locs,2);
idx=feeders(feeders(:,2)==feeder_max,1);
rank3=total_order(:,idx)*(stdV_mc(idx).^2);
rank3=[s num2cell(rank3)];
rank3=sortrows(rank3,2,'descend');

%Calculate fs ranking
unc_sum=zeros(feeders(end,2),1);
for i=1:feeders(end,2)
    unc_sum(i)=sum(stdV_mc([false;feeders(:,2)==i]));
end
[~,feeder_max]=max(unc_sum);
idx=feeders(feeders(:,2)==feeder_max,1);
rank4=total_order(:,idx)*(stdV_mc(idx).^2);
rank4=[s num2cell(rank4)];
rank4=sortrows(rank4,2,'descend');


jj=[rank rank2 rank3 rank4];