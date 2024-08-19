function [branch_samples,node_samples,zdata_samples]=sampling_call(branch,node,zdata,num_samples)

%Function that makes the call for sampling each of the inputs to sample

%% Branch sampling (R,X,B)
%R, X and B uncertainties fixed or read? 
%Fixed: set up here and passed to the sampling function (below as fixed for
%testing)
%Read: add to the branch variable a matrix with the uncertainties
%Do somenthing and only pass a matrix "X" with two columns, mean and pu above and below mean 
X=[branch.R 0.1*ones(length(branch.R),1)];
raw_samples_R=sampling(X,num_samples,1); %sampling R
X=[branch.X 0.1*ones(length(branch.X),1)];
raw_samples_X=sampling(X,num_samples,1); %Sampling X
if max(branch.B)==0
    raw_samples_B=zeros(num_samples,length(branch.B));
else
    X=[branch.B' 0.1*ones(length(branch.B),1)];
    raw_samples_B=sampling(X,num_samples,1); %Sampling B (if not considered, then do not sample)
end
%Convert matrix with samples to struct "branch_samples" to be used in the
%main code
raw_samples_R(isnan(raw_samples_R))=0;
raw_samples_X(isnan(raw_samples_X))=0;
raw_samples_B(isnan(raw_samples_B))=0;
branch_samples.R=raw_samples_R';
branch_samples.X=raw_samples_X';
branch_samples.B=raw_samples_B;

%% Node sampling (P,Q)
%P and Q uncertainties fixed or read? 
%Fixed: set up here and passed to the sampling function (below as fixed for
%testing)
%Read: add to the node variable a matrix with the uncertainties, keep in
%mind relative values
%Do somenthing and only pass a matrix "X" with two columns, mean and std, and num_samples 
X=cell2mat(node.type(:,2));
X=[X abs((1/3000)*X)];
raw_samples_P=sampling(X,num_samples,2);
raw_samples_P(isnan(raw_samples_P))=0;
X=cell2mat(node.type(:,3));
X=[X abs((1/300)*X)];
raw_samples_Q=sampling(X,num_samples,2);
raw_samples_Q(isnan(raw_samples_Q))=0;

%Convert matrix with samples to struct "node_samples" to be used in the
%main code
node_samples=node;
node_samples.type=cell(size(node.type,1),3,num_samples);
P=reshape(raw_samples_P',[size(node.type,1),1,num_samples]);
Q=reshape(raw_samples_Q',[size(node.type,1),1,num_samples]);
node_samples.type(:,2,:)=num2cell(P);
node_samples.type(:,3,:)=num2cell(Q);

%% Measurements sampling (V,Pi,Qi,Pij,Qij,Iij(?))
zdata_samples=sampling([zdata(:,2) zdata(:,6)],num_samples,2)';






end