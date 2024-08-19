function [samples]=sampling(X,N,pdf)

%Function to sample normally distributed variables
%X: Matrix with two columns, mean and std. number of rows=number of
%variables
%N: Number of samples
%Samples: matrix with the sampled data Nxd (# of samples x number of
%variables)
%pdf: Uniform (1) or normal (2) distribution

%% Sobol sequence to sample within U(0,1)
d=size(X,1);
p=sobolset(d);
p=scramble(p,'MatousekAffineOwen');
samplesU=p(1:N,:);
% figure(),plot(samplesU(:,1),samplesU(:,2),'*')

%% Sampling in the pdf
samples=zeros(N,d);
switch pdf
    case 1 %Uniform
        for i=1:d
            samples(:,i)=unifinv(samplesU(:,i),X(i,1)*(1-X(i,2)),X(i,1)*(1+X(i,2)));
        end
    case 2 %Normal
        for i=1:d
            samples(:,i)=norminv(samplesU(:,i),X(i,1),X(i,2));
        end
    case 3 %Discrete Uniform
        for i=1:d
            samples(:,i)=unidinv(samplesU(:,i),2)-1;
        end
end

end