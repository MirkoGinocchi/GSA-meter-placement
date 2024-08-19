function branch2 = network_data_corruption(branch, net)

std_devR = (net.unc/300)*branch.R;
std_devX = (net.unc/300)*branch.X;
branch2 = branch;

if net.type == 0
    branch2.R = branch.R + std_devR.*randn(branch.num,1);
    branch2.X = branch.X + std_devX.*randn(branch.num,1);
elseif net.type == 1
    err = randn(1);
    branch2.R = branch.R + std_devR*err;
    branch2.X = branch.X + std_devX*err;
%     branch2.R = (1-net.unc/100)*branch.R;
%     branch2.X = (1-net.unc/100)*branch.X;
end