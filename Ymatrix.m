function [Y, Bs] = Ymatrix(branch, node)

Y = zeros(node.num);
Bs = zeros(node.num);
R = branch.R;
X = branch.X;
B = branch.B;

for L = 1:branch.num              %%% computation of the terms outside the main diagonal
    i = branch.start(L);
    j = branch.end(L);
    a = R(L)/(R(L)^2+X(L)^2);
    b = -X(L)/(R(L)^2+X(L)^2);
    c = complex(a,b);
    Y(i,j) = -c;
    Y(j,i) = Y(i,j);
	Y(i,i) = Y(i,i) + B(L);
	Y(j,j) = Y(j,j) + B(L);
	Bs(i,j) = B(L);
	Bs(j,i) = B(L);
end

for m = 1:node.num
    for n = 1:node.num             %%% computation of the terms in the main diagonal
        if m ~= n
            Y(m,m) = Y(m,m) - Y(m,n);
        end
    end
end
end
