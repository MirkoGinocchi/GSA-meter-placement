function[V, I, S, num_iter,uncVm_percent] = NV_DSSE_traditional(branch, node, V, zdata, Y, Bs)

type = zdata(:,1);
z = zdata(:,2);
fbus = zdata(:,4);
tbus = zdata(:,5);
dev_std = zdata(:,6);
sigma2 = dev_std.^2;                                                       % Vector of meas covariances
sigma2(sigma2<10^-10) = 10^-10;                                            % Covariance limit (otherwise ill-conditioning problems)
vii = (type == 1);
pii = (type == 2);
qii = (type == 3);
pfi = (type == 4);
qfi = (type == 5);
ifi = (type == 6);

vi = find(vii);                                                            % Index of voltage magnitude measurements..
pi = find(pii);                                                            % Index of real power injection measurements..
qi = find(qii);                                                            % Index of reactive power injection measurements..
pf = find(pfi);                                                            % Index of real powerflow measurements..
qf = find(qfi);                                                            % Index of reactive powerflow measurements..
iampf = find(ifi);                                                         % Index of current amplitude measurements..

nvi = length(vi);                                                          % Number of Voltage amplitude measurements..
busvi = fbus(vi);                                                          % Nodes of Voltage amplitude measurements..
npi = length(pi);                                                          % Number of Real Power Injection measurements..
busppi = fbus(pi);                                                         % Nodes of Real Power Injection measurements..
npf = length(pf);                                                          % Number of Real Power Flow measurements..
fbuspf = fbus(pf);                                                         % Start Nodes of Real Power Flow measurements..
tbuspf = tbus(pf);                                                         % End Nodes of Real Power Flow Measurements..
niampf = length(iampf);                                                    % Number of Current measurements..
fbusiampf = fbus(iampf);                                                   % Start Nodes of Current amplitude measurements..
tbusiampf = tbus(iampf);                                                   % End Nodes of Current amplitude measurements..

G = real(Y);                                                               % parte reale matrice di ammettenza
B = imag(Y);                                                               % parte immaginaria matrice di ammettenza
Yabs = abs(Y);
Yphase = angle(Y);

P_inj = zdata(pii,2);
Q_inj = zdata(qii,2);
P_br = zdata(pfi,2);
Q_br = zdata(qfi,2);

Cov = diag(sigma2);
offset = node.num - 1;                                                     % offset used to place in the Jacobian the derivatives with respect to the imaginary current

%%% Definition constant Jacobian sub-matrices

% 2-3) Power Injection Measurements
H2 = zeros(npi,2*node.num-1);
H3 = zeros(npi,2*node.num-1);
for i=1:npi
    m = busppi(i);
    H2(i,1:node.num) = G(m,:);
    H2(i,node.num+1:end) = -B(m,2:end);
    H3(i,1:node.num) = B(m,:);
    H3(i,node.num+1:end) = G(m,2:end);
end

% 4-5) Power Flow Measurements
H4 = zeros(npf,2*node.num-1);
H5 = zeros(npf,2*node.num-1);
for i=1:npf
    m = fbuspf(i);
    n = tbuspf(i);
    H4(i,m) = -G(m,n);
    H4(i,n) = G(m,n);
    H5(i,m) = -B(m,n);
    H5(i,n) = B(m,n);
    H4(i,n + offset) = -B(m,n);
    H5(i,n + offset) = G(m,n);
    if (m > 1)
        H4(i,m + offset) = B(m,n);
        H5(i,m + offset) = -G(m,n);
    end
end
    
W = Cov\eye(size(Cov));                                                    % Weighting matrix
W = sparse(W);
State = [V.real; V.imag(2:end)];                                                   % State vector

num_iter = 0;                                                              % initialization number of iterations
epsilon = 5;                                                               % initialization variable used to check the convergence of the SE algorithm

while (epsilon > 10^-6 && num_iter <= 100)
    
    %%% Computation equivalent current measurements for pwr inj and branch pwr;
    Ir_inj = (P_inj.*V.real(busppi)+Q_inj.*V.imag(busppi))./(V.mag(busppi).^2);
    Ix_inj = (P_inj.*V.imag(busppi)-Q_inj.*V.real(busppi))./(V.mag(busppi).^2);
    z(pi,1) = Ir_inj;
    z(qi,1) = Ix_inj;                                              
    
    Ir_br = (P_br.*V.real(fbuspf)+Q_br.*V.imag(fbuspf))./(V.mag(fbuspf).^2);
    Ix_br = (P_br.*V.imag(fbuspf)-Q_br.*V.real(fbuspf))./(V.mag(fbuspf).^2);
    z(pf,1) = Ir_br;
    z(qf,1) = Ix_br;
        
    %%% 1) Voltage Magnitude Measurements
    h1 = V.mag(busvi); 
    H1 = zeros(nvi,2*node.num-1);
    for i=1:nvi
        m = busvi(i);
        H1(i,m) = cos(V.phase(m));
        if busvi(i) > 1 
            H1(i,m + offset) = sin(V.phase(m));
        end
    end
    
    %%% 2-3) Power Injection Measurements
    h2 = H2*State;
    h3 = H3*State;
    
    %%% 4-5) Power Flow Measurements
    h4 = H4*State;          
    h5 = H5*State; 
        
    %%% 6) Current Magnitude Measurements
    if num_iter == 0
        h6 = [];
        H6 = [];
        z2 = z;
        z(iampf,:) = [];
        W2 = W;
        W(iampf,:) = [];
        W(:,iampf) = [];
    else
        h6_re = zeros(niampf,1);
        h6_im = zeros(niampf,1);
        h6 = zeros(niampf,1);
        H6 = zeros(niampf,2*node.num-1);
        for i = 1:niampf
            m = fbusiampf(i);
            n = tbusiampf(i);
            h6_re(i,1) = Yabs(m,n)*((V.real(n)-V.real(m))*cos(Yphase(m,n)) +  (V.imag(m)-V.imag(n))*sin(Yphase(m,n)));
            h6_im(i,1) = Yabs(m,n)*((V.real(n)-V.real(m))*sin(Yphase(m,n)) +  (V.imag(n)-V.imag(m))*cos(Yphase(m,n)));
            h6(i,1) = sqrt(h6_re(i,1)^2 + h6_im(i,1)^2);
            H6(i,m) = - Yabs(m,n)*(cos(Yphase(m,n))*h6_re(i,1)+sin(Yphase(m,n))*h6_im(i,1))/h6(i,1);
            H6(i,n) = - H6(i,m);
            H6(i,n + offset) = Yabs(m,n)*(cos(Yphase(m,n))*h6_im(i,1)-sin(Yphase(m,n))*h6_re(i,1))/h6(i,1);
            if (m > 1)
                H6(i,m + offset) = - H6(i,n + offset);
            end
        end
    end
    if num_iter == 1 && niampf>0
        z = z2;
        W = W2;
    end
        
    %%%  fine misure
    
    H = [H1; H2; H3; H4; H5; H6];                                     % creation Jacobian matrix
    H = sparse(H);
    h = [h1; h2; h3; h4; h5; h6];                                     % vector of all the meas functions
    r = z-h;                                                               % computation meas residuals
    g = H'*W*r;                                                            % right term normal equations
    Gm = H'*W*H;                                                           % Gain matrix
    
    Delta_X = Gm\g;                                                        % solution normal equations
    State = State + Delta_X;                                               % update state variables
    epsilon = max(abs(Delta_X));                                           % computation convergence threshold
    
    V.real = State(1:node.num);
    V.imag(2:end) = State(node.num+1:end);
    V.complex = complex(V.real,V.imag);
    V.mag = abs(V.complex);
    V.phase = angle(V.complex);
    
    num_iter = num_iter + 1;

end

% K=(H/Gm)*H'*W;
% Ss=eye(size(K,1))-K;
% omega=Ss/W;
cov = inv(Gm);
covVrx = [cov(1:node.num,1:node.num),zeros(node.num,1),cov(1:node.num,node.num+1:end); zeros(1,2*node.num); cov(node.num+1:end,1:node.num),zeros(node.num-1,1),cov(node.num+1:end,node.num+1:end)];
covVmp = CovXmp_from_CovXrx(covVrx, V);
varVm = diag(covVmp(1:node.num,1:node.num));    %vector of variances of the voltage magnitudes
stddevVm = sqrt(varVm);                         %vector of standard deviation of the voltage magnitudes
uncVm_percent = 300*stddevVm./V.mag;            %vector of percent uncertainties (with coverage factor 3)
for i = 1:branch.num
    from = branch.start(i);
    to = branch.end(i);
    Irx1(i,1) = - (V.complex(from) - V.complex(to))*Y(from,to) + 1i*Bs(from,to)*V.complex(from);
    Irx2(i,1) = (V.complex(from) - V.complex(to))*Y(from,to) + 1i*Bs(to,from)*V.complex(to);
end
Ir1 = real(Irx1);
Ix1 = imag(Irx1);

Ir2 = real(Irx2);
Ix2 = imag(Irx2);

I.br1.complex = Irx1;
I.br1.mag = abs(Irx1);
I.br1.phase = angle(Irx1);
I.br1.real = Ir1;
I.br1.imag = Ix1;

I.br2.complex = Irx2;
I.br2.mag = abs(Irx2);
I.br2.phase = angle(Irx2);
I.br2.real = Ir2;
I.br2.imag = Ix2;

for i = 1:node.num
    to = branch.end == i;
    from = branch.start == i;
    I.inj.real(i,1) = sum(Ir2(to)) + sum(Ir1(from));
    I.inj.imag(i,1) = sum(Ix2(to)) + sum(Ix1(from));
end

I.inj.complex = complex(I.inj.real, I.inj.imag);
I.inj.mag = abs(I.inj.complex);
I.inj.phase = angle(I.inj.complex);
S.inj.complex = V.complex .* conj(I.inj.complex);
S.inj.real = real(S.inj.complex);
S.inj.imag = imag(S.inj.complex);
S.br1.complex = V.complex(branch.start) .* conj(I.br1.complex);
S.br1.real = real(S.br1.complex);
S.br1.imag = imag(S.br1.complex);
S.br2.complex = V.complex(branch.end) .* conj(-I.br1.complex);
S.br2.real = real(S.br2.complex);
S.br2.imag = imag(S.br2.complex);
