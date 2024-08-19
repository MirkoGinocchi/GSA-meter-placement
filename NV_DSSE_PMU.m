function[V, I, S, num_iter] = NV_DSSE_PMU(branch, node, V, zdata, Y, Bs)

type = zdata(:,1);
z = zdata(:,2);
fbus = zdata(:,4);
tbus = zdata(:,5);
dev_std = zdata(:,6);
sigma2 = dev_std.^2;
sigma2(sigma2<10^-10) = 10^-10;
vii = (type == 1);
pii = (type == 2);
qii = (type == 3);
pfi = (type == 4);
qfi = (type == 5);
ifi = (type == 6);  
syncvi = (type == 7);  
syncii = (type == 9);  

vi = find(vii);                                                            % Index of voltage magnitude measurements..
pi = find(pii);                                                            % Index of real power injection measurements..
qi = find(qii);                                                            % Index of reactive power injection measurements..
pf = find(pfi);                                                            % Index of real powerflow measurements..
qf = find(qfi);                                                            % Index of reactive powerflow measurements..
iampf = find(ifi);                                                         % Index of current amplitude measurements..
syncv = find(syncvi);                                                      % Index of synchro voltage amplitude measurements..
synci = find(syncii);                                                      % Index of synchro current amplitude measurements..

nvi = length(vi);                                                          % Number of Voltage amplitude measurements..
npi = length(pi);                                                          % Number of Real Power Injection measurements..
busppi = fbus(pi);                                                         % Nodes of Real Power Injection measurements..
npf = length(pf);                                                          % Number of Real Power Flow measurements..
fbuspf = fbus(pf);                                                         % Start Nodes of Real Power Flow measurements..
tbuspf = tbus(pf);                                                         % End Nodes of Real Power Flow measurements..
niampf = length(iampf);                                                    % Number of Current measurements..
nsyncvamp = length(syncv);                                                 % Number of Synchro Voltage amplitude measurements.. 
nsyncvangle = length(syncv);                                               % Number of Synchro Voltage angle measurements.. 
nsyncv = nsyncvamp + nsyncvangle;
bussyncvamp = fbus(syncv);                                                 % Nodes of Synchro Voltage amplitude measurements..    
nsynciamp = length(synci);                                                 % Number of Synchro Current amplitude measurements.. 
fbussynciamp = fbus(synci);                                                % Start Nodes of Synchro Current amplitude measurements..
tbussynciamp = tbus(synci);                                                % End Nodes of Synchro Current amplitude measurements..

pistart = nvi + 1;                                                         % Index of start of Real Power Injection measurements..
qistart = pistart + npi;                                                   % Index of start of Reactive Power Injection measurements..
pfstart = qistart + npi;                                                   % Index of start of Real Power Flow measurements..
qfstart = pfstart + npf;                                                   % Index of start of Reactive Power Flow measurements..
iampfstart = qfstart + npf;                                                % Index of start of Current measurements..
syncvstart = iampfstart + niampf;                                          % Index of start of Synchro Voltage measurements..
syncistart = syncvstart + nsyncv;                                          % Index of start of Synchro Current measurements.. 

G = real(Y);                                                               % real part admittance matrix
B = imag(Y);                                                               % imaginary part admittance matrix

P_inj = zdata(pii,2);
Q_inj = zdata(qii,2);
P_br = zdata(pfi,2);
Q_br = zdata(qfi,2);

hidx = 1 : length(sigma2);                                                 % column index for the definition of the sparse covariance matrix
vidx = 1 : length(sigma2);                                                 % row index for the definition of the sparse covariance matrix
covariances = sigma2;
offset = node.num;

%%% Definition of the constant Jacobian matrices

%%% 2-3) Power Injection Measurements
H2 = zeros(npi,2*node.num);
H3 = zeros(npi,2*node.num);
for i = 1:npi
    m = busppi(i);
    H2(i,1:node.num) = G(m,:);
    H2(i,node.num+1:end) = -B(m,:);
    H3(i,1:node.num) = B(m,:);
    H3(i,node.num+1:end) = G(m,:);
end

%%% 4-5) Power Flow Measurements
H4 = zeros(npf,2*node.num);
H5 = zeros(npf,2*node.num);
for i=1:npf
    m = fbuspf(i);
    n = tbuspf(i);
    H4(i,m) = -G(m,n);
    H4(i,n) = G(m,n);
    H4(i,m + offset) = B(m,n);
    H4(i,n + offset) = -B(m,n);
    H5(i,m) = -B(m,n);
    H5(i,n) = B(m,n);
    H5(i,m + offset) = -G(m,n);
    H5(i,n + offset) = G(m,n);
end

%%% 7-8) Voltage PMU Measurements
H78 = zeros(2*nsyncvamp,2*node.num);
for i=1:nsyncvamp
    idx = syncvstart + 2*(i-1);
    vamp = z(idx);
    vangle = z(idx+1);
    z(idx) = vamp * cos(vangle);
    z(idx + 1) = vamp * sin(vangle);
    rot_mat = [cos(vangle) , -vamp*sin(vangle); sin(vangle), vamp*cos(vangle)];
    covar_syncv_rect = rot_mat * [sigma2(idx), 0; 0, sigma2(idx+1)] * rot_mat';
    hidx = [hidx, idx, idx + 1];
    vidx = [vidx, idx+1, idx];
    covariances(idx) = covar_syncv_rect(1,1);
    covariances(idx+1) = covar_syncv_rect(2,2);
    covariances = [covariances; covar_syncv_rect(1,2); covar_syncv_rect(2,1)];
    index = 2*i-1;
    stateidx = bussyncvamp(i);
    H78(index, stateidx) = 1;
    H78(index+1, stateidx + offset) = 1;
end

%%% 9-10) Current PMU Measurements
H910 = zeros(2*nsynciamp,2*node.num);
for i=1:nsynciamp                
    idx = syncistart + 2*(i-1);           
    iamp = z(idx);
    iangle = z(idx+1);
    z(idx) = iamp * cos(iangle);
    z(idx + 1) = iamp * sin(iangle);
    rot_mat = [cos(iangle) , -iamp*sin(iangle); sin(iangle), iamp*cos(iangle)];
    covar_synci_rect = rot_mat * [sigma2(idx), 0; 0, sigma2(idx+1)] * rot_mat';
    hidx = [hidx, idx, idx + 1];
    vidx = [vidx, idx+1, idx];
    covariances(idx) = covar_synci_rect(1,1);
    covariances(idx+1) = covar_synci_rect(2,2);
    covariances = [covariances; covar_synci_rect(1,2); covar_synci_rect(2,1)];
    m = fbussynciamp(i);
    n = tbussynciamp(i);
    idx = 2*i-1;  
    H910(idx,m) = - G(m,n);
    H910(idx,n) = G(m,n);
    H910(idx+1,m) = - B(m,n);
    H910(idx+1,n) = B(m,n);
    H910(idx,m + offset) = B(m,n);
    H910(idx,n + offset) = - B(m,n);
    H910(idx+1,m + offset) = - G(m,n);
    H910(idx+1,n + offset) = G(m,n);
end

H = [H2; H3; H4; H5; H78; H910];
H = sparse(H);
Cov = sparse(hidx, vidx, covariances);
Cov = full(Cov);
W = Cov\eye(size(Cov));
W = sparse(W);
State = [V.real; V.imag];
Gm = H'*W*H;

num_iter = 0;
epsilon = 5;

while (epsilon > 10^-6 && num_iter <= 100)
    
    Ir_inj = (P_inj.*V.real(busppi)+Q_inj.*V.imag(busppi))./(V.mag(busppi).^2);
    Ix_inj = (P_inj.*V.imag(busppi)-Q_inj.*V.real(busppi))./(V.mag(busppi).^2);
    z(pi,1) = Ir_inj;
    z(qi,1) = Ix_inj;
    
    Ir_br = (P_br.*V.real(fbuspf)+Q_br.*V.imag(fbuspf))./(V.mag(fbuspf).^2);
    Ix_br = (P_br.*V.imag(fbuspf)-Q_br.*V.real(fbuspf))./(V.mag(fbuspf).^2);
    z(pf,1) = Ir_br;
    z(qf,1) = Ix_br;
        
    h = H*State;
    r = z-h;
    g = H'*W*r;
    
    Delta_X = Gm\g;
    State = State + Delta_X;              % update state variables
    epsilon = max(abs(Delta_X));
    
    V.real = State(1:node.num);
    V.imag = State(node.num+1:end);
    V.complex = complex(V.real,V.imag);
    V.mag = abs(V.complex);
    V.phase = angle(V.complex);
    
    num_iter = num_iter + 1;
end

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