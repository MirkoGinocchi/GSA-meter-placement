function [V, I, S, num_iter] = Power_flow_NV(branch, node, Y, Bs) 
  
State = [ones(node.num,1); zeros(node.num,1)];

G = real(Y);              
B = imag(Y);

z = zeros(2*(node.num),1);
h = zeros(2*(node.num),1);
H = zeros(2*(node.num));

for i = 1:node.num
    j = 2*(i-1)+1;
    type = node.type{i,1};
    switch type
        case 'slack'
            z(j) = node.type{i,2}*cos(node.type{i,3});
            z(j+1) = node.type{i,2}*sin(node.type{i,3});
            H(j,i) = 1;
            H(j+1,i+node.num) = 1;
        case 'PQ'
            H(j, 1:node.num) = G(i,:);
            H(j, node.num+1:end) = -B(i,:);
            H(j+1, 1:node.num) = B(i,:);
            H(j+1, node.num+1:end) = G(i,:);
        case 'PV'
            z(j) = node.type{i,2};
            z(j+1) = node.type{i,3};
    end
end

epsilon = 5;
Vabs = ones(node.num,1);
Vtheta = zeros(node.num,1);
Vr = Vabs;
Vx = Vtheta;
num_iter = 0;

while epsilon > 10^-7 && num_iter<200
    for i = 1:node.num
        j = 2*(i-1)+1;
        type = node.type{i,1};
        switch type
            case 'slack'
                h(j) = H(j,:)*State;
                h(j+1) = H(j+1,:)*State;
            case 'PQ'
                z(j) = (node.type{i,2}*Vr(i)+node.type{i,3}*Vx(i))./(Vabs(i)^2);
                z(j+1) = (node.type{i,2}*Vx(i)-node.type{i,3}*Vr(i))./(Vabs(i)^2);
                h(j) = H(j,:)*State;
                h(j+1) = H(j+1,:)*State;
            case 'PV'
                h(j) = Vr(i)*(G(i,:)*Vr-B(i,:)*Vx) + Vx(i)*(G(i,:)*Vx + B(i,:)*Vr);
                h(j+1) = Vabs(i);
                H(j, 1:node.num) = Vr(i).*G(i,:) + Vx(i).*B(i,:);
                H(j, i) = H(j,i) + G(i,:)*Vr - B(i,:)*Vx;
                H(j, node.num+1: end) = Vx(i).*G(i,:) - Vr(i).*B(i,:);
                H(j, i + node.num) = H(j,i+node.num) + G(i,:)*Vx + B(i,:)*Vr;
                H(j+1, i) = cos(Vtheta(i));
                H(j+1, i+node.num) = sin(Vtheta(i));
        end
    end
    
    r = z-h;
%     Hinv = inv(H);
    Delta_State = H\r;
    State = State + Delta_State;
    epsilon = max(abs(Delta_State));
    
    Vr = State(1:node.num);
    Vx = State(node.num+1:end);
    
    Vcomplex = complex(Vr,Vx);
    Vabs = abs(Vcomplex);
    Vtheta = angle(Vcomplex);
    
    num_iter = num_iter + 1;
end
V.complex.true_val = Vcomplex;
V.mag.true_val = Vabs;
V.phase.true_val = Vtheta;
V.real.true_val = Vr;
V.imag.true_val = Vx;

for i = 1:branch.num
    from = branch.start(i);
    to = branch.end(i);
    Irx1(i,1) = - (Vcomplex(from) - Vcomplex(to))*Y(from,to) + 1i*Bs(from,to)*Vcomplex(from);
    Irx2(i,1) = (Vcomplex(from) - Vcomplex(to))*Y(from,to) + 1i*Bs(to,from)*Vcomplex(to);
end
Ir1 = real(Irx1);
Ix1 = imag(Irx1);

Ir2 = real(Irx2);
Ix2 = imag(Irx2);

I.br1.complex.true_val = Irx1;
I.br1.mag.true_val = abs(Irx1);
I.br1.phase.true_val = angle(Irx1);
I.br1.real.true_val = Ir1;
I.br1.imag.true_val = Ix1;

I.br2.complex.true_val = Irx2;
I.br2.mag.true_val = abs(Irx2);
I.br2.phase.true_val = angle(Irx2);
I.br2.real.true_val = Ir2;
I.br2.imag.true_val = Ix2;

for i = 1:node.num
    to = branch.end == i;
    from = branch.start == i;
    I.inj.real.true_val(i,1) = - sum(Ir1(to)) + sum(Ir1(from));
    I.inj.imag.true_val(i,1) = - sum(Ix1(to)) + sum(Ix1(from));
end

I.inj.complex.true_val = complex(I.inj.real.true_val, I.inj.imag.true_val);
I.inj.mag.true_val = abs(I.inj.complex.true_val);
I.inj.phase.true_val = angle(I.inj.complex.true_val);
S.inj.complex.true_val = V.complex.true_val .* conj(I.inj.complex.true_val);
S.inj.real.true_val = real(S.inj.complex.true_val);
S.inj.imag.true_val = imag(S.inj.complex.true_val);
S.br1.complex.true_val = V.complex.true_val(branch.start) .* conj(I.br1.complex.true_val);
S.br1.real.true_val = real(S.br1.complex.true_val);
S.br1.imag.true_val = imag(S.br1.complex.true_val);
S.br2.complex.true_val = V.complex.true_val(branch.end) .* conj(-I.br1.complex.true_val);
S.br2.real.true_val = real(S.br2.complex.true_val);
S.br2.imag.true_val = imag(S.br2.complex.true_val);

