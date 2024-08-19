function [V, I, S, num_iter, uncVm_percent] = NV_DSSE_call(branch, node, zdata, Y, Bs, V, I, S)

type = zdata(:,1);                                                         % type of measurements

Vest.mag = ones(node.num,1);                                               % initialization of voltage magnitude and phase angle vectors;
Vest.phase = zeros(node.num,1);
Vest.real = Vest.mag;
Vest.imag = Vest.phase;

trad_code = 0;                                                             % code used to identify if there are conventional measurements in the input set (other real measurements other than PMUs)
PMU_code = 0;                                                              % code used to identify if there are PMU measurements in the input set

Vmag_measi = (type == 1);
Vmag_meas = find(Vmag_measi);                                              % identification of voltage magnitude meas
if length(Vmag_meas) >= 1
    trad_code = 1;
end
Vsync_measi = (type == 7);
Vsync_meas = find(Vsync_measi);                                            % identification of PMU voltage meas
if length(Vsync_meas) >= 1
    PMU_code = 2;
end

estimator_code = trad_code + PMU_code;                                     % code used to identify if: 1) there are only conventional meas; 2) there are only PMU meas; 3) there are both conventional and PMU meas.
                                           
switch estimator_code
    case 1
        [Vest, Iest, Sest, num_iter,uncVm_percent] = NV_DSSE_traditional(branch, node, Vest, zdata, Y, Bs);    % estimator used if only conventional measurements are present
%         K=full(K);
    case 2
        [Vest, Iest, Sest, num_iter] = NV_DSSE_PMU(branch, node, Vest, zdata, Y, Bs);            % estimator used if only PMU measurements are present
    case 3
        [Vest, Iest, Sest, num_iter] = NV_DSSE_hybrid(branch, node, Vest, zdata, Y, Bs);         % estimator used if both conventional and PMU measurements are present
end

V.complex.est = Vest.complex;
V.mag.est = Vest.mag;
V.phase.est = Vest.phase;
V.real.est = Vest.real;
V.imag.est = Vest.imag;

I.br1.complex.est = Iest.br1.complex;
I.br1.mag.est = Iest.br1.mag;
I.br1.phase.est = Iest.br1.phase;
I.br1.real.est = Iest.br1.real;
I.br1.imag.est = Iest.br1.imag;

I.inj.complex.est = Iest.inj.complex;
I.inj.mag.est = Iest.inj.mag;
I.inj.phase.est = Iest.inj.phase;
I.inj.real.est = Iest.inj.real;
I.inj.imag.est = Iest.inj.imag;

S.br1.complex.est = Sest.br1.complex;
S.br1.real.est = Sest.br1.real;
S.br1.imag.est = Sest.br1.imag;

S.br2.complex.est = Sest.br2.complex;
S.br2.real.est = Sest.br2.real;
S.br2.imag.est = Sest.br2.imag;

S.inj.complex.est = Sest.inj.complex;
S.inj.real.est = Sest.inj.real;
S.inj.imag.est = Sest.inj.imag;
