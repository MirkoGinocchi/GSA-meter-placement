function [zdata] = Meas_data_generation(zdatatrue, errn)

ztrue = zdatatrue(:,2);
dev_std = zdatatrue(:,6);

z = ztrue + dev_std.*errn;   % addition of random noise to the true value according to uncertainty distribution
zdata = zdatatrue;    % only one change with respect to zdatatrue -> the measurement values are changed; true values are replaced by corrupted values to consider the uncertainty of the measurements
zdata(:,2) = z;


end

