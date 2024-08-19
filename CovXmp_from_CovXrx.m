function covXmp = CovXmp_from_CovXrx(covXrx, X)
% This function computes the covariance matrix of magnitude and phase
%  angles of a quantity from its real and imaginary parts
%
% SYNTAX:
%   covXmp = CovXmp_from_CovXrx(X, covXrx)
%
% INPUTS:
%   covXrx : covariance matrix of the real and imaginary components of X
%   X : structure of the quantity whose uncertainty has to be computed,
%        which includes real, imaginary, magnitude, phase angles and  
%        complex vectors of X
%
% OUTPUTS:
%   covXmp : covariance matrix of the magnitude and phase angle components 
%             of X.

%% Covariance matrix calculation (uncertainty propagation)
a = diag(cos(X.phase));                         % derivative of magnitude with respect to the real component
b = diag(sin(X.phase));                         % derivative of magnitude with respect to the imaginary component
c = diag(sin(X.phase)./abs(X.mag));             % derivative of phase angle with respect to the real component
d = diag(cos(X.phase)./abs(X.mag));             % derivative of phase angle with respect to the imaginary component
der = [a, b; c, d];                         % overall derivative matrix
covXmp = der * covXrx * der';          % output covariance matrix of real and imaginary phasors (uncertainty propagation law)