%PAMODEL - Parameter fitting for PA gain compression model (H1 line)
%Fits parameters alpha, beta, gamma, delta in f(x) = 10*log10(alpha/(sqrt(beta + exp(2*(x-gamma)/delta))))
%using observed data from datasheet for different PA cases (ZX60-V62+, ZX60-V63+)
%
% Syntax:  run the script
%
% Inputs:
%    pa_case - 'V62+' or 'V63+' (string), selects PA model parameters
%
% Outputs:
%    alpha, beta, gamma, delta      - Fitted model parameters
%    alpha_hat, beta_hat, gamma_hat, delta_hat - Rounded fitted parameters
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: lsqcurvefit

% Project: Parameter fitting for PA gain compression model
%
% Thanks:
%    Denis Gilbert (2004). M-file Header Template
%    (https://www.mathworks.com/matlabcentral/fileexchange/4908-m-file-header-template),
%    MATLAB Central File Exchange. 
% 
% Author: Germain PHAM
% Affiliation: C2S, Télécom Paris, IP Paris
% email: 
% September 2025; Last revision: 

%------------- BEGIN CODE --------------

% Solve for alpha, beta, gamma, delta in f(x) = 10*log10(alpha/(sqrt(beta + exp(2*(x-gamma)/delta))))

% Example context: ZX60-V62+
% Small signal gain = 15.0 dB
% Output Power @ 1 dB compression = 19.0 dBm => Input = 19-14.0 = 5 dBm

% Example context: ZX60-V63+
% Small signal gain = 20 dB
% Output Power @ 1 dB compression = 17.8 dBm => Input = 17.8-19.0 = -1.2 dBm

pa_case = 'V63+';

if contains(pa_case, 'V62+')
    disp('Using parameters for ZX60-V62+');
    gain = 15; % Small signal gain
    p1db_out = 19; % Output Power @ 1 dB compression
elseif contains(pa_case, 'V63+')
    disp('Using parameters for ZX60-V63+');
    gain = 20; % Small signal gain
    p1db_out = 17.8; % Output Power @ 1 dB compression
else
    error('Unknown PA case. Please specify either V62+ or V63+.');
end

% Observed data ("manually" infered from datasheet)
p1db_in = p1db_out - (gain-1); % Input Power @ 1 dB compression
xdata = [p1db_in-35, p1db_in-6, p1db_in, p1db_in+1];
ydata = [gain, gain, gain-1, gain-2];

% Model function
model = @(params, x) 10*log10(params(1) ./ (sqrt(params(2) + exp(2*(x - params(3)) / params(4)))));

% Initial guess for [alpha, beta, gamma, delta]
init_params = [20, 1, 2, 2];

% Use lsqcurvefit to fit parameters
opts = optimset('Display','off');
params_fit = lsqcurvefit(@(p,x) model(p,x), init_params, xdata, ydata, [], [], opts);

% Display results
alpha = params_fit(1);
beta = params_fit(2);
gamma = params_fit(3);
delta = params_fit(4);

fprintf('alpha = %.6f\n', alpha);
fprintf('beta  = %.6f\n', beta);
fprintf('gamma = %.6f\n', gamma);
fprintf('delta = %.6f\n', delta);

% rounded fitted parameters
[alpha_hat, beta_hat, gamma_hat, delta_hat] = deal(round(alpha,2), round(beta,2), round(gamma,2), round(delta,2));
fprintf('Rounded parameters:\n');
fprintf('alpha_hat = %.2f\n', alpha_hat);
fprintf('beta_hat  = %.2f\n', beta_hat);
fprintf('gamma_hat = %.2f\n', gamma_hat);
fprintf('delta_hat = %.2f\n', delta_hat);

% Plot results
x_fit = linspace(-40, 20, 100);
y_fit = model(params_fit, x_fit);
y_fit_hat = model([alpha_hat, beta_hat, gamma_hat, delta_hat], x_fit);
% figure;
plot(xdata, ydata, 'ro', 'MarkerFaceColor', 'r'); hold on;
plot(x_fit, y_fit, 'b-');
plot(x_fit, y_fit_hat, 'g--');

xlabel('x');
ylabel('f(x)');
title('Parameter Fitting for f(x)');
legend('Data', 'Fitted Curve', 'Rounded Fitted Curve','Location','Best');
grid on;



