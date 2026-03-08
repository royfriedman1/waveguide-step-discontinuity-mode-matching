clc; clear; close all;

%% Physical parameters
H1 = 7.9e-3;      % Height of region 1 [m]
H2 = 10.6e-3;     % Height of region 2 [m]
y0 = -(H2 - H1)/2;     % Vertical alignment offset
c  = 3e8;         % Speed of light [m/s]

epsr1 = 2;        % Region 1 (input)
epsr2 = 1;        % Region 2 (output)

%% Cutoff frequency check for TE1
fc1 = c / (2 * H1 * sqrt(epsr1));
fc2 = c / (2 * H2 * sqrt(epsr2));
fprintf('Cutoff frequency Region 1 (TE1): %.2f GHz\n', fc1/1e9);
fprintf('Cutoff frequency Region 2 (TE1): %.2f GHz\n\n', fc2/1e9);

%% Frequency sweep setup
f_start = max(fc1, fc2) * 1.05;      % start just above the higher cutoff
f_end   = 25.00e9;
fvec    = linspace(f_start, f_end, 30);  % Frequency vector [Hz]

N       = 100;                             % Number of modes
n       = (1:N).';
S11     = zeros(size(fvec));
S21     = zeros(size(fvec));
PowerCheck = zeros(size(fvec));
Deviation = zeros(size(fvec));
epsilon = 1e-6;

%% Loop over frequency
for idx = 1:numel(fvec)
    f  = fvec(idx);
    k0 = 2*pi*f / c;

    % Compute transverse and propagation constants
    gamma_n = n*pi / H1;
    gamma_m = n*pi / H2;

    beta_n  = sqrt(k0.^2 * epsr1 - gamma_n.^2);
    beta_m  = sqrt(k0.^2 * epsr2 - gamma_m.^2);

    % Ensure only propagating modes have real beta
    is_real_n = real(beta_n) > 0;
    is_real_m = real(beta_m) > 0;
    beta_n(~is_real_n) = 1i*abs(beta_n(~is_real_n));
    beta_m(~is_real_m) = 1i*abs(beta_m(~is_real_m));

    % Normalize with power (beta/epsr)
    D1 = diag(beta_n ./ epsr1);
    D2 = diag(beta_m ./ epsr2);

    %% Overlap matrix I(m,n)
    I = zeros(N, N);
    for mIdx = 1:N
        for nIdx = 1:N
            integrand = @(y) sin(gamma_n(nIdx)*y) .* sin(gamma_m(mIdx)*(y - y0));
            I(mIdx, nIdx) = (2/H2) * integral(integrand, 0, H1, 'ArrayValued', true);
        end
    end
    M1 = I; M2 = I;

    %% Matrix formulation
    tmp = M1.' * D2;
    A   = [ -M2,              eye(N);
             D1 + epsilon*eye(N), tmp ];
    A   = A + epsilon * eye(2*N);

    %% RHS and solve
    a   = zeros(N,1); a(1) = 1;            % Only TE1 incident
    rhs = [M2 * a;
           D1 * a];
    x   = A \ rhs;

    b = x(1:N);
    c_vec = x(N+1:end);

    S11(idx) = b(1);
    S21(idx) = c_vec(1);
     % Consider only propagating modes for power
    valid_n = real(beta_n) > 0;
    valid_m = real(beta_m) > 0;
    power_in  = sum(abs(a(valid_n)).^2);  % Incident mode power (usually 1)
    power_ref = sum(abs(b(valid_n)).^2);
    power_trn = sum(abs(c_vec(valid_m)).^2);

    PowerCheck(idx) = power_ref + power_trn;
    Deviation(idx)  = abs(PowerCheck(idx) - power_in);

end

%% Display maximum deviation
fprintf('Max power conservation deviation: %.4f\n\n', max(Deviation));

%% Plot S‑parameters
figure;
plot(fvec/1e9, 20*log10(abs(S11)), 'b-', 'LineWidth', 1.5); hold on;
plot(fvec/1e9, 20*log10(abs(S21)), 'r-', 'LineWidth', 1.5);
xlabel('Frequency [GHz]');
ylabel('Magnitude [dB]');
legend('|S_{11}|','|S_{21}|','Location','Best');
title('Mode-Matching |S_{11}| & |S_{21}| vs Frequency');
grid on;

%% Plot power conservation
figure;
plot(fvec/1e9, PowerCheck, 'k', 'LineWidth', 1.5); hold on;
yline(1, '--r', 'Ideal Power'); 
xlabel('Frequency [GHz]');
ylabel('|S_{11}|^2 + |S_{21}|^2');
title('Power Conservation Check');
ylim([0.9 1.1]); grid on;

%% Plot deviation
figure;
plot(fvec/1e9, Deviation, 'm', 'LineWidth', 1.5);
xlabel('Frequency [GHz]');
ylabel('Deviation from 1');
title('Deviation from Power Conservation');
grid on;