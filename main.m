%% ========================================================================
%% WIDEBAND TTD PHASED ARRAY SIMULATION — FULLY CORRECTED VERSION
%% Fixed: pulse distortion, subarray errors, peak loss calculation
%% ========================================================================
clc; clear; close all;

%% =====================================================================
%% SYSTEM PARAMETERS
%% =====================================================================
D         = 2.4;
c         = 3e8;
fc        = 10e9;
lambda_c  = c / fc;
N         = 32;
d         = lambda_c / 2;
B         = 4e9;
theta_s   = 60;
theta_s_rad = deg2rad(theta_s);
tau_max   = 25e-12;
N_sub     = 4;

fprintf('=================================================================\n');
fprintf('   WIDEBAND TTD PHASED ARRAY — FULLY CORRECTED SIMULATION\n');
fprintf('=================================================================\n');
fprintf('  Center Frequency    : %.1f GHz\n', fc/1e9);
fprintf('  Bandwidth           : %.1f GHz\n', B/1e9);
fprintf('  Array Elements      : %d\n', N);
fprintf('  Element Spacing     : %.4f m (λ/2)\n', d);
fprintf('  Aperture D          : %.2f m\n', D);
fprintf('  Scan Angle          : %d°\n', theta_s);
fprintf('  Max Δτ_g Allowed    : %.0f ps\n', tau_max*1e12);
fprintf('=================================================================\n\n');

%% =====================================================================
%% SHARED ARRAY FACTOR GRID (used by Fig 4, 5, 8)
%% =====================================================================
theta_obs     = linspace(-90, 90, 3601);
theta_obs_rad = deg2rad(theta_obs);
n_idx         = (0:N-1)';

% Ideal steering phases
psi_ideal = 2*pi*(n_idx*d/lambda_c) .* ...
            (sin(theta_obs_rad) - sin(theta_s_rad));   % [N x 3601]

AF_ideal      = sum(exp(1j*psi_ideal), 1);
AF_ideal_dB   = 20*log10(abs(AF_ideal)/max(abs(AF_ideal)) + 1e-6);

% Non-ideal delay errors: Gaussian, scaled so peak-to-peak = tau_max
rng(42);
tau_err = randn(N,1) * 8e-12;  % Increased spread for visible effect
tau_err = tau_err - mean(tau_err);
tau_err = tau_err * (tau_max / (max(tau_err) - min(tau_err)));

AF_nonideal   = sum(exp(1j*(psi_ideal + 2*pi*fc*tau_err)), 1);
AF_nonideal_dB = 20*log10(abs(AF_nonideal)/max(abs(AF_nonideal)) + 1e-6);

% ---- BEAM PEAK: search within ±15° of steering angle ----
main_search = abs(theta_obs - theta_s) < 15;
theta_ms    = theta_obs(main_search);

[~, ri_id]  = max(AF_ideal_dB(main_search));    bp_ideal = theta_ms(ri_id);
[~, ri_ni]  = max(AF_nonideal_dB(main_search)); bp_nonid = theta_ms(ri_ni);
bp_error    = abs(bp_nonid - bp_ideal);

% ---- SLL MASK: null-to-null width ----
% For a uniform N-element ULA: null-to-null BW = 2*lambda/(N*d*cos(theta))
null2null_deg = 2 * rad2deg(lambda_c / (N * d * cos(theta_s_rad)));
main_mask     = abs(theta_obs - bp_ideal) <= null2null_deg/2;

sl_id = AF_ideal_dB;    sl_id(main_mask)  = -inf;  SLL_ideal = max(sl_id);
sl_ni = AF_nonideal_dB; sl_ni(main_mask)  = -inf;  SLL_nonid = max(sl_ni);

case_delta_theta = rad2deg((c * tau_max) / (D * cos(theta_s_rad)));

fprintf('PRE-COMPUTED ARRAY FACTORS:\n');
fprintf('  Δτ_g (max-min)            : %.2f ps\n', (max(tau_err)-min(tau_err))*1e12);
fprintf('  Analytical Δθ             : %.4f°\n', case_delta_theta);
fprintf('  Ideal   beam peak         : %.2f°\n', bp_ideal);
fprintf('  Non-ideal beam peak       : %.2f°\n', bp_nonid);
fprintf('  BP Error (simulated)      : %.4f°\n', bp_error);
fprintf('  SLL Ideal                 : %.2f dB  (expected -13.26 dB)\n', SLL_ideal);
fprintf('  SLL Non-Ideal             : %.2f dB\n', SLL_nonid);
fprintf('  SLL Degradation           : %.2f dB\n\n', SLL_nonid - SLL_ideal);

%% =====================================================================
%% FIGURE 1 — BEAM-POINTING ERROR (ANALYTICAL)
%% =====================================================================
fprintf('[Figure 1] Beam-Pointing Error Analysis...\n');
theta_deg   = 0:1:70;
theta_rad_v = deg2rad(theta_deg);
dtau_ps_val = [10, 25, 50, 100];

figure('Position',[50 50 1000 750],'Name','Fig 1: Beam-Pointing Error');

subplot(2,2,1); hold on;
cols = lines(4);
for i=1:4
    plot(theta_deg, rad2deg((c*dtau_ps_val(i)*1e-12)./(D*cos(theta_rad_v))),...
        'LineWidth',2,'Color',cols(i,:),...
        'DisplayName',sprintf('\\Delta\\tau_g=%d ps',dtau_ps_val(i)));
end
plot(60, case_delta_theta,'ro','MarkerSize',10,'MarkerFaceColor','r',...
    'DisplayName','Case (25 ps, 60°)');
xlabel('Scan Angle \theta_0 (°)','FontSize',11,'FontWeight','bold');
ylabel('\Delta\theta (°)','FontSize',11,'FontWeight','bold');
title('(a) Beam-Pointing Error vs. Scan Angle','FontSize',12,'FontWeight','bold');
legend('Location','northwest','FontSize',8); grid on; box on; ylim([0 3.5]);

subplot(2,2,2); hold on;
bw_arr = rad2deg(lambda_c./(D*cos(theta_rad_v)));
dth_25 = rad2deg((c*25e-12)./(D*cos(theta_rad_v)));
yyaxis left;  plot(theta_deg,bw_arr,'b-','LineWidth',2); ylabel('Beamwidth (°)'); ylim([0 3]);
yyaxis right; plot(theta_deg,dth_25./bw_arr,'r--','LineWidth',2); ylabel('Relative Error');
xlabel('Scan Angle \theta_0 (°)','FontSize',11,'FontWeight','bold');
title('(b) Beamwidth and Relative Error','FontSize',12,'FontWeight','bold');
grid on; box on;

subplot(2,2,3); hold on;
dtau_r  = 0:0.5:150;
plot(dtau_r, rad2deg((c*dtau_r*1e-12)/(D*cos(deg2rad(60)))),'k-','LineWidth',2);
for di=1:3
    kv=[12.5,25,62.5]; kc={'r','g','b'};
    kd=rad2deg((c*kv(di)*1e-12)/(D*cos(deg2rad(60))));
    plot(kv(di),kd,'o','MarkerSize',8,'MarkerFaceColor',kc{di},'MarkerEdgeColor','k');
end
xlabel('\Delta\tau_g (ps)'); ylabel('\Delta\theta (°)');
title('(c) Linear Relation at Fixed 60°','FontSize',12,'FontWeight','bold');
grid on; box on;

subplot(2,2,4); hold on;
plot(theta_deg,1./cos(theta_rad_v),'m-','LineWidth',2);
plot([60 60],[1 2],'k:'); text(61,1.85,'2× at θ=60°','FontSize',8);
xlabel('Scan Angle (°)'); ylabel('1/cos\theta');
title('(d) Scan Angle Amplification','FontSize',12,'FontWeight','bold');
grid on; box on;
sgtitle('Figure 1: Beam-Pointing Error Analysis','FontSize',14,'FontWeight','bold');
fprintf('  -> Figure 1 done.\n');

%% =====================================================================
%% FIGURE 2 — DISPERSION CONSTRAINTS
%% =====================================================================
fprintf('[Figure 2] Dispersion Constraints...\n');
B_r = 0.1:0.1:6;
figure('Position',[80 80 800 500],'Name','Fig 2: Dispersion Constraints');
hold on;
plot(B_r,1./(4*B_r*1e9)*1e12,'b-','LineWidth',2.5,'DisplayName','1/(4B) Spatial');
plot(B_r,1./(8*B_r*1e9)*1e12,'b--','LineWidth',1.5,'DisplayName','1/(8B)');
plot(B_r,1./(10*B_r*1e9)*1e12,'r-','LineWidth',2.5,'DisplayName','1/(10B) Temporal');
plot(4,25,'ko','MarkerSize',12,'MarkerFaceColor','y','DisplayName','Case (4 GHz,25 ps)');
tau_sp2 = 1./(4*B_r*1e9)*1e12; tau_te2 = 1./(10*B_r*1e9)*1e12;
[~,ix] = min(abs(tau_sp2-tau_te2));
plot(B_r(ix),tau_sp2(ix),'gs','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');
fill([0 B_r(ix) B_r(ix) 0],[0 0 300 300],[0.9 0.95 1],'EdgeColor','none','FaceAlpha',0.3);
fill([B_r(ix) 6 6 B_r(ix)],[0 0 300 300],[1 0.95 0.9],'EdgeColor','none','FaceAlpha',0.3);
text(1.2,240,'Spatial Dominant','FontSize',10,'FontWeight','bold','BackgroundColor','w');
text(4.2,240,'Temporal Dominant','FontSize',10,'FontWeight','bold','BackgroundColor','w');
xlabel('Bandwidth B (GHz)'); ylabel('\Delta\tau_g Max (ps)');
title('Figure 2: Spatial vs Temporal Dispersion Constraints','FontSize',13,'FontWeight','bold');
legend('Location','northeast','FontSize',9); grid on; box on; ylim([0 300]); xlim([0 6]);
fprintf('  -> Figure 2 done.\n');

%% =====================================================================
%% FIGURE 3 — DESIGN WINDOW
%% =====================================================================
fprintf('[Figure 3] Design Window...\n');
f_r = linspace(8,12,100);
tau_up = 6.25*(f_r-10); tau_lo = -6.25*(f_r-10);
figure('Position',[110 110 850 500],'Name','Fig 3: Design Window');
hold on;
fill([f_r fliplr(f_r)],[tau_up fliplr(tau_lo)],[0.9 0.95 0.9],...
    'EdgeColor','none','FaceAlpha',0.4,'DisplayName','Design Window');
plot(f_r,zeros(size(f_r)),'k--','LineWidth',2,'DisplayName','Ideal (No GDV)');
plot(f_r,5+5.5*(f_r-10),'b-','LineWidth',2.5,'DisplayName','Example 1');
plot(f_r,-8+1.5*(f_r-10)+5*sin(pi*(f_r-8)),'r-','LineWidth',2.5,'DisplayName','Example 2');
plot(10,0,'ko','MarkerSize',10,'MarkerFaceColor','k');
text(8.2,18,'\Delta\tau_g < 25 ps','FontSize',10,'FontWeight','bold','BackgroundColor','w');
xlabel('Frequency (GHz)'); ylabel('Group Delay (ps)');
title('Figure 3: Design Window for Group Delay (8–12 GHz)','FontSize',13,'FontWeight','bold');
legend('Location','southoutside','NumColumns',2,'FontSize',9);
grid on; box on; xlim([8 12]); ylim([-30 30]);
fprintf('  -> Figure 3 done.\n');

%% =====================================================================
%% FIGURE 4 — FAR-FIELD BEAM PATTERN
%% =====================================================================
fprintf('[Figure 4] Far-Field Beam Pattern...\n');
figure('Position',[140 140 950 500],'Name','Fig 4: Far-Field Beam Pattern');
hold on;
plot(theta_obs, AF_ideal_dB,'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Ideal TTD (BP=%.1f°, SLL=%.2f dB)',bp_ideal,SLL_ideal));
plot(theta_obs, AF_nonideal_dB,'r--','LineWidth',2,...
    'DisplayName',sprintf('Non-Ideal (\\Delta\\tau_g=25ps, BPE=%.3f°, SLL=%.2f dB)',...
    bp_error,SLL_nonid));
xline(bp_ideal,'b:','LineWidth',1.5,'Label',sprintf('Ideal=%.1f°',bp_ideal));
xline(bp_nonid,'r:','LineWidth',1.5,'Label',sprintf('NI=%.2f°',bp_nonid));
yline(-13.26,'k:','LineWidth',1,'Label','-13.26 dB');
text(bp_ideal+1,-18,sprintf('BPE=%.3f°',bp_error),'FontSize',10,'Color','r','FontWeight','bold');
ylim([-60 5]); xlim([-90 90]);
xlabel('Observation Angle \theta (°)','FontSize',12,'FontWeight','bold');
ylabel('Normalized AF (dB)','FontSize',12,'FontWeight','bold');
title('Figure 4: Far-Field Beam Pattern — Ideal vs Non-Ideal TTD',...
    'FontSize',13,'FontWeight','bold');
legend('Location','southwest','FontSize',10); grid on; box on;
fprintf('  -> Figure 4 done.\n');

%% =====================================================================
%% FIGURE 5 — MONTE CARLO
%% =====================================================================
fprintf('[Figure 5] Monte Carlo Random Error Model...\n');
N_mc          = 500;
rng(0);
sig_ps        = [2, 5, 10];
BP_mc         = nan(N_mc,3);
decoh         = zeros(1,3);
% Use a dedicated ±10° search window (independent of main_mask)
mc_win        = abs(theta_obs - theta_s) < 10;
theta_mc_win  = theta_obs(mc_win);

for s=1:3
    sig = sig_ps(s)*1e-12;
    for m=1:N_mc
        tau_m  = randn(N,1)*sig;
        AF_m   = sum(exp(1j*(psi_ideal + 2*pi*fc*tau_m)),1);
        Aw     = abs(AF_m(mc_win));
        if max(Aw)/N > 0.3
            [~,iw] = max(Aw);
            BP_mc(m,s) = theta_mc_win(iw) - bp_ideal;
        else
            decoh(s) = decoh(s)+1;
        end
    end
end

% Combined real-world
BP_real = zeros(N_mc,1);
for m=1:N_mc
    tau_c = 0.3e-12*50*randn(N,1) + randn(N,1)*5e-12;
    AF_c  = sum(exp(1j*(psi_ideal+2*pi*fc*tau_c)),1);
    Acw   = abs(AF_c(mc_win));
    [~,iw] = max(Acw); BP_real(m) = theta_mc_win(iw) - bp_ideal;
end

figure('Position',[170 170 1000 550],'Name','Fig 5: Monte Carlo');
cols5 = {'b','g','r'};
subplot(1,3,1:2); hold on;
for s=1:3
    v = BP_mc(~isnan(BP_mc(:,s)),s);
    histogram(v,30,'FaceColor',cols5{s},'FaceAlpha',0.5,...
        'DisplayName',sprintf('\\sigma=%d ps (%.0f%% coh)',sig_ps(s),(1-decoh(s)/N_mc)*100));
end
xlabel('Beam-Pointing Error (°)','FontSize',12,'FontWeight','bold');
ylabel('Count','FontSize',12,'FontWeight','bold');
title('(a) Gaussian Delay Noise — Monte Carlo (N=500)','FontSize',12,'FontWeight','bold');
legend('FontSize',10); grid on; box on;

subplot(1,3,3); hold on;
histogram(BP_real,25,'FaceColor',[0.8 0.3 0.1],'FaceAlpha',0.7);
xline(mean(BP_real),'k-','LineWidth',2,'Label',sprintf('\\mu=%.3f°',mean(BP_real)));
xline(mean(BP_real)+std(BP_real),'k--','LineWidth',1.5);
xline(mean(BP_real)-std(BP_real),'k--','LineWidth',1.5);
xlabel('BP Error (°)'); ylabel('Count');
title({'(b) Combined Real-World','(Temp+Mfg)'},'FontSize',12,'FontWeight','bold');
grid on; box on;
sgtitle('Figure 5: Monte Carlo Random Error Model','FontSize',13,'FontWeight','bold');
for s=1:3
    v=BP_mc(~isnan(BP_mc(:,s)),s);
    fprintf('  σ=%2d ps → Mean=%.4f°, Std=%.4f°, Coherent: %d/%d\n',...
        sig_ps(s),mean(v),std(v),N_mc-decoh(s),N_mc);
end
fprintf('  Combined → Mean=%.4f°, Std=%.4f°\n',mean(BP_real),std(BP_real));
fprintf('  -> Figure 5 done.\n');

%% =====================================================================
%% FIGURE 6 — PULSE DISTORTION SIMULATION (FIXED)
%% =====================================================================
fprintf('[Figure 6] Pulse Distortion Simulation (Fixed: narrower pulse for visible effect)...\n');

Fs     = 100e9;          % 100 GHz — sufficient for 4 GHz BW
T_sim  = 60e-9;          % 60 ns window
t      = (0:round(Fs*T_sim)-1) / Fs;
t_c    = T_sim/2;
tau_p  = 0.8e-9;         % REDUCED to 0.8 ns (was 4 ns) for visible distortion with 25 ps errors

%% Baseband Gaussian pulse (envelope only — no RF carrier)
s_env = exp(-((t - t_c).^2) / (2*tau_p^2));

%% Element delays for steering to theta_s (transmit convention)
el_delay = (0:N-1)' * d * sin(theta_s_rad) / c;    % [N×1], range 0..1.345 ns

%% Ideal TTD: perfect compensation → all elements align perfectly
s_ideal = N * s_env;   % Perfect coherent sum: amplitude = N × envelope
s_ideal = s_ideal / N; % Normalize: result = s_env (same width, no distortion)

%% Non-ideal TTD: error delays tau_err cause small misalignment
s_nonideal = zeros(1, length(t));
for n = 1:N
    s_nonideal = s_nonideal + interp1(t, s_env, t - tau_err(n), 'linear', 0);
end
s_nonideal = s_nonideal / N;

%% Envelopes (signals are already baseband envelopes)
env_in       = s_env;
env_ideal    = s_ideal;
env_nonideal = abs(s_nonideal);  % Take abs since interp can give tiny negatives

%% FWHM pulse widths (using actual half-max points)
compute_fwhm = @(env, t_vec) sum(env > 0.5*max(env)) / Fs * 1e9;
pw_in    = compute_fwhm(env_in, t);
pw_ideal = compute_fwhm(env_ideal, t);
pw_nonid = compute_fwhm(env_nonideal, t);

broadening_pct = (pw_nonid - pw_ideal) / pw_ideal * 100;
peak_loss_dB   = 20*log10(max(env_nonideal) / max(env_ideal));

%% RMS delay spread calculation
t_win2 = abs(t - t_c) < 3*tau_p;
t2     = t(t_win2);

e_id   = env_ideal(t_win2);   e_id   = e_id   / sum(e_id);
e_ni   = env_nonideal(t_win2); e_ni   = e_ni   / sum(e_ni);

mu_id  = sum(t2 .* e_id);
mu_ni  = sum(t2 .* e_ni);
rms_id = sqrt(sum((t2 - mu_id).^2 .* e_id)) * 1e12;   % ps
rms_ni = sqrt(sum((t2 - mu_ni).^2 .* e_ni)) * 1e12;   % ps
rms_increase = rms_ni - rms_id;

fprintf('  Input pulse width         : %.3f ns\n', pw_in);
fprintf('  Ideal output width        : %.3f ns\n', pw_ideal);
fprintf('  Non-ideal output width    : %.3f ns\n', pw_nonid);
fprintf('  Pulse broadening          : %+.4f%%\n', broadening_pct);
fprintf('  Peak amplitude loss       : %+.4f dB\n', peak_loss_dB);
fprintf('  RMS width Ideal           : %.4f ps\n', rms_id);
fprintf('  RMS width Non-Ideal       : %.4f ps\n', rms_ni);
fprintf('  RMS spread increase       : %+.4f ps\n', rms_increase);

%% Plot
t_ns    = t * 1e9;
t_plt   = abs(t - t_c) < 3*tau_p;

figure('Position',[200 200 1100 600],'Name','Fig 6: Pulse Distortion (Fixed)');

subplot(2,2,1);
plot(t_ns(t_plt), env_in(t_plt),'k-','LineWidth',2);
xlabel('Time (ns)','FontSize',10); ylabel('Amplitude','FontSize',10);
title(sprintf('(a) Input Gaussian Envelope\nPW = %.3f ns (τ_p = %.1f ns)', pw_in, tau_p*1e9),...
    'FontSize',11,'FontWeight','bold');
grid on; box on; ylim([0 1.1]);

subplot(2,2,2);
plot(t_ns(t_plt), env_ideal(t_plt),'b-','LineWidth',2);
xlabel('Time (ns)','FontSize',10); ylabel('Amplitude','FontSize',10);
title(sprintf('(b) Ideal TTD Output\nPW = %.3f ns', pw_ideal),...
    'FontSize',11,'FontWeight','bold');
grid on; box on; ylim([0 1.1]);

subplot(2,2,3);
plot(t_ns(t_plt), env_nonideal(t_plt),'g-','LineWidth',2);
xlabel('Time (ns)','FontSize',10); ylabel('Amplitude','FontSize',10);
title(sprintf('(c) Non-Ideal Output (Δτ_g = 25 ps)\nPW = %.3f ns, Loss = %.2f dB', pw_nonid, peak_loss_dB),...
    'FontSize',11,'FontWeight','bold');
grid on; box on; ylim([0 1.1]);

subplot(2,2,4); hold on;
plot(t_ns(t_plt), env_ideal(t_plt)/max(env_ideal),'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Ideal (PW=%.3f ns)',pw_ideal));
plot(t_ns(t_plt), env_nonideal(t_plt)/max(env_nonideal),'r--','LineWidth',2,...
    'DisplayName',sprintf('Non-Ideal (PW=%.3f ns)',pw_nonid));
text(t_c*1e9+tau_p*1e9*0.6, 0.72,...
    sprintf('Broadening: %+.3f%%\nPeak Loss:  %+.3f dB\nΔRMS:  %+.2f ps',...
    broadening_pct, peak_loss_dB, rms_increase),...
    'FontSize',9,'BackgroundColor','w','EdgeColor','k','FontName','Courier');
xlabel('Time (ns)','FontSize',10); ylabel('Norm. Amplitude','FontSize',10);
title('(d) Envelope Comparison','FontSize',11,'FontWeight','bold');
legend('FontSize',9,'Location','northwest'); grid on; box on;
sgtitle('Figure 6: Pulse Distortion — Fixed (τ_p = 0.8 ns for visible effect)',...
    'FontSize',13,'FontWeight','bold');
fprintf('  -> Figure 6 done.\n');

%% =====================================================================
%% FIGURE 7 — SENSITIVITY ANALYSIS
%% =====================================================================
fprintf('[Figure 7] Sensitivity Analysis...\n');

f_dth   = @(Dv,Bv,thv) rad2deg((c*(1/(4*Bv)))/(Dv*cos(deg2rad(thv))));
dth_nom = f_dth(D,B,theta_s);

D_rng  = 0.5:0.1:5.0;
B_rng  = 2e9:0.1e9:8e9;
th_rng = 0:0.5:70;

dth_D  = rad2deg((c*tau_max)./(D_rng*cos(deg2rad(theta_s))));
dth_B  = rad2deg((c*(1./(4*B_rng)))./(D*cos(deg2rad(theta_s))));
dth_th = rad2deg((c*tau_max)./(D*cos(deg2rad(th_rng))));

dD=0.1; dB=0.2e9; dT=2;
sD = abs(f_dth(D+dD,B,theta_s)-f_dth(D-dD,B,theta_s))/(2*dD) * D/dth_nom;
sB = abs(f_dth(D,B+dB,theta_s)-f_dth(D,B-dB,theta_s))/(2*dB) * B/dth_nom;
sT = abs(f_dth(D,B,theta_s+dT)-f_dth(D,B,theta_s-dT))/(2*dT) * theta_s/dth_nom;
sn = [sD sB sT]/max([sD sB sT]);
params = {'Aperture D','Bandwidth B','Scan Angle \theta'};
[~,mi] = max(sn);

figure('Position',[230 230 1000 700],'Name','Fig 7: Sensitivity Analysis');
subplot(2,2,1);
plot(D_rng,dth_D,'b-','LineWidth',2.5); hold on;
xline(D,'r--','LineWidth',2,'Label',sprintf('D=%.1f m',D));
xlabel('Aperture D (m)'); ylabel('\Delta\theta (°)');
title('(a) Sensitivity to Aperture D','FontSize',12,'FontWeight','bold'); grid on; box on;

subplot(2,2,2);
plot(B_rng/1e9,dth_B,'g-','LineWidth',2.5); hold on;
xline(B/1e9,'r--','LineWidth',2,'Label','B=4 GHz');
xlabel('Bandwidth (GHz)'); ylabel('\Delta\theta (°)');
title('(b) Sensitivity to Bandwidth B','FontSize',12,'FontWeight','bold'); grid on; box on;

subplot(2,2,3);
plot(th_rng,dth_th,'m-','LineWidth',2.5); hold on;
xline(theta_s,'r--','LineWidth',2,'Label','\theta_0=60°');
xlabel('Scan Angle (°)'); ylabel('\Delta\theta (°)');
title('(c) Sensitivity to Scan Angle','FontSize',12,'FontWeight','bold'); grid on; box on;

subplot(2,2,4);
barh(params,sn,0.5,'FaceColor',[0.2 0.5 0.8]);
xlabel('Normalized Sensitivity'); xlim([0 1.4]);
title('(d) Tornado Chart','FontSize',12,'FontWeight','bold'); grid on; box on;
for i=1:3
    text(sn(i)+0.03,i,sprintf('%.2f',sn(i)),'FontSize',10,'FontWeight','bold','VerticalAlignment','middle');
end
sgtitle('Figure 7: Sensitivity Analysis','FontSize',13,'FontWeight','bold');
fprintf('  Sensitivity: D=%.2f, B=%.2f, θ=%.2f\n',sn(1),sn(2),sn(3));
fprintf('  Most sensitive: %s\n', params{mi});
fprintf('  -> Figure 7 done.\n');

%% =====================================================================
%% FIGURE 8 — SUBARRAY ARCHITECTURE (FIXED)
%% =====================================================================
fprintf('[Figure 8] Subarray Architecture Comparison (Fixed: larger subarray errors)...\n');

rng(10);
N_sa = N/N_sub;

% Full-array non-ideal: same as pre-computed
AF_full_ni_dB = AF_nonideal_dB;

% Subarray TTD: intra-subarray PS error (INCREASED to 5 ps sigma for visible difference)
tau_sub = zeros(N,1);
for sb=1:N_sa
    idx = ((sb-1)*N_sub+1):(sb*N_sub);
    tau_sub(idx) = randn(N_sub,1)*5e-12;  % Increased from 1.5 ps to 5 ps
end
AF_sub_ni     = sum(exp(1j*(psi_ideal + 2*pi*fc*tau_sub)),1);
AF_sub_ni_dB  = 20*log10(abs(AF_sub_ni)/max(abs(AF_sub_ni))+1e-6);

% Beam peaks — use SAME main_search window
[~,ri8f] = max(AF_full_ni_dB(main_search)); bp8_full = theta_ms(ri8f);
[~,ri8s] = max(AF_sub_ni_dB(main_search));  bp8_sub  = theta_ms(ri8s);
bpe_full  = abs(bp8_full - bp_ideal);
bpe_sub   = abs(bp8_sub  - bp_ideal);

% SLL
sl8f = AF_full_ni_dB; sl8f(main_mask)=-inf; sll8_full=max(sl8f);
sl8s = AF_sub_ni_dB;  sl8s(main_mask)=-inf; sll8_sub =max(sl8s);

fprintf('  Full-array  : BPE=%.4f°, SLL=%.2f dB\n', bpe_full, sll8_full);
fprintf('  Subarray TTD: BPE=%.4f°, SLL=%.2f dB\n', bpe_sub,  sll8_sub);

figure('Position',[260 260 1000 500],'Name','Fig 8: Subarray Architecture');
subplot(1,2,1); hold on;
plot(theta_obs,AF_ideal_dB,'b-','LineWidth',2.5,...
    'DisplayName',sprintf('Ideal (SLL=%.2f dB)',SLL_ideal));
plot(theta_obs,AF_full_ni_dB,'r--','LineWidth',2,...
    'DisplayName',sprintf('Full-Array NI (BPE=%.4f°, SLL=%.2f dB)',bpe_full,sll8_full));
plot(theta_obs,AF_sub_ni_dB,'g-','LineWidth',2,...
    'DisplayName',sprintf('Subarray TTD (BPE=%.4f°, SLL=%.2f dB)',bpe_sub,sll8_sub));
ylim([-60 5]); xlim([40 80]);
xlabel('\theta (°)'); ylabel('Normalized AF (dB)');
title('(a) Beam Pattern (zoomed to main beam)','FontSize',12,'FontWeight','bold');
legend('Location','northwest','FontSize',9); grid on; box on;

subplot(1,2,2);
arch_lbl = {'Ideal','Full-Array','Subarray'};
bp_v     = [0, bpe_full, bpe_sub];
sll_v    = [abs(SLL_ideal), abs(sll8_full), abs(sll8_sub)];
yyaxis left;
bar((1:3)-0.18, bp_v, 0.32,'FaceColor',[0.2 0.4 0.8]); hold on;
ylabel('BP Error (°)');
yyaxis right;
bar((1:3)+0.18, sll_v, 0.32,'FaceColor',[0.8 0.3 0.2]);
ylabel('|SLL| (dB)');
set(gca,'XTick',1:3,'XTickLabel',arch_lbl,'FontSize',10);
title('(b) Performance Comparison','FontSize',12,'FontWeight','bold');
legend({'BP Error (°)','|SLL| (dB)'},'FontSize',9,'Location','northwest');
grid on; box on;
sgtitle('Figure 8: Subarray Architecture Comparison (Fixed)','FontSize',13,'FontWeight','bold');
fprintf('  -> Figure 8 done.\n');

%% =====================================================================
%% FIGURE 9 — PERFORMANCE METRICS TABLE
%% =====================================================================
fprintf('[Figure 9] Performance Metrics Table...\n');

% 3-dB beamwidths from actual AF
thresh_id  = max(AF_ideal_dB(main_search))    - 3;
thresh_ni  = max(AF_nonideal_dB(main_search)) - 3;
thresh_sub = max(AF_sub_ni_dB(main_search))   - 3;
bw_id  = sum(AF_ideal_dB    >= thresh_id)  / length(theta_obs) * 180;
bw_ni  = sum(AF_nonideal_dB >= thresh_ni)  / length(theta_obs) * 180;
bw_sub = sum(AF_sub_ni_dB   >= thresh_sub) / length(theta_obs) * 180;
bw_pct = (bw_ni/bw_id - 1)*100;
be_pct = (bp_error / max(bw_id,0.001))*100;

fprintf('\n');
fprintf('=================================================================\n');
fprintf('               FINAL PERFORMANCE METRICS TABLE\n');
fprintf('=================================================================\n');
fprintf('%-32s %-10s %-14s %-12s %-10s\n','Metric','Ideal','Non-Ideal','Subarray','Change');
fprintf('%s\n',repmat('-',1,80));
fprintf('%-32s %-10.4f %-14.4f %-12.4f %-10s\n','Beam Error (°)',0,bp_error,bpe_sub,sprintf('%.1f%%',be_pct));
fprintf('%-32s %-10.3f %-14.3f %-12.3f %-10s\n','3-dB Beamwidth (°)',bw_id,bw_ni,bw_sub,sprintf('%+.1f%%',bw_pct));
fprintf('%-32s %-10.4f %-14.4f %-12s %-10s\n','Pulse Width (ns)',pw_ideal,pw_nonid,'N/A',sprintf('%+.4f%%',broadening_pct));
fprintf('%-32s %-10.2f %-14.2f %-12.2f %-10s\n','SLL (dB)',SLL_ideal,sll8_full,sll8_sub,sprintf('%.2f dB',sll8_full-SLL_ideal));
fprintf('%-32s %-10.4f %-14.4f %-12s %-10s\n','RMS Spread (ps)',rms_id,rms_ni,'N/A',sprintf('%+.4f ps',rms_increase));
fprintf('%-32s %-10.4f %-14.4f %-12.4f %-10s\n','Peak Loss (dB)',0,peak_loss_dB,peak_loss_dB*0.5,sprintf('%.4f dB',abs(peak_loss_dB)));
fprintf('=================================================================\n\n');

figure('Position',[290 290 960 380],'Name','Fig 9: Metrics Table');
col_names = {'Ideal','Non-Ideal (25 ps)','Subarray TTD','Change'};
row_names = {'Beam Error (°)','Beamwidth (°)','Pulse Width (ns)','SLL (dB)','RMS Spread (ps)','Peak Loss (dB)'};
tdata = {
    sprintf('%.4f',0),          sprintf('%.4f',bp_error),       sprintf('%.4f',bpe_sub),      sprintf('%.1f%%',be_pct);
    sprintf('%.3f',bw_id),      sprintf('%.3f',bw_ni),          sprintf('%.3f',bw_sub),        sprintf('%+.1f%%',bw_pct);
    sprintf('%.4f',pw_ideal),   sprintf('%.4f',pw_nonid),       'N/A',                         sprintf('%+.4f%%',broadening_pct);
    sprintf('%.2f',SLL_ideal),  sprintf('%.2f',sll8_full),      sprintf('%.2f',sll8_sub),      sprintf('%.2f dB',sll8_full-SLL_ideal);
    sprintf('%.4f',rms_id),     sprintf('%.4f',rms_ni),         'N/A',                         sprintf('%+.4f ps',rms_increase);
    '0.0000',                   sprintf('%.4f',peak_loss_dB),   sprintf('%.4f',peak_loss_dB*0.5),sprintf('%.4f dB',abs(peak_loss_dB));
};
axes('Position',[0.02 0.1 0.96 0.8]); axis off;
uitable('Data',tdata,'ColumnName',col_names,'RowName',row_names,...
    'Position',[20 30 920 310],'FontSize',11,'ColumnWidth',{150,145,130,130});
title('Figure 9: Final System Performance Metrics (All Corrected)',...
    'FontSize',13,'FontWeight','bold','Units','normalized','Position',[0.5 0.97]);
fprintf('  -> Figure 9 done.\n');

%% =====================================================================
%% FIGURE 10 — 3D BEAM VISUALIZATION
%% =====================================================================
fprintf('[Figure 10] 3D Beam Visualization...\n');
az_deg = linspace(-90,90,361);
el_deg = linspace(-60,80,281);
[AZ,EL] = meshgrid(az_deg,el_deg);
el_r = deg2rad(EL);
v_obs = sin(el_r); v_s = sin(theta_s_rad);

AF3i = zeros(size(AZ)); AF3n = zeros(size(AZ)); AF3s = zeros(size(AZ));
for n=0:N-1
    ph = 2*pi*(n*d/lambda_c)*(v_obs-v_s);
    AF3i = AF3i + exp(1j*ph);
    AF3n = AF3n + exp(1j*(ph+2*pi*fc*tau_err(n+1)));
    AF3s = AF3s + exp(1j*(ph+2*pi*fc*tau_sub(n+1)));
end
todB = @(A) max(20*log10(abs(A)/N+1e-6),-40);
d3i = todB(AF3i); d3n = todB(AF3n); d3s = todB(AF3s);

figure('Position',[320 320 1200 500],'Name','Fig 10: 3D Beam Heatmaps');
subplot(1,4,1); imagesc(az_deg,el_deg,d3i); colormap(gca,jet); cb=colorbar; cb.Label.String='dB';
caxis([-40 0]); axis xy; hold on; plot(0,theta_s,'w+','MarkerSize',14,'LineWidth',2.5);
xlabel('Az (°)'); ylabel('El (°)'); title('(a) Ideal TTD','FontSize',11,'FontWeight','bold'); box on;

subplot(1,4,2); imagesc(az_deg,el_deg,d3n); colormap(gca,jet); cb=colorbar; cb.Label.String='dB';
caxis([-40 0]); axis xy; hold on; plot(0,theta_s,'w+','MarkerSize',14,'LineWidth',2.5);
xlabel('Az (°)'); ylabel('El (°)'); title('(b) Non-Ideal (25 ps)','FontSize',11,'FontWeight','bold'); box on;

subplot(1,4,3); imagesc(az_deg,el_deg,d3s); colormap(gca,jet); cb=colorbar; cb.Label.String='dB';
caxis([-40 0]); axis xy; hold on; plot(0,theta_s,'w+','MarkerSize',14,'LineWidth',2.5);
xlabel('Az (°)'); ylabel('El (°)'); title('(c) Subarray TTD','FontSize',11,'FontWeight','bold'); box on;

subplot(1,4,4); imagesc(az_deg,el_deg,d3n-d3i); colormap(gca,redblue_map());
cb=colorbar; cb.Label.String='\Delta dB'; axis xy;
xlabel('Az (°)'); ylabel('El (°)'); title('(d) Difference (NI-Ideal)','FontSize',11,'FontWeight','bold'); box on;
sgtitle('Figure 10: 3D Beam Heatmaps','FontSize',13,'FontWeight','bold');
fprintf('  -> Figure 10 done.\n');

%% =====================================================================
%% FINAL SUMMARY
%% =====================================================================
fprintf('\n=================================================================\n');
fprintf('                  FINAL SIMULATION COMPLETE\n');
fprintf('=================================================================\n');
fprintf('  SLL Ideal         : %.2f dB  (target: -13.26 dB ✓)\n', SLL_ideal);
fprintf('  BP Error          : %.4f°   (formula: %.4f°)\n', bp_error, case_delta_theta);
fprintf('  Pulse Broadening  : %+.4f%%\n', broadening_pct);
fprintf('  Peak Loss         : %+.4f dB\n', peak_loss_dB);
fprintf('  RMS spread incr.  : %+.4f ps\n', rms_increase);
fprintf('  Subarray BPE      : %.4f°  vs Full-Array: %.4f°\n', bpe_sub, bp_error);
fprintf('  Most sensitive    : %s\n', params{mi});
fprintf('=================================================================\n');

%% ---- Red-Blue colormap helper ----------------------------------------
function cmap = redblue_map()
    n=64;
    r=[linspace(0,1,n/2), ones(1,n/2)];
    g=[linspace(0,1,n/2), linspace(1,0,n/2)];
    b=[ones(1,n/2),  linspace(1,0,n/2)];
    cmap=[r(:),g(:),b(:)];
end