astyle()
%Nonlinear activation function also restricts voltages to be on [-1, 1]
function hov = h(v)
    hov = 0 + (-20.*v.^7 + 70.*v.^6 - 84.*v.^5 + 35.*v.^4) .* (v <= 1 & v >= 0) + 1 .* (v > 1);
end

%v is vector of voltages (1=L, 2=R) v = [R1, E1, L1, C1, R2, E2, L2, C2]
function daderiv = dvdt(t, v, GLC)
    %GR=GE=GL=GC=3.5
    GR = 3.5; 
    %GR=3.5 GE=0.875 GL=0.35 GC=3.5
    %GTC = tonic drive to C - variation will 
    GTC = 3.5;
    %Other Tonic Drives
    GT = [3.5; 0.875; 0.35; GTC; 3.5; 0.875; 0.35; 3.5]; 
    %Resting membrane potential
    VR = 0;
    %Tonic synaptic pote
    VT = 1;
    %synaptic potentials
    %VR = VE = 1 VL = VC = -1
    Vsyn = [1; 1; -1; -1; 1; 1; -1; -1]; 
    %2nd figs, graph 
    %Synaptic weight matrix, Gsyn(i, j) = Gj->i
    Gsyn = [0, 55, 0, 0, 0, 0, 0, 55;
            15, 0, 0, 0, 0, 0, 0, 35;
            5.5, 35, 0, 0, 0, 0, 0, 35;
            7, 35, GLC, 0, 0, 0, 0, 35;
            0, 0, 0, 55, 0, 55, 0, 0;
            0, 0, 0, 35, 15, 0, 0, 0;
            0, 0, 0, 35, 5.5, 35, 0, 0;
            0, 0, 0, 35, 7, 35, GLC, 0
    ];
    daderiv = zeros(8, 1);
    %i is neuron computing deriv for
    for i=1:8
        %could I do this without a for loop? Yes. Did I break it and was
        %way too close to the deadline so decided to make it slower? also
        %yes.
        sum_term = 0;
        for j = 1:8
            %Sum over all inputs into neuron i
            sum_term = sum_term + Gsyn(i, j)*h(v(j))*(Vsyn(j) - v(i));
        end
        daderiv(i) = GR * (VR - v(i)) + GT(i) * (VT - v(i)) + sum_term;   
    end
end

%returns vectors with GLC values and stable limit cycle frequencies
function [GLC_vals, limit_cycle_freqs, lc_mins, lc_maxes] = calc_limit_cycle_freqs(GLC_lb, GLC_ub)
    %Initial conditions implied from figures in paper
    v0 = [0.96, .8, .95, .15, -0.65, -0.65, -0.65, -0.55];
    %v0 = [0,0,0,0,0,0,0,0];
    %how many sample of the parameter being varied, plot of frequencies is
    %discrete in nature, so doesnt need to be too high
    param_res = 100;
    %GLC values to input
    GLC_range = linspace(GLC_lb, GLC_ub, param_res);
    %Number of time steps for solver
    %Need to keep timesteps constant for freq analysis bc it needs a constantly
    %sampled signal to extract accurate frequency
    solv_res = 1000;
    %trial time
    T = 5;
    timeline = linspace(0, T, solv_res);
    %To store output values from astyle_fft()
    freqy = zeros(param_res, 1);
    lc_mins = zeros(param_res, 1);
    lc_maxes = zeros(param_res, 1);
    for i=1:param_res
        %Wrapper function that sets GLC and makes passable solver
        curr_dvdt = @(t, v) dvdt(t, v, GLC_range(i));
        [~, vot] = ode45(curr_dvdt, timeline, v0);
        %Left C voltages 
        C1 = vot(:, 4);
        %Pass C to fft function to find the 
        freq_stuff = astyle_fft(C1, T);
        %principle frequency
        freqy(i) = freq_stuff{3};
        lc_mins(i) = min(C1(500:end));
        lc_maxes(i) = max(C1(500:end));
    end
    limit_cycle_freqs = freqy;
    GLC_vals = GLC_range;
end

function [GLC_range, stable_eq_pts] = find_stable_eq()
    GLC_ub = 200;
    GLC_lb = 0;
    param_res = 100;
    GLC_range = linspace(GLC_lb, GLC_ub, param_res);
    bard = [0.96, .8, .95, .15, -0.65, -0.65, -0.65, -0.55];
    ics = [ones(1, 8); 0.2*ones(1, 8); ones(1, 8)];
    num_ics = size(ics, 1);
    tspan = [0, 5];
    stable_eq_pts = zeros(num_ics, param_res, 1);
    for i=1:num_ics
        curr_ics = ics(i, :)';
        for k=1:param_res
            curr_dvdt = @(t, v) dvdt(t, v, GLC_range(k));
            [~, curvot] = ode45(curr_dvdt, tspan, curr_ics);
            pt = curvot(end, 4);
            curr_ics(4) = pt;
            stable_eq_pts(i, k) = pt;
        end
    end
end
T=5;
%[GLC_range, stable_pts] = find_stable_eq();
% save("stable_pts.mat", "stable_pts");
load("stable_pts.mat", "stable_pts");
% figure;
% subplot(2, 1, 1);
% plot(GLC_range, stable_pts);
% legend(["1", "2", "3"])

%load("limit_cycle_freqs_GLC.mat", "freqy");

% subplot(2, 1, 2);
% plot(GLC_range, freqy); hold on
% xline(24.2424, "--r", "LineWidth",4); hold off
% title("WITHLINE");
 
%[GLC_range, freqy, lc_mins, lc_maxes] = calc_limit_cycle_freqs(0, 200);
%Computation takes way too long
%save("limit_cycle_freqs_GLC.mat", "freqy", "GLC_range", "lc_mins", "lc_maxes");

load("limit_cycle_freqs_GLC.mat", "freqy", "GLC_range", "lc_mins", "lc_maxes");

%Model with standard params, default GLC=35
v0 = [0.96, .8, .95, .15, -0.65, -0.65, -0.65, -0.55];
std_tspan = [0 T];
standard_model = @(t, v) dvdt(t, v, 35);
[stan_t, stan_v] = ode45(standard_model, std_tspan, v0);

figure('Position', [0, 0, 1200, 800], 'Name','Response Suite');
subplot(3, 2, 1);
%left neurons from standard model
plot(stan_t, stan_v(:, 1:4));
legend(["R", "E", "L", "C"]);
title("Standard Parameters - LEFT");
xlabel("time (s)")
ylabel("Voltage")
axis auto

subplot(3, 2, 3);
%right neurons from standard model
plot(stan_t, stan_v(:, 5:8));
legend(["R", "E", "L", "C"]);
title("Standard Parameters RIGHT");
xlabel("time (s)")
ylabel("Voltage")
axis auto

peaks = zeros(4, 1);
%All neurons have the same frequency response, just different phase, so
%only need to calculate for one of them to get freq for all
subplot(3, 2, 5);
%pretend to care abt C specifically
C1 = stan_v(:, 4);
freq_test = astyle_fft(C1, T);
freq_ax = freq_test{1};
amp_spectra = freq_test{2};
peak_freq = freq_test{3};
plot(freq_ax, amp_spectra); 
hold off
xlim([0 10])
title_string_freq = sprintf("Frequency Response; Peak = %.2fHz", peak_freq);
xlabel("Frequency (Hz)");
ylabel("Relative Strength");
title(title_string_freq);

subplot(3, 2, 2);
C1 = stan_v(:, 4);
C2 = stan_v(:, 8);
plot(stan_t, C1); hold on
plot(stan_t, C2); hold off
title("Left/Right C");
legend(["Left", "Right"]);
xlabel("time (s)");
ylabel("Voltage");

subplot(3, 2, 4);
plot(GLC_range(13:83), freqy(13:83)); hold on
xline(165.657, "--r", "LineWidth",2);
xline(24.2424, "--r", "LineWidth",2); hold off
xlim([0 200]);
axis auto
title("Limit Cycle Frequencies");
xlabel("(L->C) Synaptic Conductance");
ylabel("Freq. (Hz)")
xlim([0 200]);


subplot(3, 2, 6);
hold on
lx = [24.2424, 30.303];
ly = [.28266, .16374];
plot(lx, ly, "b");
xline(24.2424, "--r", "LineWidth",2); 
xline(165.657, "--r", "LineWidth",2);
plot(GLC_range(10:83), lc_mins(10:83), "g"); 
plot(GLC_range(83:end), lc_mins(83:end), "b");
plot(GLC_range(83:end), lc_maxes(83:end), "b");
plot(GLC_range(1:13), stable_pts(:, 1:13), "b");
plot(GLC_range(16:end), stable_pts(3, 16:end), "b");
plot(GLC_range(13:83), lc_maxes(13:83), "Color","g"); 
title("Bifurcation Diagram");
xlabel("(L->C) Synaptic Conductance");
ylabel("Voltage of C")
h1 = plot(NaN, NaN, 'b'); 
h2 = plot(NaN, NaN, 'g'); 
legend([h1, h2], {'Stable Fixed', 'Stable Limit Cycle'});
hold off
axis auto




