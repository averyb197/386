% loop of t -> do once per t/100 times per trial:
% - calculate v(t)
% - then calc del(t)
% - then update w from tau=0:t -> subloop 
%delta candidate --- d(t) = r(t) v(t) - v(t-1)

%scalar: t; arr: u, w_arr
%output 1d array: w
function wtau = update_w(t, u, w_arr, deltat)
    epsilon = 0.05;
    for idx=1:t
        w_arr(idx) = w_arr(idx) + epsilon*deltat*u(t-idx+1);
    end
    wtau = w_arr;
end

%take in scalar: v, t, T; arr: r
%output scalar delta(t)

%input scalar: t; array: r, v
function rdeltat = rdelta(t, r, v)
    % D = 0.05;
    % deltat = r(t) + v(t+1) - v(t) + D;
    % rdeltat = max(deltat, D);
    rdeltat = r(t) + v(t+1) - v(t);
end

%input scalar t; array: w, u
%output scalar v = v(t)
function vee = vfunc(t, w, u)
    vsum = 0;
    if t==101
    else
    for tau=1:t 
        % disp('start sum block')
        % w(tau)
        % u(t-tau+1)
        % disp(t-tau+1)
        % disp(u(t-tau+1))
        %disp(w(tau)*u(t-tau+1));
        vsum = vsum + w(tau)*u(t-tau+1);
        %disp('end sum block')
    end
    end
    vee = vsum;
end

%%%output single scalar v, w array, updated up to t
%calc v(t)
%calc delta(t)
%update all w from 0:t
%output scalar v, w array
function vwout = take_step(t, w, r, u, v)
    real_tau = 5;
    T = 100; %total trial time
    voft = vfunc(t, w, u); %scalar v=v(t)
    %doft = delfunky(t, r, voft, T); %scalar delta(t)
    doft = rdelta(t, r, v);

    new_w = update_w(t, u, w, doft);
    %scalar voft, dt needs to be added to array in trial, w is perfect len
    v(t) = voft;
    vwout = {voft, new_w, doft};
end

function trial_data = run_trial(w, r, u)
    vout = zeros(101, 1);
    dovt = zeros(101, 1);
    w_surf = zeros(101, 101); %100x100
    for t=1:100
        step_data = take_step(t, w, r, u, vout);
        vout(t) = step_data{1}; %v(t)
        step_data{1};
        w = step_data{2};
        w_surf(:, t) = step_data{2}; %output is whole w vector
        dovt(t) = step_data{3};%scalar d(t)
    end
    trial_data = {vout, dovt, w_surf};
end

real_tau = 5;
r = arrayfun(@(t) ((t>=50) * ((t - 50)/real_tau)*exp(-(t - 50)/real_tau)), 1:101);
u = arrayfun(@(t) (t==1)*1 , 1:101);
wi = zeros(101, 1);

full_trial = run_trial(wi, r, u);
t = linspace(0, 100, 101);
vovert = full_trial{1};
delt_ovt = full_trial{2};
w_trial = full_trial{3};
last_w = w_trial(:, end-1);

w_surf = zeros(101, 3000);
next_w = wi;
delta_surf = zeros(101, 3000);
v_surf = zeros(101, 3000);
for i=1:3000
    ith_info = run_trial(next_w, r, u);
    ith_delta = ith_info{2};
    ith_v = ith_info{1};
    v_surf(:, i) = ith_v;
    delta_surf(:, i) = ith_delta;
    ith_w = ith_info{3};
    next_w = ith_w(:, end-1);
    w_surf(:, i) = next_w;
end

figure;
subplot(3, 1, 1);
mesh(delta_surf);
title("Delta Surface");

subplot(3, 1, 2);
mesh(w_adapt);
title("W surface")

subplot(3, 1, 3);
mesh(v_surf);
title("V Surface");

figure;
subplot(3,1, 1);
plot(t, v_surf(:, 1)); hold on
plot(t, v_surf(:, end));
hold off
legend(["Trial 1", "Trial 3000"])
title("V")
xlabel("time")
ylabel("trials")
zlabel("V")

subplot(3, 1, 2);
plot(t, w_surf(:, 1)); hold on
plot(t, w_surf(:, end));
hold off
legend(["Trial 1", "Trial 3000"])
title("W")
xlabel("time")
ylabel("trials")
zlabel("Delta")

subplot(3, 1, 3);
plot(t, delta_surf(:, 1)); hold on
plot(t, delta_surf(:, end));
hold off
legend(["Trial 1", "Trial 3000"])
title("Delta");
xlabel("time")
ylabel("trials")
zlabel("Delta")








%function run trial -> take 100 steps, input w array, output entire v array
%over trial
    


 
    
































































    
    