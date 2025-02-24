function hov = h(v)
    hov = 0 + (-20.*v.^7 + 70.*v.^6 - 84.*v.^5 + 35.*v.^4) .* (v <= 1 & v >= 0) + 1 .* (v > 1);
end

function daderiv = dvdt(t, v)
    %GR=GE=GL=GC=3.5
    %GR = 3.5 * ones(8, 1); 
    GR = 3.5; %* ones(8, 1);

    % GR=3.5 GE=0.875 GL=0.35 GC=3.5
    GT = [3.5; 0.875; 0.35; 3.5; 3.5; 0.875; 0.35; 3.5]; 
    VR = 0;%zeros(8, 1); 
    VT = 1;%;ones(8, 1); 
    %VR = VE = 1 VL = VC = -1
    Vsyn = [1; 1; -1; -1; 1; 1; -1; -1]; 

    % EL=EC=LC=CE=CL=CC=35
    % CR=ER=55
    % RE=15
    % RL=5.5
    % RC=7
    % R E L C R E L C
    % Gsyn = [ 
    %     0,  15, 5.5, 0,     0,  0,  0,   7;
    %     55, 0,  35,  0,     0,  0,  0,   35;
    %     0,  0,  0,   0,     0,  0,  0,   35;
    %     55, 35, 35,  0,     0,  0,  0,   35;
    %     0,  0,  0,   7,     0,  15, 5.5, 0;
    %     0,  0,  0,   35,    55, 0,  35,  0;
    %     0,  0,  0,   35,    0,  0,  0,   0;
    %     0,  0,  0,   35,    55, 35, 35,  0;
    % ];
    %Gsyn = Gsyn';

    Gsyn = [0,   55, 0,  0,  0, 0,  0, 55;
            15,  0,  0,  0,  0, 0,  0, 35;
            5.5, 35, 0,  0,  0, 0,  0, 35;
            7,   35, 35, 0,  0, 0,  0, 35;
            0,   0,  0,  55, 0  55, 0, 0;
            0,   0,  0,  35, 55, 0, 0, 0;
            0,   0,  0,  35, 5.5, 35, 0, 0;
            0, 0, 0, 35, 7, 35, 35, 0  
    ]';

    hv = h(v); 
    daderiv = zeros(8, 1);

    for i=1:8
        %sum over all of its connections, use (:, i) all other -> i, 
        % (i, :) == i -> all others 
        syns = sum(Gsyn(:, i) .*hv .* (Vsyn - v(i)));
        daderiv(i) = GR * (VR - v(i)) + GT(i) * (VT - v(i)) + syns;
    end
end

v0 = [0.96, .8, .95, .15, -0.65, -0.65, -0.65, -0.55];

timeline = [0 3];

[tt, vot] = ode45(@dvdt, timeline, v0);

figure;
subplot(2, 1, 1);
plot(tt, vot(:, 1:4));
legend(["R", "E", "L", "C"]);
title("LEFT");

subplot(2, 1, 2);
plot(tt, vot(:, 5:8));
legend(["R", "E", "L", "C"]);
title("RIGHT");




% profile off;
% profile viewer



