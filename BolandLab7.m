sigma_arr = [8 16; 8 32;4 32];

function dwdt = daODE(t, w, K)
    %make w col vec into matrix to live as true self
    %NO Q BECAUSE IT IS MEANINGLESS AND WE HATE IT
    W = reshape(w, 512, 2);
    dwdtMAT = K*W;
    %flatten back to column vec
    dwdt = dwdtMAT(:);
end

function freq_info = fourier_anal(vectro_tull) %put in a vector to analyze frequency
    vectro_tull = vectro_tull - mean(vectro_tull); %remove DC offset
    vectro_tull = detrend(vectro_tull); % get rid of stuff near 0 -> long term trends
    %vectro_tull = highpass(vectro_tull, 1, 512);
    freq_spec = fft(vectro_tull); %frequency spectra
    n = length(vectro_tull);  %num sample

    norm_fs = abs(freq_spec/n); %normalized frequency spectra, divide by number of elemens for fourier reasons
    fs_full = norm_fs(1:floor(n/2)+1); %only want 1st half bc real signal, 2nd half redundant
    fs_full(2:end-1) = 2*fs_full(2:end-1); % double that parts we care about
    
    [maxi, peaki] = max(fs_full); %maximum freq, idx ie amplitude y i
    
    f = 512; %call freq 512 even though not really time dependent, needs to be f/2 = max sample rate to avoid fuckery
    freqax = (0:floor(n/2)) * (f / n); %axis for freq graph, divide by 2 bc fourier transform is symmetric and only 
    % care abt positive for a real signal
    amps = fs_full; %renaming bc i hate myself
    max_freq = maxi; % see above 
    freq_ax = freqax; % see above 
    freq_info = {freqax, amps, maxi, peaki}; % python list but maybe better?
end

%Method for generating karyotype looking plots
%subtract the absolute value of col2 from col1. If col has a larger
%magnitude, then the overall value is negative, which I consider right eye
%dominance, and vice versa where positive is left eye dominance. Then I
%normalize the 1D array to use for generating the blur plot

function sweet = response_suite(sigma_arr, num_steps, blur_plot_thicc)
    timeline = linspace(0, num_steps, 100);
    W0 = randn(1024,1); % need W to be flat to treat like 1d vec for ode45
    W0_pre_blur = reshape(W0, 512, 2);
    %make as horrendous as possible
    W0_blur = arrayfun(@(i) (W0_pre_blur(i, 1) < W0_pre_blur(i, 2)) * (abs(W0_pre_blur(i, 1)) - abs(W0_pre_blur(i, 2))) + ...
                          (W0_pre_blur(i, 1) > W0_pre_blur(i, 2)) * (-(abs(W0_pre_blur(i, 1)) - abs(W0_pre_blur(i, 2)))), 1:512);

    w_ind = arrayfun(@(i) i, 1:512); % index for weights
    [kcol, krow] = ndgrid(1:512, 1:512); %to use for making Ws
    
    normed_W0_blur = (W0_blur - min(W0_blur)) / (max(W0_blur) - min(W0_blur)); %normalize W0
    grayWi = repmat(normed_W0_blur, blur_plot_thicc, 1); %repeat matrix blur_plot_thicc # of times along 1 axis
    
    %figure;
    figure('Position', [0, 0, 1400, 800], 'Name','Response Suite');
    
    subplot(3, 4, 1);
    imshow(grayWi, 'Colormap',gray);
    title("Initial Weights, Random Normal");

    fig_count = 0;
    wi = arrayfun(@(i) i, 1:512); % index for weights

    for i=1:3
        sigmas = sigma_arr(i, :);
        K = arrayfun(@(i, j) (1/sigmas(1)^2) * exp(-((i - j)^2) / (2 * sigmas(1)^2)) - ...
        (1/sigmas(2)^2) * exp(-((i - j)^2) / (2 * sigmas(2)^2)), krow, kcol); %make K based on current sigmas

        [t, wot] = ode45(@(t,w) daODE(t, w, K), timeline, W0);
        w_final = wot(end, :); % final weight values
        reshape_w = reshape(w_final, 512, 2); %reshape w into actual shape 512x2 

        last_w = abs(reshape_w(:, 1)) -  abs(reshape_w(:, 2)); %take absolute difference, neg for right dom, pos for left dom
        last_w_norm = (last_w - min(last_w)) / (max(last_w) - min(last_w)); % normalize absolute diff
        w_blur = last_w_norm; %instead of writing code effectively, just change var names

        blur_norm = (w_blur - min(w_blur)) / (max(w_blur) - min(w_blur)); %normalize again bc too lazy to change car names
        alt_gray = repmat(blur_norm', blur_plot_thicc, 1); %repeat 

        %output = {frequency axis, amplitude spectra, peak frequency, peak
        %freq amplitude}
        fft_info = fourier_anal(w_blur);

        freq_ax = fft_info{1};
        amp_spectra = fft_info{2};
        peak_freq = fft_info{4};

        [KV, KD] = eig(K); %eigs of K

        eig_fft = fourier_anal(KV(:, end)); %get freq of principle eigenvector

        eig_spectra = eig_fft{2};
        peak_eig_freq = eig_fft{4};

        subplot(3, 4 , fig_count+1);
        imshow(grayWi, 'Colormap',gray); %USE IMSHOW USE IMSHOW USE IMSHOW +++++ DO NOT USE PLOT HOE
        title("Initial Weights, Random Normal");
        c = colorbar;
        c.Ticks = [c.Limits(1), c.Limits(2)]
        c.TickLabels = {'Right', 'Left'}
        xlim([0, 512]);
        xlabel('Cortical Cells')
        axis normal

        subplot(3, 4, fig_count+2);
        last_w_norm = (last_w - min(last_w)) / (max(last_w) - min(last_w));
        gray_wf = repmat(blur_norm', blur_plot_thicc, 1); %REP THE TRANSPOSE IDIOT 
        imshow(gray_wf, 'Colormap',gray);
        ftitle = sprintf('Final Weights \\sigma_1 = %d \\sigma_2 = %d', sigmas(1), sigmas(2));
        title(ftitle, 'Interpreter','tex');
        xlabel('Cortical Cells')
        xlim([0, 512]);
        axis normal

        subplot(3, 4, fig_count+3);
        bor_f = arrayfun(@(i) 1 * ((abs(reshape_w(i, 1)) > abs(reshape_w(i, 2)))) + 0 , 1:length(reshape_w));
        plot(wi, bor_f);
        title('Ocular Dominance');
        xlim([0, 512]);
        axis normal

        subplot(3, 4, fig_count+4);
        plot(freq_ax, amp_spectra); hold on %plot w spectra
        plot(freq_ax, eig_spectra); hold off %plot 
        title_string = sprintf('Spectra, Peak W freq = %.2f \n Principle Eigenvector Peak = %.2f', peak_freq, peak_eig_freq);
        set(legend, 'FontSize', 8);
        title(title_string);
        axis normal
        xlim([0, 50]);
        xlabel("Relative Frequency");
        ylabel("Frequency Strength");
        legend({'Normalized Abs. Weight Diff.','Principle Eigenvector'},'Location','northeast','Orientation','vertical');

        fig_count = fig_count + 4; % do it again

    end

end

response_suite(sigma_arr, 500, 100)

        


                




    