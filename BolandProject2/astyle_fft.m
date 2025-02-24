%Input vector with sample and then duration of sample
function freq_info = astyle_fft(vectro_tull, T) %put in a vector to analyze frequency
    freq_spec = fft(vectro_tull); %frequency spectra
    n = length(vectro_tull);  %num sample
    norm_fs = abs(freq_spec/n); %normalized frequency spectra, divide by number of elemens for fourier reasons
    fs_full = norm_fs(1:floor(n/2)+1); %only want 1st half bc real signal, 2nd half redundant
    fs_full(2:end-1) = 2*fs_full(2:end-1); % double that parts we care about
    [~, peaki] = max(fs_full); %maximum freq, idx ie amplitude
    f = n/T; % avoid fuckery f = num_samp/time interval
   % n = length(fs_full);  %num sample
    freqax = (0:floor(n/2)) * (f / n); %axis for freq graph, divide by 2 bc fourier transform is symmetric and only care abt 
    % positive for a real signal
    amps = fs_full; %renaming bc i hate myself
    %disp(length(amps));
    max_freq = freqax(peaki); % see above 
    freq_ax = freqax; % see above 
    %MAXI HAS NO FRIENDS 
    freq_info = {freqax, amps, max_freq, peaki}; % python list but maybe worse?
end
% 
% %TO MAKE THE PLOT I KNOW YOU WANT TO MAKE: DO PLOT(1, 2) AND 3 IS THE VALUE
% %OF PEAK FREQ
