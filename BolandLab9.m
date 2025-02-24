sigma = 2;
num_samples = 1000;
time_tau = 10;
timeline = 1:num_samples;
D_true = timeline .* exp(-timeline/time_tau);
%multiply by 2 because variance = sigma^2 = 4 so sigma=2
s = 2*randn(num_samples, 1);

%use full for 0 padding when index outside bounds treat as 0 when <0 or 1000< like
%true filter should
r_true_full = conv(D_true, s, 'same');
r_true = r_true_full(1:1000);

%Qrs*1000 = int tau=0:inf [r(t)s(t+tau)] = sum tau=0:1000 r(t)s(t+tau) =
%conv(r(t), flip(s)) --- flip since s(t+tau) same as reversing order bc of
%something????
%%(t+tau) goes from t:1000+t maybe????
%D(t) = Qrs(-t)/sigma^2 
Qrs_real = conv(r_true, flip(s), "same")/1000;
dest = Qrs_real/2;
D_est = zeros(1000, 1);

D_est(1:100) = dest(1:100);
D_est = D_est';
%new data to check filter
test_s = 2*randn(num_samples, 1);
test_r_true_full = conv(D_true, s, 'same');
test_r_true = test_r_true_full(1:1000);

test_r_est_full = conv(D_est, s, 'same');
test_r_est = test_r_est_full(1:1000);

test_abs_error = abs((test_r_est - test_r_true)) ./ abs(test_r_true);
median_error = median(test_abs_error)*100;
mean_abs_error = mean(test_abs_error);

num_over_50 = sum(test_abs_error > .5);
prop_50 = (num_over_50 / 1000) * 100;

error = sum( (test_r_est - test_r_true).^2 )/1000;

figure('Position', [0, 0, 1400, 800], 'Name','Response Suite');
subplot(4, 1, 1);
plot(1:1000, D_true); hold on
plot(1:1000, D_est);
hold off
title("Filter");
%legend(["True", "Estimated", 'Wrong']);

subplot(4, 1, 2);
plot(1:1000, test_s);
title("Twitch Intensity");

subplot(4, 1, 3);
plot(1:1000, test_r_true); hold on
plot(1:1000, test_r_est);
hold off
title("Responses");
%legend(["True", "Estimated", "Wrong"]);


subplot(4, 1, 4);
axis off;
title("Stats");
stats = sprintf("Trial Data:\nMedian Error: %.2f%% \n Prop. Error >50%%: %.d/1000 = %.2f%%\n Mean Absolute Error %.2f\nMean Square Error: %.2f", median_error, num_over_50, prop_50 , 100*mean_abs_error, error);
%text(0.5, 0.5, stats, "FontSize", 20)
text(0.5, 0.4, stats, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, "FontWeight", "bold");


