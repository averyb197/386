%for every valid X value 0.1:1:
%   make a 3vec of xi
%   subtract 3vec from x
%   compute ri with a row from n
%
%round ris to 0.1
%count number of entire rounded ri rows that equal target
%

num_trials = 10000;
n1 = randn(num_trials, 3);
n2 = randn(num_trials, 3);

lower = 0;
upper = 2;
%increment from lower to upper, step 0.1
x_cand = lower:0.1:upper;
%number of x values to check
stupidity = size(x_cand);

% will store the total number of matches found for each x distance
tamp_match_count = zeros(size(x_cand));
true_match_count = zeros(size(x_cand));
xis = 0.22*[1, 2, 3];

%for each x to test
for i=1:stupidity(2)
    x = x_cand(i);
    %matlab docs say this is valid to do to get mean of 1, to me it makes
    %way more sense to just i dont know, add like 2 lines into the function
    %definition that way you can just pass the mean as an argument like you
    %do in every other programming language known to man but whatev ig.
    n1 = .5 + randn(num_trials, 3);
    n2 = .5 + randn(num_trials, 3);
    targets = [0.8, 0.8, 0.9; 0.6, 0.6, 0.7];
    
    %1 dec place
    %10000x3
    ri_tamp = round(exp(-abs(x-xis)./n1), 1);
    ri_true = round(exp(-abs(x-xis)./n2), 1);
    
    %all returns 1 if all elements are nonzero, or 0 is one or more is
    %zeros, so everytime a row perfectly matches the given target, return 1
    %check along 2 axis
    tamper_matches = all(ri_tamp==targets(1,:), 2);
    true_matches = all(ri_true==targets(2,:), 2);
    
    %Record count all the perfect matches
    tamp_match_count(i) = sum(tamper_matches);
    true_match_count(i) = sum(true_matches);
end

[max_tamp, max_idx] = max(tamp_match_count);
[max_true, tru_idx] = max(true_match_count);

da_title = sprintf("Likelihood Response: Tampered Max = %.2fpc; True Max = %.2fpc", x_cand(max_idx), x_cand(tru_idx));

figure;
plot(x_cand, tamp_match_count); hold on
plot(x_cand, true_match_count);
hold off
legend(["Tampered", "True"]);
title(da_title);
xlabel("Distance From Star (Parsec)");
ylabel("# Matching Target");

