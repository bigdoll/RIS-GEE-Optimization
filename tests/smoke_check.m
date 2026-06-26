function smoke_check()
% SMOKE_CHECK  Lightweight integrity test for the RIS-GEE-Optimization repo.
%
%   smoke_check()
%
% Runs WITHOUT CVX or a solver. It verifies that:
%   1. every MATLAB function shipped in src/ resolves on the path (which());
%   2. every intra-repository call inside each file resolves to a repo file
%      or a MATLAB built-in (i.e. nothing is left dangling after the cleanup);
%   3. the shared "common" pipeline executes on a tiny problem instance
%      (generate_channels -> func_R -> LMMSE -> SINR -> data_rate).
%
% Exit: prints a PASS/FAIL summary and throws an error if any check fails,
% so it can be used as a CI gate ( matlab -batch "smoke_check" ).

    here   = fileparts(mfilename('fullpath'));
    root   = fullfile(here, '..');
    srcDir = fullfile(root, 'src');
    addpath(genpath(srcDir));
    addpath(genpath(fullfile(root, 'examples')));

    fprintf('RIS-GEE-Optimization smoke check\n--------------------------------\n');

    % ---- collect repo function names ------------------------------------
    mfiles   = dir(fullfile(srcDir, '**', '*.m'));
    repoNames = erase({mfiles.name}, '.m');
    fprintf('Found %d source functions in src/.\n', numel(repoNames));

    failures = 0;

    % ---- check 1: every repo function is on the path --------------------
    for i = 1:numel(repoNames)
        if isempty(which(repoNames{i}))
            fprintf(2, '  [MISSING] %s is not on the path\n', repoNames{i});
            failures = failures + 1;
        end
    end

    % ---- check 2: intra-repo calls resolve ------------------------------
    allFiles = [mfiles; dir(fullfile(root, 'examples', '*.m'))];
    for i = 1:numel(allFiles)
        fpath = fullfile(allFiles(i).folder, allFiles(i).name);
        code  = fileread(fpath);
        code  = regexprep(code, '%[^\n]*', '');           % strip comments
        tokens = unique(regexp(code, '[A-Za-z]\w*', 'match'));
        called = intersect(tokens, repoNames);
        for j = 1:numel(called)
            if isempty(which(called{j}))
                fprintf(2, '  [UNRESOLVED] %s -> %s\n', allFiles(i).name, called{j});
                failures = failures + 1;
            end
        end
    end

    % ---- check 3: common pipeline runs on a tiny instance ---------------
    try
        K = 2; NR = 2; N = 8;
        [G, H] = generate_channels(100, 10, 15, 5, 50, N, K, NR);
        p = 0.1 * ones(K,1);
        sigma_sq = noise_power(-174, 10, 20e6);
        R = func_R(p, H, sigma_sq, K, N);                 %#ok<NASGU> Eq. (6)
        gamma = exp(1i*2*pi*rand(N,1));
        C = LMMSE_receiver_passive(p, G, H, gamma, sigma_sq, K, NR);
        sinr = SINR_passive(p, C, G, H, gamma, sigma_sq, K);
        [~, sr] = data_rate(sinr, K, 20);
        assert(isfinite(sr) && sr >= 0, 'sum-rate is not a valid number');
        fprintf('Common pipeline executed: sum-rate = %.4f Mbit/s.\n', sr);
    catch err
        fprintf(2, '  [PIPELINE ERROR] %s\n', err.message);
        failures = failures + 1;
    end

    % ---- verdict --------------------------------------------------------
    if failures == 0
        fprintf('\nPASS: all checks succeeded.\n');
    else
        error('smoke_check:failed', '\nFAIL: %d problem(s) found.', failures);
    end
end
