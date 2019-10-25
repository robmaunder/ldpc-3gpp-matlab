function plot_SNR_vs_A(A, R, BG, Modulation, rv_id_sequence, iterations, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
% PLOT_SNR_VS_A Plots Signal to Noise Ratio (SNR) required to achieve a
% particular Block Error Rate (BLER) as a function of block length, for
% LDPC codes.
%
%   target_BLER should be a real scalar, in the range (0, 1). The
%   simulation of each coding rate will continue until the BLER plot
%   reaches this value.
%
%   EsN0_start should be a real row vector, having the same length as the
%   vector of coding rates. Each value specifies the Es/N0 SNR to begin at
%   for the simulation of the corresponding coding rate.
%
%   EsN0_delta should be a real scalar, having a value greater than 0.
%   The Es/N0 SNR is incremented by this amount whenever
%   target_block_errors number of block errors has been observed for the
%   previous SNR. This continues until the BLER reaches target_BLER.
%
%   seed should be an integer scalar. This value is used to seed the random
%   number generator, allowing identical results to be reproduced by using
%   the same seed. When running parallel instances of this simulation,
%   different seeds should be used for each instance, in order to collect
%   different results that can be aggregated together.
%
%   See also MAIN_BLER_VS_SNR and MAIN_FAR
%
% Copyright © 2017 Robert G. Maunder. This program is free software: you
% can redistribute it and/or modify it under the terms of the GNU General
% Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version. This
% program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

% Default values
if nargin == 0
    A = (1000:1000:8000);
    R = 1/3;
    BG = 1;
    Modulation = 'QPSK';
    rv_id_sequence = [0];
    iterations = 50;
    target_block_errors = 100;
    target_BLER = 1e-2;
    EsN0_start = -2;
    EsN0_delta = 0.1;
    seed = 0;
end

% Seed the random number generator
rng(seed);

% Create a figure to plot the results.
figure
axes1 = axes;
title(['3GPP New Radio LDPC code, BG',num2str(BG),', ',Modulation,', AWGN, iterations = ',num2str(iterations),', errors = ',num2str(target_block_errors)]);
ylabel('Required E_s/N_0 [dB]');
xlabel('A');
grid on
hold on
drawnow

hMod = NRModulator('Modulation',Modulation);
hDemod = NRDemodulator('Modulation',Modulation);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)');

% Consider each encoded block length in turn
for R_index = 1:length(R)
    
    % Create the plot
    plots(R_index) = plot(nan,'Parent',axes1);
    set(plots(R_index),'XData',A);
    legend(cellstr(num2str(R(1:R_index)', 'R=%0.2f')),'Location','eastoutside');
    
    EsN0s = nan(1,length(A));
    
    % Open a file to save the results into.
    filename = ['results/SNR_vs_A_',num2str(target_BLER),'_',num2str(R(R_index)),'_',num2str(BG),'_',Modulation,'_',num2str(iterations),'_',num2str(target_block_errors),'_',num2str(seed)];
    fid = fopen([filename,'.txt'],'w');
    if fid == -1
        error('Could not open %s.txt',filename);
    end
    
    
    % Consider each information block length in turn
    for A_index = 1:length(A)
        
        found_start = false;
        
        % Skip any combinations of block lengths that are not supported
        try
            % Initialise the BLER and SNR
            BLER=1;
            prev_BLER = nan;
            EsN0 = EsN0_start-EsN0_delta;
            G = round((A(A_index))/R(R_index)/hMod.Q_m)*hMod.Q_m;
            
            
            hEnc = NRLDPCEncoder('A',A(A_index),'BG',BG,'G',G,'Q_m',hMod.Q_m);
            hDec = NRLDPCDecoder('A',A(A_index),'BG',BG,'G',G,'Q_m',hMod.Q_m,'I_HARQ',1,'iterations',iterations);
            
            % Loop over the SNRs
            while BLER > target_BLER
                prev_EsN0 = EsN0;
                EsN0 = EsN0 + EsN0_delta;
                
                hChan.SNR = EsN0;
                hDemod.Variance = 1/10^(EsN0/10);
                
                
                % Start new counters
                block_error_count = 0;
                block_count = 0;
                
                keep_going = true;
                
                % Continue the simulation until enough block errors have been simulated
                while keep_going && block_error_count < target_block_errors
                    
                    
                    
                    a = round(rand(hEnc.A,1));
                    
                    a_hat = [];
                    rv_id_index = 1;
                    reset(hDec); % Reset the incremental redundancy buffer
                    
                    while isempty(a_hat) && rv_id_index <= length(rv_id_sequence)
                        
                        hEnc.rv_id = rv_id_sequence(rv_id_index);
                        hDec.rv_id = rv_id_sequence(rv_id_index);
                        
                        g = step(hEnc,a);
                        tx = step(hMod, g);
                        rx = step(hChan, tx);
                        g_tilde = step(hDemod, rx);
                        a_hat = step(hDec, g_tilde);
                        
                        rv_id_index = rv_id_index + 1;
                    end
                    
                    
                    
                    if found_start == false && ~isequal(a,a_hat)
                        keep_going = false;
                        
                        block_error_count = 1;
                        block_count = 1;
                    else
                        found_start = true;
                        
                        % Determine if we have a block error
                        if ~isequal(a, a_hat)
                            block_error_count = block_error_count+1;
                        end
                        
                        % Accumulate the number of blocks that have been simulated
                        % so far
                        block_count = block_count+1;
                    end
                end
                prev_BLER = BLER;
                BLER = block_error_count/block_count;
            end
        catch ME
            if strcmp(ME.identifier, 'ldpc_3gpp_matlab:UnsupportedParameters')
                warning('ldpc_3gpp_matlab:UnsupportedParameters','The requested combination of parameters is not supported. %s', getReport(ME, 'basic', 'hyperlinks', 'on' ));
                continue
            else
                rethrow(ME);
            end
        end
        % Use interpolation to determine the SNR where the BLER equals the target
        EsN0s(A_index) = interp1(log10([prev_BLER, BLER]),[prev_EsN0,EsN0],log10(target_BLER));
        
        % Plot the SNR vs A results
        set(plots(R_index),'YData',EsN0s);
        
        xlim auto;
        xl = xlim;
        xlim([floor(xl(1)), ceil(xl(2))]);
        
        drawnow;
        
        fprintf(fid,'%d\t%f\n',A(A_index),EsN0s(A_index));
        
        
    end
    
    % Close the file
    fclose(fid);
    
end
