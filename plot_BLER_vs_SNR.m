function plot_BLER_vs_SNR(A, R, BG, Modulation, rv_id_sequence, iterations, target_block_errors, target_BLER, EsN0_start, EsN0_delta, seed)
% PLOT_BLER_VS_SNR Plots Block Error Rate (BLER) versus Signal to Noise
% Ratio (SNR) for 3GPP New Radio LDPC code.
%   target_block_errors should be an integer scalar. The simulation of each
%   SNR for each coding rate will continue until this number of block
%   errors have been observed. A value of 100 is sufficient to obtain
%   smooth BLER plots for most values of A. Higher values will give
%   smoother plots, at the cost of requiring longer simulations.
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

% Default values
if nargin == 0
    A = 3842;
    R = 1/3;
    BG = 2;
    Modulation = 'QPSK';
    rv_id_sequence = [0];
    iterations = 8;
    target_block_errors = 3;
    target_BLER = 1e-3;
    EsN0_start = 0;
    EsN0_delta = 0.5;
    seed = 0;
end

% Seed the random number generator
rng(seed);


hMod = NRModulator('Modulation',Modulation);
hDemod = NRDemodulator('Modulation',Modulation);
hChan = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (SNR)');

% Consider each base graph in turn
for BG_index = 1:length(BG)
    for R_index = 1:length(R)
        
        % Create a figure to plot the results.
        figure
        axes1 = axes('YScale','log');
        title(['3GPP New Radio LDPC code, R = ',num2str(R(R_index)),', BG',num2str(BG(BG_index)),', ',Modulation,', AWGN, iterations = ',num2str(iterations),', errors = ',num2str(target_block_errors)]);
        ylabel('BLER');
        xlabel('E_s/N_0 [dB]');
        ylim([target_BLER,1]);
        hold on
        drawnow
        
        % Consider each information block length in turn
        for A_index = 1:length(A)
            
            % Create the plot
            plot1 = plot(nan,'Parent',axes1);
            legend(cellstr(num2str(A(1:A_index)', 'A=%d')),'Location','southwest');
            
            % Counters to store the number of bits and errors simulated so far
            block_counts=[];
            block_error_counts=[];
            EsN0s = [];
            
            % Open a file to save the results into.
            filename = ['results/BLER_vs_SNR_',num2str(A(A_index)),'_',num2str(R(R_index)),'_',num2str(BG(BG_index)),'_',Modulation,'_',num2str(iterations),'_',num2str(target_block_errors),'_',num2str(EsN0_start),'_',num2str(seed)];
            fid = fopen([filename,'.txt'],'w');
            if fid == -1
                error('Could not open %s.txt',filename);
            end
            
            % Initialise the BLER and SNR
            BLER = 1;
            EsN0 = EsN0_start;
            
            found_start = false;
            
            % Skip any encoded block lengths that generate errors
            try
                
                G = round((A(A_index))/R(R_index)/hMod.Q_m)*hMod.Q_m;

                
                
                hEnc = NRLDPCEncoder('A', A(A_index),'BG',BG(BG_index),'G',G,'Q_m',hMod.Q_m)
                hDec = NRLDPCDecoder('A', A(A_index),'BG',BG(BG_index),'G',G,'Q_m',hMod.Q_m,'I_HARQ',1,'iterations',iterations);
                                
                
                
                % Loop over the SNRs
                while BLER > target_BLER
                    hChan.SNR = EsN0;
                    hDemod.Variance = 1/10^(EsN0/10);
                    
                    % Start new counters
                    block_counts(end+1) = 0;
                    block_error_counts(end+1) = 0;
                    EsN0s(end+1) = EsN0;
                    
                    keep_going = true;
                    
                    % Continue the simulation until enough block errors have been simulated
                    while keep_going && block_error_counts(end) < target_block_errors
                        
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
                            BLER = 1;
                        else
                            found_start = true;
                            
                            % Determine if we have a block error
                            if ~isequal(a,a_hat)
                                block_error_counts(end) = block_error_counts(end) + 1;
                            end
                            
                            % Accumulate the number of blocks that have been simulated
                            % so far
                            block_counts(end) = block_counts(end) + 1;
                            
                            % Calculate the BLER and save it in the file
                            BLER = block_error_counts(end)/block_counts(end);
                            
                            % Plot the BLER vs SNR results
                            set(plot1,'XData',EsN0s);
                            set(plot1,'YData',block_error_counts./block_counts);
                            drawnow
                        end
                    end
                    
                    if BLER < 1
                        fprintf(fid,'%f\t%e\n',EsN0,BLER);
                    end
                    
                    % Update the SNR, ready for the next loop
                    EsN0 = EsN0 + EsN0_delta;
                    
                end
            catch ME
                if strcmp(ME.identifier, 'ldpc_3gpp_matlab:UnsupportedParameters')
                    warning('ldpc_3gpp_matlab:UnsupportedParameters','The requested combination of parameters is not supported. %s', getReport(ME, 'basic', 'hyperlinks', 'on' ));
                    continue
                else
                    rethrow(ME);
                end
            end
            
            % Close the file
            fclose(fid);
        end
    end
end




