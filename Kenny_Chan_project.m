%% Desired Filter Specs
F_sample = 1000000; %Input sample rate = 1MHz
Fp = 100000; %Pass-band edge = 100kHz
Fs = 125000; %Stop-band edge = 125kHz
A = 60; %Minimum stop-band attenuation = 60dB
delta_p = 0.1; %Maximum passband ripple = 0.1dB
group_delay = 0.5; %Pass-band group delay variation = +-0.5samples
SNR = 72; %Minumum output SNR = 72dB
 
%% Convert units
wp = 2*Fp/F_sample; %normalized to pi
ws = 2*Fs/F_sample; %normalized to pi
delta_s = 1/A;
delta_p_decimal = -10^(-0.05*delta_p) + 1; %convert to decimal
delta_s_decimal = 1/(10^((1/20)*A));  

%% Generate an estimate filter order
EXTENDEDRANGE = 6;
[PM_estimate,f,a,w] = firpmord([Fp Fs],[1 0],[delta_p_decimal delta_s_decimal],F_sample); %Estimate the order
N = PM_estimate + EXTENDEDRANGE; %Provide some margin to optimize cost

%% Begin optimizing for lowest cost that meets specs
RippleMet = true;
RipplePostQuantMet = true;
SNRMet = true;
cost_index = [];
order_index = [];
bits_index = [];
SNR_index = [];
recorded_SNR = 0;

while(RippleMet && RipplePostQuantMet && SNRMet)
    %% Decrement
    N = N - 1;
    
    %% Create Ideal Filter
    FIR_ideal = firpm(N,f,a,w);  %Perform Parks-McClellan filter design
    [H_ideal, w_ideal] = freqz(FIR_ideal); %Plot ideal filter frequency response

    %% Verify Ripple before Quantization
    [p_max_ripple_ideal,s_max_ripple_ideal] = find_max_ripple(H_ideal,w_ideal,wp*pi,ws*pi); %Verify ripple
    if p_max_ripple_ideal >= delta_p_decimal || s_max_ripple_ideal >= delta_s_decimal
        fprintf('Ripple does not meet specifications for N = %i\n', N);
        RippleMet = false;
        break;
    end

    %% Estimating number of bits
    margin_p = delta_p_decimal - p_max_ripple_ideal;
    margin_s = delta_s_decimal - s_max_ripple_ideal;
    if margin_s < margin_p
        n = log(3*sqrt(2)*sqrt(length(FIR_ideal))/margin_s)/log(2);
    else
        n = log(3*sqrt(2)*sqrt(length(FIR_ideal))/margin_p)/log(2);
    end
    n = ceil(n) + 1;

    %% Finding the minimal number of bits to meet SNR and ripple
    num_loops = 0;
    recorded_SNR = 0;
    while(RipplePostQuantMet && SNRMet)
        num_loops = num_loops + 1;
        n = n - 1;
        %% Quantize coefficients
        FIR_quant_coef = round(FIR_ideal,n);
        [H_quant_coef,w_quant_coef] = freqz(FIR_quant_coef);

        %% Verify Ripple after Quantization
        [p_max_ripple_quant_coef,s_max_ripple_quant_coef] = find_max_ripple(H_quant_coef,w_quant_coef,wp*pi,ws*pi); %Verify ripple
        if p_max_ripple_quant_coef >= delta_p_decimal || s_max_ripple_quant_coef >= delta_s_decimal
            RipplePostQuantMet = false;
            break;
        end

        %% Simulate Sine input and Verify SNR
        F_input = 3000;
        sig_size = 20000;
        input = 0.5*cos(2*pi*F_input/F_sample*(1:sig_size));
        input = [input,zeros(1,length(FIR_quant_coef)-1)];  %Pad the input with 0's based on FIR_quant_coef size
        y = zeros(1,sig_size);
        for i=1:sig_size
            y(i) = generate(input((sig_size-i+1):(length(FIR_quant_coef)+sig_size-i)),FIR_quant_coef,n);
        end
        r = snr(y);
        if r < SNR
            SNRMet = false;
            break;
        else
            recorded_SNR = r;
        end
    end
    %% Record results and Reset
    if (num_loops > 1) %only if the filter could meet specs at least once
        RipplePostQuantMet = true;
        SNRMet = true;
        n = n + 1; %Keep the number of bits that meets specs, not the one that broke them
        order_index = [order_index, N];
        bits_index = [bits_index, n];
        SNR_index = [SNR_index, recorded_SNR];
        if(rem(N,2)==0) %Cost for even filter length
            cost_index = [cost_index, (2*(n)^2 + 23*n)*N/2 - 9*n];
        else %Cost for odd filter length
            cost_index = [cost_index, (2*(n)^2 + 23*n)*(N-1)/2 + 2*n^2 + 5*n];
        end
    end
    
end

%% Find the optimal order and cost
min_index = 0;
min_bits = 0;
min_cost = inf;
min_SNR = 0;
for i = 1:length(cost_index)
    if(cost_index(i) < min_cost)
        min_index = order_index(i);
        min_bits = bits_index(i);
        min_cost = cost_index(i);
        min_SNR = SNR_index(i);
    end
end
fprintf('Optimal filter has order %i, %i+1 bits, %i SNR, and %i hardware cost.\n', min_index-1,min_bits,min_SNR,min_cost);

%% Plot scan results
%figure(1);
%plot(order_index, cost_index);
%title("Cost vs Filter Length");

%figure(2);
%plot(order_index, SNR_index);
%title("SNR vs Filter Length");

%figure(3);
%plot(order_index, bits_index);
%title("Number of Bits vs Filter Length");

%% Get optimal filter's performance and graphs
optimal_f_passband = 0;
optimal_f_stopband = pi;
optimal_f_attenuation = 0;
optimal_f_passripple = 0;

FIR_ideal = firpm(min_index,f,a,w);
[H_ideal, w_ideal] = freqz(FIR_ideal); %Plot ideal filter frequency response
figure(1)
title("Unquantized Filter Frequency Reponse")
freqz(FIR_ideal);
grid
for i = 1:length(w_ideal) %Scan for passband edge
    if( abs(H_ideal(i)) < 1-delta_p_decimal)
        optimal_f_passband = w_ideal(i)/pi;
        break;
    end
end
for i = 1:length(w_ideal) %Scan for stopband edge
    if( abs(H_ideal(length(w_ideal)+1-i)) > delta_s_decimal)
        optimal_f_stopband = w_ideal(length(w_ideal)+1-i)/pi;
        break;
    end
end
[p_max_ripple_ideal,s_max_ripple_ideal] = find_max_ripple(H_ideal,w_ideal,wp*pi,ws*pi);
optimal_f_passband_ripple = 20*log10(p_max_ripple_ideal+1); %Convert decimal ripple to dB
optimal_f_A = -20*log(s_max_ripple_ideal)/log(10);   %Convert decimal ripple to attenuation
optimal_f_grpdelay = grpdelay(FIR_ideal);
figure(2)
title("Unquantized Filter Impulse Reponse")
impz(FIR_ideal);
grid
figure(3)
title("Unquantized Filter Pole-Zero Plot")
zplane(FIR_ideal);
grid

%% Get quantized filter's performance and graphs
quant_f_passband = 0;
quant_f_stopband = pi;
quant_f_attenuation = 0;
quant_f_passripple = 0;

FIR_quant_coef = round(FIR_ideal,n);
[H_quant_coef,w_quant_coef] = freqz(FIR_quant_coef);
figure(4)
title("Quantized Filter Frequency Reponse")
freqz(FIR_quant_coef);
grid
for i = 1:length(w_quant_coef) %Scan for passband edge
    if( abs(H_quant_coef(i)) < 1-delta_p_decimal)
        quant_f_passband = w_quant_coef(i)/pi;
        break;
    end
end
for i = 1:length(w_quant_coef) %Scan for stopband edge
    if( abs(H_quant_coef(length(w_quant_coef)+1-i)) > delta_s_decimal)
        quant_f_stopband = w_quant_coef(length(w_quant_coef)+1-i)/pi;
        break;
    end
end
[p_max_ripple_ideal,s_max_ripple_ideal] = find_max_ripple(H_quant_coef,w_quant_coef,wp*pi,ws*pi);
quant_f_passband_ripple = 20*log10(p_max_ripple_ideal+1); %Convert decimal ripple to dB
quant_f_A = -20*log(s_max_ripple_ideal)/log(10);   %Convert decimal ripple to attenuation
quant_f_grpdelay = grpdelay(FIR_quant_coef);
figure(5)
title("Quantized Filter Impulse Reponse")
impz(FIR_quant_coef);
grid
figure(6)
title("Quantized Filter Pole-Zero Plot")
zplane(FIR_quant_coef);
grid





%% Defining helper functions
function y = generate(x,h,n) %Calculates summation(x[n]h[n]) and quantizes to match filter architecture
    i = 1;
    y = 0;
    N = length(x);
    if(rem(N,2)==0) %For even filter length    
        while i <= N/2
            y = y + quant(h(i)*(x(i)+x(N-i)),1/2^n);
            i = i + 1;
        end
    else    %For odd filter length
        while i <= (N-1)/2
            y = y + quant(h(i)*(x(i)+x(N-i)),1/2^n);
            i = i + 1;
        end
        y = y + h((N-1)/2 + 1)*x((N-1)/2 + 1);
    end
end
function [p_max_ripple,s_max_ripple] = find_max_ripple(H,w,wp,ws) %Finds the maximum ripple in pass and stop band
    p_max_ripple = 0;
    s_max_ripple = 0;
    i=1;
    %Parse through H from 0 to passband, calculate ripples, and save the max
    while(w(i) < wp)
        ripple = abs(1-abs(H(i)));
        if ripple > p_max_ripple
            p_max_ripple = ripple;
        end
        i = i + 1;
    end
    %Skip the transition band
    while(w(i) < ws)
        i = i + 1;
    end
    %Parse through H from stopband to pi, calculate ripples, and save the max
    while(i < length(w))
        ripple = abs(H(i));
        if ripple > s_max_ripple
            s_max_ripple = ripple;
        end
        i = i + 1;
    end
    
end
