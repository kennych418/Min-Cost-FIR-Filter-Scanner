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
[PM_estimate,f,a,w] = firpmord([Fp Fs],[1 0],[delta_p_decimal delta_s_decimal],F_sample); %Estimate the order
PM_estimate = PM_estimate + 11; %Provide some margin to optimize

%% Create Ideal Filter
FIR_ideal = firpm(PM_estimate,f,a,w);  %Perform Parks-McClellan filter design
[H_ideal, w_ideal] = freqz(FIR_ideal); %Plot ideal filter frequency response

%% Verify Ripple
[p_max_ripple_ideal,s_max_ripple_ideal] = find_max_ripple(H_ideal,w_ideal,wp*pi,ws*pi); %Verify ripple
if p_max_ripple_ideal >= delta_p_decimal || s_max_ripple_ideal >= delta_s_decimal
    fprintf('Ripple does not meet specifications for N = %i\n', PM_estimate);
end

%% Calculating number of bits
margin_p = delta_p_decimal - p_max_ripple_ideal;
margin_s = delta_s_decimal - s_max_ripple_ideal;
if margin_s < margin_p
    n = log(3*sqrt(2)*sqrt(length(FIR_ideal))/margin_s)/log(2);
else
    n = log(3*sqrt(2)*sqrt(length(FIR_ideal))/margin_p)/log(2);
end
n = ceil(n);

%% Quantize coefficients
FIR_quant_coef = round(FIR_ideal,n);
[H_quant_coef,w_quant_coef] = freqz(FIR_quant_coef);

%% Verify Ripple
[p_max_ripple_quant_coef,s_max_ripple_quant_coef] = find_max_ripple(H_quant_coef,w_quant_coef,wp*pi,ws*pi); %Verify ripple
if p_max_ripple_quant_coef >= delta_p_decimal || s_max_ripple_quant_coef >= delta_s_decimal
    fprintf('Ripple does not meet specifications for N = %i\n', PM_estimate);
end

%% Simulate Sine input
F_input = 3000;
sig_size = 20000;
input = 0.5*cos(2*pi*F_input/F_sample*(1:sig_size));

%% Verify SNR


%% Calculate Hardware Cost
cost = 2*(n)^2 + 18*n;
