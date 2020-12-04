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
delta_p_decimal = -10^(-0.05*delta_p) + 1; %convert to decimal
delta_s_decimal = 1/(10^((1/20)*A));  
%k = delta_p/delta_s;

%% Create and analyze ideal filter
FIR_ideal = designfilt('lowpassfir','PassbandFrequency',wp,'StopbandFrequency',ws,'PassbandRipple',delta_p,'StopbandAttenuation',60,'DesignMethod','equiripple');
specs_ideal = info(FIR_ideal);
[H_ideal, w_ideal] = freqz(FIR_ideal); %Plot ideal filter frequency response
[p_max_ripple_ideal,s_max_ripple_ideal] = find_max_ripple(H_ideal,w_ideal,wp*pi,ws*pi); %Verify ripple
zplane(FIR_ideal); %Plot poles and zeros

%% Calculating number of bits
margin_p = delta_p_decimal - p_max_ripple_ideal;
margin_s = delta_s_decimal - s_max_ripple_ideal;
if margin_s < margin_p
    n = log(3*sqrt(2)*sqrt(length(FIR_ideal.Coefficients))/margin_s)/log(2);
else
    n = log(3*sqrt(2)*sqrt(length(FIR_ideal.Coefficients))/margin_p)/log(2);
end
n = ceil(n);

%% Quantize coefficients
FIR_quant_coef = FIR_ideal.Coefficients;
FIR_quant_coef = round(FIR_quant_coef,n);

[H_quant_coef,w_quant_coef] = freqz(FIR_quant_coef);
[p_max_ripple_quant_coef,s_max_ripple_quant_coef] = find_max_ripple(H_quant_coef,w_quant_coef,wp*pi,ws*pi); %Verify ripple
zplane(FIR_quant_coef); %Plot poles and zeros
plot(FIR_quant_coef);   %plot impulse response

%% Calculate Hardware Cost
cost = 2*(n)^2 + 18*n;
%% Estimate filter order
%Kaiser_estimate = ceil((-20*log10(sqrt(delta_s*delta_p))-13)/(14.6*abs(ws-wp)/(2*pi)));
%Bellanger_estimate = ceil(-2*log10(10*delta_s*delta_p)/(3*abs(ws-wp)/(2*pi))) - 1;

%{
%% Optimization
[PM_estimate,f,a,w] = firpmord([Fp Fs],[1 0],[delta_p delta_s],F_sample);
%f = [0,wp,ws,1]; %Define the stopband and passband
%a = [1,1,0,0]; %Define the desired amplitude in each band
%w = [1/k, 1]; %Define the weights in each band
%if rem(Kaiser_estimate,2) == 0 %Make sure estimate is odd for type I
%    Kaiser_estimate = Kaiser_estimate+1;
%end
[FIR, err]= firpm(PM_estimate+4,f,a,w);  %Perform Parks-McClellan filter design
%FIR = round(FIR, 5);   %Round filter coefficients to 4 decimal places (4+1 bits)
%while err < delta_p %Find lowest order filter that meets specs
%    Kaiser_estimate = Kaiser_estimate - 1;
%    [FIR, err]= firpm(Kaiser_estimate,f,a);
%end
%fprintf('%i\n',Kaiser_estimate + 1); %Best filter order
figure(1);
freqz(FIR);
%}