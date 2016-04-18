function feat = FFT_featFn(signal, fs);

Y           = fft(signal);
L           = length(signal);
f           = fs *(0:(L/2))/L;
P2          = abs(Y/L);
P1          = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

tmp = abs(f-5);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(f-15);
[idx2 idx2] = min(tmp); %index of closest value
FFT_5_15    = mean(P1(idx1:idx2));

tmp = abs(f-20);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(f-25);
[idx2 idx2] = min(tmp); %index of closest value
FFT_20_25   = mean(P1(idx1:idx2));

tmp = abs(f-75);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(f-115);
[idx2 idx2] = min(tmp); %index of closest value
FFT_75_115  = mean(P1(idx1:idx2));

tmp = abs(f-125);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(f-160);
[idx2 idx2] = min(tmp); %index of closest value
FFT_125_160 = mean(P1(idx1:idx2));

tmp = abs(f-160);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(f-175);
[idx2 idx2] = min(tmp); %index of closest value
FFT_160_175 = mean(P1(idx1:idx2));

feat = [FFT_5_15 FFT_20_25 FFT_75_115 FFT_125_160 FFT_160_175];
