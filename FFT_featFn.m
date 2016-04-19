function feat = FFT_featFn(signal, fs, window, overlap);

for i = 1:length(signal);
[s w t] = spectrogram(signal{i}, window*fs, overlap*fs, [], fs, 'yaxis');


tmp = abs(w-5);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(w-15);
[idx2 idx2] = min(tmp); %index of closest value
FFT_5_15    = mean(s(idx1:idx2, :),1);

tmp = abs(w-20);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(w-25);
[idx2 idx2] = min(tmp); %index of closest value
FFT_20_25   = mean(s(idx1:idx2, :),1);
 
tmp = abs(w-75);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(w-115);
[idx2 idx2] = min(tmp); %index of closest value
FFT_75_115  = mean(s(idx1:idx2, :),1);

tmp = abs(w-125);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(w-160);
[idx2 idx2] = min(tmp); %index of closest value
FFT_125_160 = mean(s(idx1:idx2, :),1);

tmp = abs(w-160);
[idx1 idx1] = min(tmp); %index of closest value
tmp = abs(w-175);
[idx2 idx2] = min(tmp); %index of closest value
FFT_160_175 = mean(s(idx1:idx2, :),1);

feat{:,i} = [FFT_5_15' FFT_20_25' FFT_75_115' FFT_125_160' FFT_160_175'];
end
