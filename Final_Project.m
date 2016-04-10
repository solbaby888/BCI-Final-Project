% Final Project

%% Import Data

% Subject 1
sesh_sub1_1 = IEEGSession('I521_A0012_D001', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub1_2 = IEEGSession('I521_A0012_D002', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub1_3 = IEEGSession('I521_A0012_D003', 'solbaby888', 'sol_ieeglogin.bin');
nr_1        = ceil((sesh_sub1_1.data(1).rawChannels(11).get_tsdetails.getEndTime)/...
                1e6*sesh_sub1_1.data(1).sampleRate);
% Subject 2
sesh_sub2_1 = IEEGSession('I521_A0013_D001', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub2_2 = IEEGSession('I521_A0013_D002', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub2_3 = IEEGSession('I521_A0013_D003', 'solbaby888', 'sol_ieeglogin.bin');
nr_2        = ceil((sesh_sub2_1.data(1).rawChannels(11).get_tsdetails.getEndTime)/...
                1e6*sesh_sub2_1.data(1).sampleRate);

% Subject 3
sesh_sub3_1 = IEEGSession('I521_A0014_D001', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub3_2 = IEEGSession('I521_A0014_D002', 'solbaby888', 'sol_ieeglogin.bin');
sesh_sub3_3 = IEEGSession('I521_A0014_D003', 'solbaby888', 'sol_ieeglogin.bin');
nr_3        = ceil((sesh_sub2_3data(1).rawChannels(11).get_tsdetails.getEndTime)/...
                1e6*sesh_sub2_1.data(1).sampleRate);

%{
    TL_Comment: Loads pretty fast

%}

%% Check for 60 Hz noise

Y = fft(Chan_11);
L = length(Chan_11);
f = fs*(0:(L/2))/L;
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
figure;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')