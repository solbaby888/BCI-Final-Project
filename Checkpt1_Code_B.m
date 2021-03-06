    % Checkpoint 1

%{ 
    % Working only with Subject 1 but with all the channels of subject 1
    
    Table of Contents (By Section):
    % Import Data
    % Apply 60 Hz filter
    % Checkpt 1

%}

%% Import Data
%{
    Each subject has ECoG data set, a glove data set, and testing data set.
    
%}

% Subject 1
sesh_sub1_1 = IEEGSession('I521_A0012_D001', 'karanm', 'kar_ieeglogin.bin');
sesh_sub1_2 = IEEGSession('I521_A0012_D002', 'karanm', 'kar_ieeglogin.bin');
sesh_sub1_3 = IEEGSession('I521_A0012_D003', 'karanm', 'kar_ieeglogin.bin');
nr_1        = ceil((sesh_sub1_1.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
                1e6*sesh_sub1_1.data(1).sampleRate);

% Testing with all channels of Subject 1
no_of_channels_ECOG = size(sesh_sub1_1.data(1).rawChannels, 2);
for i = 1:no_of_channels_ECOG;
    ECoG_Sub1_Chan{:, i} = sesh_sub1_1.data(1).getvalues(1:nr_1, i);   
end;

% Sampling Rate
fs_Sub1 = sesh_sub1_1.data(1).sampleRate;

%% Apply 60 Hz filter and Bandpass filter
%{
    Will eventually have to save as a function to run on all channels
%}

% 60 Hz Notch Filter
f0        = 60;         % notch frequency
fn        = fs_Sub1/2;  % Nyquist frequency
freqRatio = f0/fn;      % ratio of notch freq. to Nyquist freq.

notchWidth = 0.1;       % width of the notch

%#Compute zeros
zeross = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )]; %two 's' on purpose

%#Compute poles
poles = (1-notchWidth) * zeross;

figure;
zplane(zeross.', poles.');

b = poly( zeross ); % Get moving average filter coefficients
a = poly( poles );  % Get autoregressive filter coefficients

%#filter signal
for i = 1:no_of_channels_ECOG;
    ECoG_Sub1_Chan_filt{i} = filtfilt(b, a, ECoG_Sub1_Chan{i});  
end;

% Bandpass filter
A_stop1 = 60;		% Attenuation in the first stopband = 60 dB
F_stop1 = 0.1;		% Edge of the stopband = 8400 Hz
F_pass1 = 0.15;     % Edge of the passband = 10800 Hz
F_pass2 = 200;      % Closing edge of the passband = 15600 Hz
F_stop2 = 201;      % Edge of the second stopband = 18000 Hz
N_order = 100;

BandPassSpecObj = ...
   fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2', ...
		N_order, F_stop1, F_pass1, F_pass2, F_stop2, fs_Sub1)
    
   Hd = design(BandPassSpecObj,'equiripple');
   
for i = 1:no_of_channels_ECOG;
   ECoG_Sub1_Chan1_filt_2{i} = filtfilt(Hd.Numerator,1,ECoG_Sub1_Chan_filt{i});
end;
%{
    TL_comment: Not sure if this is the best way to filter the data. Should check if
    filtering between .15 and 200 is right. Need to check if we should only
    use channels 48, 63, 47, 64, 61
%}

%% Checkpt 1
%{ 
    window of 100 ms; overlap of 50 ms

    Some comments: I accidently decided to calculate some other features.
    We can remove them for now for the first checkpoint. Stupid extra work.
    
%}

% ====== Function handle to calculate number of windows ================= %
    NumWins = @(xLen, fs, winLen, winDisp) 1+floor((xLen-fs*winLen)/(winDisp*fs));  
    % winLen and winDisp are in secs, xLen are in samples
% ======================================================================= %

% Number of Windows for Subj1
window_time = 100*10^-3;                                 % (secs)
overlap     = 50*10^-3;                                  % (secs)
windowLen   = window_time * fs_Sub1;                     % number of pts per window 
L           = length(ECoG_Sub1_Chan1_filt_2{1});
NoW         = NumWins(L, fs_Sub1, window_time, overlap); 

% Calculate frequency domain avg power features over windows
%# Calculates spectrogram and respective averages for assignment-specified frequencies
Freq_domain_Feat     = FFT_featFn(ECoG_Sub1_Chan1_filt_2, fs_Sub1, window_time, overlap); % windows X 5 Matrix

% Calculate time domain avg voltage 
kernnel = repmat(1/(window_time*fs_Sub1), 1, window_time*fs_Sub1);

for i=1:no_of_channels_ECOG;
    time_avg_volt_temp = conv(ECoG_Sub1_Chan1_filt_2{1,i}, kernnel, 'valid');
    time_avg_volt(:, i) = time_avg_volt_temp(1:overlap*fs_Sub1:end);  
end
% Recall: overlap is 50 ms. 

% Feature Matrix
for i=1:no_of_channels_ECOG
    Feat_Mat{i} = [time_avg_volt(:,i) Freq_domain_Feat{1,i}];    % time window X no. of feature
end
% Downsampling dataglove traces                 
%# Cell array glove positions (1:5)                   
nr                  = ceil((sesh_sub1_2.data(1).rawChannels(1).get_tsdetails.getEndTime)/...
                        1e6*sesh_sub1_2.data(1).sampleRate);
for i = 1:5;
    Glovedata_Sub1{i} = sesh_sub1_2.data(1).getvalues(1:nr, i);
end;

SR_dataglove        = sesh_sub1_2.data(1).sampleRate;
for i = 1:5;
    Glovedata_Sub1_ds{i}  =  round(decimate(Glovedata_Sub1{i},  50)); 
end;

% WARNING: may want to remove rounding. Ignoring 1st point of decimate

% /(SR_dataglove*10^-3)
% Linear Regression Prediction
numoffeat       = 6;
numofprev_win   = 3;
noofchan        = size(sesh_sub1_1.data.rawChannels, 2);     % number of channels
n_of_R          = NoW - numofprev_win;                       % number of windows for regression
p_of_R          = noofchan * numoffeat * numofprev_win;            
R_mat           = zeros(n_of_R, p_of_R);                     % six features per window
curr_pt         = 3;

for i = 1:n_of_R;
        curr_pt = 1 + curr_pt;
        for j = 1:noofchan;
            R_idx1 = (j-1)* numoffeat * numofprev_win + 1;
            R_idx2 = R_idx1 + numoffeat * numofprev_win - 1;
            R_mat(i, R_idx1:R_idx2) = reshape(Feat_Mat{j}(curr_pt - 3:curr_pt - 1, :)', [1, 18]);
        end;
end; 

% Adding the first columns of ones
R_ones     = ones(length(R_mat),1);
R_mat      = [R_ones R_mat];
R_pad      = zeros(numofprev_win,min(size(R_mat))-1);

% Pad
%# position 2

for j=1:no_of_channels_ECOG
        % Prev window (-1) 
         R_idx1 = ((j-1)* numoffeat * numofprev_win + 13);
         R_idx2 = R_idx1 + numoffeat - 1;
         R_pad(2, R_idx1:R_idx2) = Feat_Mat{j}(1,:);
end

%# position 3

for j=1:no_of_channels_ECOG
        % Prev window (-2) 
         R_idx1 = (j-1)* numoffeat * numofprev_win + 7;
         R_idx2 = R_idx1 +(numoffeat * (numofprev_win-1))-1;
         R_pad(3, R_idx1:R_idx2) = reshape(Feat_Mat{j}(1:2, :)', [1, 12]);
end
 R_pad_ones        = ones(3,1);
 R_pad             = [R_pad_ones R_pad];
 R_matrix          = zeros(NoW,min(size(R_mat)));
 R_matrix(1:3,:)   = R_pad;
 R_matrix(4:end,:) = R_mat;

% Weights and prediction for each finger

for i = 1:5
    f{i}         = mldivide(R_matrix'*R_matrix, R_matrix' * Glovedata_Sub1_ds{i}(2:end));
    pred{i}      = R_matrix*f{i};
end

% spline interpolation
%# the interpolation isnt working yet..needs work
x  = (overlap*fs_Sub1:overlap*fs_Sub1:L-overlap*fs_Sub1);
xx = (1:1:L);

pred_upsample = [];

for i = 1:5
    pred_upsample(:,i) = spline(x,pred{i}, xx);
end
   
for i = 1:5;

   Glovedata_mat_Sub1(:,i) = Glovedata_Sub1{i};

end;

Training_Correlation = corr(Glovedata_mat_Sub1(:,1), round(pred_upsample(:,1)));
    
% need to find average testing corellation over 4 fingers per subject
