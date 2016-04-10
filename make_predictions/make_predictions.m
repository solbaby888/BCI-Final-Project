function [predicted_dg] = make_predictions(test_ecog)

% Inputs: test_ecog - 3 x 1 cell array containing ECoG for each subject, where test_ecog{i} 
% to the ECoG for subject i. Each cell element contains a N x M testing ECoG,
% where N is the number of samples and M is the number of EEG channels.
% Outputs: predicted_dg - 3 x 1 cell array, where predicted_dg{i} contains the 
% data_glove prediction for subject i, which is an N x 5 matrix (for
% fingers 1:5)

% Run time: The script has to run less than 1 hour. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following is a sample script.

% Load Model. The variable models contains weights for each person,
% each channel, and each finger.

load team_awesome_model.mat 


% Predict using linear predictor for each subject
%create cell array with one element for each subject
predicted_dg = cell(3,1);

%for each subject
for subj = 1:3 
    
    %get the testing ecog
    testset = test_ecog{subj}; 
    
    %initialize the predicted dataglove matrix
    yhat = zeros(size(testset,1),5);
    
    %for each finger
    for i = 1:5 
        
        %predict dg based on ECOG for each finger
        yhat(:,i) = testset*models{subj}(:,i);
        
    end
    predicted_dg{subj} = yhat;
end

