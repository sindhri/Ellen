data_path = 'matlab_zf/sample_uploaded_data/';
addpath(genpath('matlab_zf/code/sleepwake2/'));
varargin= {'summary','$SUMMARY','plot',true,'overwrite',false};
[out,combine] = sleepanalysis_batch_Jia(data_path,varargin{:});

%perl does not like space. Don't know how it was taken care of
%farnam:
%1. restwake_bott to pass the parameters
%2. restwake_slurm
%3. sleepanalysis_batch(_jia)
%4. sleepanalysis(2), for individual data and plots
%5. sleepwake for summary
%6. I don't think sleepanalysis_combine was used in the pipeline
%it's part of sleepanalysis_batch