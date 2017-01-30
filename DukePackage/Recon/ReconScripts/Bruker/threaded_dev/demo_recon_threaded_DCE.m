%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A demonstration of basic Dynamic Contrast Enhanced (DCE) reconstruction
% using radial data using a simple sliding window reconstruction.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% code_path='/Users/james/Desktop/DCE_proto';
% run([code_path '/Duke-CIVM-MRI-Tools/setup.m']);
% run([code_path '/GE-MRI-Tools/setup.m']);
% run([code_path '/Non-Cartesian-Reconstruction/setup.m'])
% u_dir='/Users/james/';

mex -largeArrayDims mex_thread_sparseMultiply.c;
mex -largeArrayDims mex_thread_calcSparseGridMatrix.c

%% Slow (but decent) Reconstruction parameters
overgridding = 3*[1 1 1]; % Use at least 2. Turn "crop" off on reconObject, then increase this until you dont see signal at the edge of FOV, then turn "crop" back on and use that oversampling
kern_sigma = 0.3*[1 1 1];  % This is a key parameter that tradesoff SNR and resolution (making sharpness smaller will blurr the object, but increase SNR and vice versa)
kern_extent = 6*[1 1 1].*kern_sigma; % 9 is a good value
nThreads = 12;
verbose = 0;
nDcfIter = 15;
cropVolume = 1;
deapodize = 1;
saveFullVol = 1;


%% Read the header file to get scan info
% This could be done nicer by reading the header, etc. - I was lazy and hard-coded
nPts = 64;
nCoils = 2;
nAcq = 4*3;
samplesPerAcq = 25735;
mat_size = 2*nPts*[1 1 1];
overgrid_mat_size = ceil(mat_size.*overgridding);
overgridding = overgrid_mat_size./mat_size;

%% Sliding window parameters
keysPerWindow = 1980*13*nPts; % 3 keys of data (I like 10-15 for this dataset)
windowStep = round(keysPerWindow/13); % 1/25th of window

load('data.mat');
load('traj.mat');

%% Perform sliding window Reconstruction
totalSamples = nPts*samplesPerAcq*nAcq; % all keys, acqs
windowStartIdxs = 1:windowStep:(totalSamples-keysPerWindow-1);
nWindows = length(windowStartIdxs);

tmpVol = zeros(overgrid_mat_size);
tmpData = zeros([1 keysPerWindow]);
tmpComplexVol = complex(zeros(overgrid_mat_size));
tmpComplexData = complex(zeros([1 keysPerWindow]));
sosComplexVol = complex(zeros(overgrid_mat_size));
dcf = complex(zeros([1 keysPerWindow]));
if(cropVolume)
    cropSosVol = complex(zeros(mat_size));
end

if(saveFullVol)
    warning('About to make a really big matrix!');
    if(cropVolume)
        fullVolume = complex(zeros([mat_size nWindows]));
    else
        fullVolume = complex(zeros([overgrid_mat_size nWindows]));
    end
end

tic;
startingTime = toc
for iWin = 1:nWindows
    disp(['Reconstructing ' num2str(iWin) '/' num2str(nWindows) ' traj subsets']);
    
    % Make a trajectory for each window
    windowKeys = windowStartIdxs(iWin)+[1:keysPerWindow]-1;
    windowTraj = squeeze(traj(:,windowKeys));
    
    % Construct system model for this trajectory

    A = SparseGridder(windowTraj,mat_size,kern_sigma,kern_extent,overgridding,nThreads);
       
       
    % Option 1 - dcf
    disp('   Computing density compensation weights');
    tmpVol = ones(A.overgridSize);
    tmpData = A.ungrid(tmpVol);
    dcf = 1./tmpData;
    for iDcfIter=1:nDcfIter
        disp(['      DCF iter ' num2str(iDcfIter) '/' num2str(nDcfIter)])
        A.grid(dcf,tmpVol);
        A.ungrid(tmpVol,tmpData);
        dcf = dcf./tmpData;
    end
    
    % Create a data matrix of all repetitions of this trajectory
    windowData = data(:,windowKeys);
    
    % Grid data from all coils and compute SOS
    sosComplexVol = complex(zeros(overgrid_mat_size)); % reset to zero
    for iCoil = 1:nCoils
        % Apply dcf
        windowData(iCoil,:)  = windowData(iCoil,:).*dcf;
        
        % Calculate gridded kspace
        tic;
        A.grid(windowData(iCoil,:),tmpComplexVol);
        gridTime = toc
        
        % Reconstruct image domain with IFFT
        tmpComplexVol = ifftshift(ifftn(tmpComplexVol));
        
        % Accumulate SOS
        sosComplexVol = sosComplexVol + tmpComplexVol.^2;
        disp(['      Finished Coil ' num2str(iCoil) '/' num2str(nCoils)]);
    end
    
    % Finish SOS
    sosComplexVol = sqrt(sosComplexVol);
    
    % Compute deapodization volume for this traj
    if(deapodize)
        disp('   Deapodizing...');
        tmpData = ~any(windowTraj,1).*dcf;
        A.grid(tmpData,tmpVol);
        tmpVol = ifftshift(ifftn(tmpVol));
        sosComplexVol = sosComplexVol./tmpVol;
    end
    
    % Crop volume
    if(cropVolume)
        disp('   Cropping volume...');
        cropSosVol = subvolume(sosComplexVol,...
            [round([0.5*(A.overgridSize-A.gridSize)+1]); ...
            round([0.5*(A.overgridSize+A.gridSize)])]);
    end
    
    % Save/store this time point
    if(~saveFullVol)
        disp('   Saving Data');
        if(cropVolume)
            nii = make_nii(abs(cropSosVol));
        else
            nii = make_nii(abs(sosComplexVol));
        end
        save_nii(nii,['reconTime_' num2str(iWin) '.nii']);
    else
        if(cropVolume)
            fullVolume(:,:,:,iWin) = cropSosVol;
        else
            fullVolume(:,:,:,iWin) = sosComplexVol;
        end
    end
    
    disp('Completed another window')
    % Compute remaining time
    currentTime = toc;
    timeSoFar = currentTime - startingTime;
    timePerIter = timeSoFar/iWin;
    totalEstTime = timePerIter*nWindows;
    remainingTime = totalEstTime - timeSoFar;
    
    % Show some progress
    disp(['Completed ' num2str(iWin) '/' num2str(nWindows) ' traj subsets, est total time =~' num2str(totalEstTime) ' sec (' num2str(timeSoFar) ' sec so far, ~' num2str(remainingTime) ' sec remaining)']);
end
endingTime = toc
if(saveFullVol)
    nii = make_nii(abs(fullVolume));
    save_nii(nii,['reconVol.nii']);
end
reconTime = toc
