function sliding_window_recon(data_buffer,opt_struct,data_in,data_work,data_out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A demonstration of basic Dynamic Contrast Enhanced (DCE) reconstruction
% using radial data using a simple sliding window reconstruction.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
% Modified into a function by James Cook.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%
% code_path='/Users/james/Desktop/DCE_proto';
% run([code_path '/Duke-CIVM-MRI-Tools/setup.m']);
% run([code_path '/GE-MRI-Tools/setup.m']);
% run([code_path '/Non-Cartesian-Reconstruction/setup.m'])
% u_dir='/Users/james/';
if ~exist('data_buffer','var')
    error('MUST HAVE A DATABUFFER TO CATCH OUTPUT');
end
if ~exist('opt_struct','var')
    error('must have options for manual mode, specify at least options.dataFile')
end
if ~isfield(opt_struct,'radial_mode')
    opt_struct.radial_mode='good';
end

if strcmp(opt_struct.radial_mode,'fast')
    %% Fast (for quick quality control) Reconstruction parameters
    overgridding = [1 1 1]; % Use at least 2. Turn "crop" off on reconObject, then increase this until you dont see signal at the edge of FOV, then turn "crop" back on and use that oversampling
    sharpness = 0.21*[1 1 1]; 
    extent = 6*sharpness;
    verbose = 0;
    nPipeIter = 2; % Use 3-15 for option 1, lots more (~75) for iterative recon. This should scale recon time linearly(ish)
    cropVolume = 1;
    saveFullVol = 1;
    deapodize = 1;
elseif strcmp(opt_struct.radial_mode,'good')
    %% Slow (but decent) Reconstruction parameters
    overgridding = [2 2 2]; % Use at least 2. Turn "crop" off on reconObject, then increase this until you dont see signal at the edge of FOV, then turn "crop" back on and use that oversampling
    sharpness = 0.21*[1 1 1];  % This is a key parameter that tradesoff SNR and resolution (making sharpness smaller will blurr the object, but increase SNR and vice versa)
    if isfield(opt_struct,'sharpness')
        sharpness=opt_struct.sharpness;
    end
    extent = 6*sharpness; % 6 is a good value
    verbose = 0;
    nPipeIter = 3; % Use 3-15 for option 1, lots more (~75) for iterative recon. This should scale recon time linearly(ish)
    cropVolume = 1;
    saveFullVol= 1;
    deapodize = 1;
else
    db_inplace('scott_grid');
end

%% Find directory containing all information
% % reconDir = [u_dir '/Desktop/B03180/'];

%% Read the header file to get scan info
nKeys = 13;% a fully sampled acq has this many keys. Constant for now.
if isfield(data_buffer.headfile,'B_vol_type') ... 
        && isempty(regexpi(data_buffer.headfile.B_vol_type,'.*keyhole.*'))
    nKeys=1;
end
if ~exist ('data_in','var')
    nPts=64;
    nCoils=4;
    nRaysPerKey=1980;
    nAcq=11;
else
    % This could be done nicer by reading the header, etc. - I was lazy and hard-coded
    nPts = data_in.ray_length;%64;
    if isfield(data_in,'ramp_points')
        nPts=nPts+data_in.ramp_points;
    end
    nCoils = data_in.ds.Sub('c');%2;%4;
    nRaysPerKey = data_in.rays_per_block;%1980;
    nAcq = data_in.ray_blocks/nKeys;%4;%11;
end
samplesPerAcq = nPts*nRaysPerKey*nKeys;
mat_size = 2*nPts*[1 1 1];
overgrid_mat_size = ceil(mat_size.*overgridding);
overgridding = overgrid_mat_size./mat_size;

%% Sliding window parameters
keysPerTimePoint = nRaysPerKey*nKeys*nPts; % 
timePointStep_sampleUnits = round(keysPerTimePoint/nKeys); % percent of keysPerTimePoint

%% Read in fid data, put in
% dataFile = [ u_dir  '/Desktop/B03180/fid'];
% fid = fopen(dataFile);
% data = fread(fid,inf,'int32');
% fclose(fid);
% data = complex(data(1:2:end),data(2:2:end)); % Make the data complex
% from interleaved complex
expected_dims=[nPts nCoils nRaysPerKey nKeys nAcq];
while expected_dims(end)==1
    expected_dims(end)=[];
end
if ~isprop(data_buffer,'radial')
    if ~isprop(data_buffer,'data')
        fid = fopen(opt_struct.dataFile);
        data = fread(fid,inf,'int32');
        fclose(fid);
        data = complex(data(1:2:end),data(2:2:end)); % Make the data actually complex
        data = reshape(data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
        data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
        data = reshape(data,[nPts nRaysPerKey*nKeys nAcq nCoils]); % Vectorize all but coil dimmension
    else
        data = reshape(data_buffer.data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
        if numel(data) ~= prod(expected_dims) ...
                || sum(size(data) ~= expected_dims)>0
            db_inplace('scott_grid','data didnt shape up');
        end
        data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
        data = reshape(data,[nPts nRaysPerKey*nKeys nAcq nCoils]); % Vectorize all but coil dimmension
    end
    data_buffer.addprop('radial');
    data_buffer.radial=data;
else
    data=data_buffer.radial;
end

%% Read in trajectory
% trajFile = [ u_dir  '/Desktop/B03180/traj'];
% fid = fopen(trajFile);
% traj = fread(fid,inf,'double');
% fclose(fid);
% traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
expected_dims=[3 nPts nRaysPerKey nKeys];
if ~isprop(data_buffer,'straj')
    if ~isprop(data_buffer,'trajectory')
        trajectory_name='traj';
        base_path=fileparts(opt_struct.dataFile);
        trajFile=[base_path '/' trajectory_name ];
        fid = fopen(trajFile);
        traj = fread(fid,inf,'double');
        fclose(fid);
        traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
        traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
        traj = reshape(traj,[nPts nRaysPerKey*nKeys 3]); % vectorize keys and Acq
    else
        if numel(data_buffer.trajectory) ~= prod(expected_dims) %...
               % || sum(size(data_buffer.trajectory) ~= expected_dims)>0
            db_inplace('sliding_window_recon','traj didnt shape up');
        end
        traj = permute(data_buffer.trajectory,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
    end
    traj = reshape(traj,[nPts nRaysPerKey*nKeys 3]); % vectorize keys and Acq
    data_buffer.addprop('straj');
    data_buffer.straj=traj;
else
    traj=data_buffer.straj;
end

% % Vectorize data and traj
% traj is rl*rpk*k*3
traj = reshape(traj,[samplesPerAcq 3])'; % vectorize all but [kx,ky,kz] dimmension
data = reshape(data,[samplesPerAcq*nAcq nCoils])'; % Vectorize all but coil dimmension

%% Perform sliding window Reconstruction
totalSamples = samplesPerAcq*nAcq; % all keys, acqs
% windowStartIdxs = 1:timePointStep_sampleUnits:(totalSamples-keysPerTimePoint-1);
windowStartIdxs = 1:timePointStep_sampleUnits:(totalSamples-keysPerTimePoint+1);
nWindows = length(windowStartIdxs);

% tmpVol = zeros(overgrid_mat_size);
% tmpData = zeros([1 keysPerTimePoint]);
% tmpComplexVol = complex(zeros(overgrid_mat_size));
% tmpComplexData = complex(zeros([1 keysPerTimePoint]));
% sosComplexVol = complex(zeros(overgrid_mat_size));
% dcf = complex(zeros([1 keysPerTimePoint]));
% if(cropVolume)
%     cropSosVol = complex(zeros(mat_size));
% end

if(saveFullVol)
    warning('About to make a really big matrix!');
    if(cropVolume)
        fullVolume = complex(zeros([mat_size nWindows]));
    else
        fullVolume = complex(zeros([overgrid_mat_size nWindows]));
    end
end

start_time=tic;
% startingTime = toc;
% parfor iWin = 1:nWindows
for iWin = 1:nWindows
%     disp(['Reconstructing ' num2str(iWin) '/' num2str(nWindows) ' traj subsets']);
    tmpVol = zeros(overgrid_mat_size);
    tmpData = zeros([1 keysPerTimePoint]);
    tmpComplexVol = complex(zeros(overgrid_mat_size));
    tmpComplexData = complex(zeros([1 keysPerTimePoint]));
    sosComplexVol = complex(zeros(overgrid_mat_size));
    dcf = complex(zeros([1 keysPerTimePoint]));
    if(cropVolume)
        cropSosVol = complex(zeros(mat_size));
    end
    
    % Make a trajectory for each window
    windowKeys = windowStartIdxs(iWin)+[1:keysPerTimePoint]-1;
    if(nAcq > 1)
        windowKeysTraj = mod(windowStartIdxs(iWin)+[1:keysPerTimePoint]-1,samplesPerAcq);
        windowKeysTraj(windowKeysTraj==0)=samplesPerAcq;
        windowTraj = squeeze(traj(:,windowKeysTraj));
    else
        windowTraj = squeeze(traj(:,windowKeys));
    end
    
    % Construct system model for this trajectory

    A = SparseGridder(windowTraj,mat_size,sharpness,extent,overgridding,opt_struct.nThreads);
       
       
    % Option 1 - dcf
%     disp('   Computing density compensation weights');
    tmpVol = ones(A.overgridSize);
    tmpData = A.ungrid(tmpVol);
    dcf = 1./tmpData;
    for iDcfIter=1:nPipeIter
%         disp(['      DCF iter ' num2str(iDcfIter) '/' num2str(nPipeIter)])

        A.grid(dcf,tmpVol);
        A.ungrid(tmpVol,tmpData);
        fig_id=disp_vol_center(tmpVol,1,100+iDcfIter);
        if fig_id>0
            set(fig_id,'Name',sprintf('kspace_dcf_i%i',iDcfIter));
        end
        dcf = dcf./tmpData;
    end
    
    % Create a data matrix of all repetitions of this trajectory
    windowData = double(data(:,windowKeys));
    
    % Grid data from all coils and compute SOS
    sosComplexVol = complex(zeros(overgrid_mat_size)); % reset to zero
    for iCoil = 1:nCoils
        % Apply dcf
        windowData(iCoil,:)  = windowData(iCoil,:).*dcf;
        
        % Calculate gridded kspace
        A.grid(windowData(iCoil,:),tmpComplexVol);
        
        % Reconstruct image domain with IFFT
        tmpComplexVol = ifftshift(ifftn(tmpComplexVol));
        fig_id=disp_vol_center(tmpComplexVol,0,200+iCoil);
        if fig_id>0
            set(fig_id,'Name',sprintf('tmpComplex_c%i',iCoil));
        end
        % Accumulate SOS
        sosComplexVol = sosComplexVol + tmpComplexVol.^2;
%         disp(['      Finished Coil ' num2str(iCoil) '/' num2str(nCoils)]);
    end
    
    % Finish SOS
    sosComplexVol = sqrt(sosComplexVol);
    fig_id=disp_vol_center(sosComplexVol,0,310);
    if fig_id>0
        set(fig_id,'Name',sprintf('sosComplex'));
    end
    % Compute deapodization volume for this traj
    if(deapodize)
%         disp('   Deapodizing...');
        tmpData = ~any(windowTraj,1).*dcf;
        A.grid(tmpData,tmpVol);
        tmpVol = ifftshift(ifftn(tmpVol));
        fig_id=disp_vol_center(tmpVol,0,311);
        if fig_id>0
            set(fig_id,'Name',sprintf('deapodize_filter'));
        end
        if sum(nonzeros(tmpVol(:)))==0
            warning('Deapoidze failure, cannot deapodize');
        else
            sosComplexVol = sosComplexVol./tmpVol;
        end
    end
    fig_id=disp_vol_center(sosComplexVol,0,312);
    if fig_id>0
        set(fig_id,'Name',sprintf('Complete_image_post_deapodize'));
    end
    % Crop volume
    if(cropVolume)
%         disp('   Cropping volume...');
        cropSosVol = subvolume(sosComplexVol,...
            [round([0.5*(A.overgridSize-A.gridSize)+1]); ...
            round([0.5*(A.overgridSize+A.gridSize)])]);
    end
    
%     % Save/store this time point
%     if(~saveFullVol)
%         disp('   Saving Data');
%         if(cropVolume)
%             nii = make_nii(abs(cropSosVol));
%         else
%             nii = make_nii(abs(sosComplexVol));
%         end
%         save_nii(nii,['reconTime_' num2str(iWin) '.nii']);
%     else
        if(cropVolume)
            fullVolume(:,:,:,iWin) = cropSosVol;
        else
            fullVolume(:,:,:,iWin) = sosComplexVol;
        end
%     end
    
    disp('Completed another window')
%     % Compute remaining time
%     currentTime = toc;
%     timeSoFar = currentTime - startingTime;
%     timePerIter = timeSoFar/iWin;
%     totalEstTime = timePerIter*nWindows;
%     remainingTime = totalEstTime - timeSoFar;
%     
%     % Show some progress
%     disp(['Completed ' num2str(iWin) '/' num2str(nWindows) ' traj subsets, est total time =~' num2str(totalEstTime) ' sec (' num2str(timeSoFar) ' sec so far, ~' num2str(remainingTime) ' sec remaining)']);
end
reconTime = toc(start_time);

fprintf('reconstruct in %f minutes\n',reconTime/60);
if ~isprop(data_buffer,'data')
    data_buffer.addprop('data');
end
if exist('fullVolume','var')
    data_buffer.data=fullVolume;
end
% figure(2);plot3(ucf.t1(:,:,1),ucf.t1(:,:,2),ucf.t1(:,:,3),'.');
% Show the reconstruction
% imslice(slidingWindowReconVol,'Sliding window');
