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
% Relies on civm_matlab_common_utils 
% retrievalable via git clone
% https://git@github.com/jamesjcook/civm_matlab_comon_utils
 
save_location='/Volumes/CivmUsers/omega/';
if ~exist(save_location,'dir')
    error('Fix save locaiton ');
end

%
% code_path='/Users/james/Desktop/DCE_proto';
% run([code_path '/Duke-CIVM-MRI-Tools/setup.m']);
% run([code_path '/GE-MRI-Tools/setup.m']);
% run([code_path '/Non-Cartesian-Reconstruction/setup.m'])
% u_dir='/Users/james/';
p=fileparts(mfilename('fullpath'));
if ~exist('data_buffer','var')
    if ~exist(fullfile(p,'sliding_window_bigtest.mat'),'file')&& ~exist('sliding_window_bigtest.mat','file')
        error('MUST HAVE A DATABUFFER TO CATCH OUTPUT');
    else
        warning('Loading test data!'); pause(2);
        if exist(fullfile(p,'sliding_window_bigtest.mat'),'file')
            load(fullfile(p,'sliding_window_bigtest.mat'));
        else
            load('sliding_window_bigtest.mat');
        end
        sliding_window_recon(data_buffer,opt_struct,data_in,data_work,data_out);return;
    end
end
if ~exist(fullfile(p,'sliding_window_bigtest.mat'),'file')
    save(fullfile(p,'sliding_window_bigtest.mat'),'data_buffer','opt_struct','data_in','data_work','data_out','-v7.3');
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

%% VARIABLE OVERRRIDES
saveFullVol = 0;

    
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
% nPts=ray_lengh
% nRaysPerKey = as stated
% nKeys is number of divisions of the trajectory
samplesPerAcq = nPts*nRaysPerKey*nKeys;
mat_size = 2*nPts*[1 1 1];
if isfield(data_in,'ramp_points')
    %nPts=nPts
    mat_size = 2*(nPts-data_in.ramp_points)*[1 1 1];
end

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
        data_buffer.data=[];
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
        data_buffer.trajectory=[];
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


if(saveFullVol)
    warning('About to make a really big matrix!');
    if(cropVolume)
        tmpVol=zeros([mat_size nWindows],'single');
    else
        tmpVol=zeros([overgrid_mat_size nWindows],'single');
    end
    fullVolume = complex(tmpVol,tmpVol);
end

start_time=tic;
%% preallocate all arrays in memory, so we can have constant mem useage.
tmpVol = zeros(overgrid_mat_size);% consider magBlank as we'll use this later
%tmpComplexVol = complex(zeros(overgrid_mat_size));% oddly slow
tmpComplexVol = complex(tmpVol,tmpVol);
%sosComplexVol = complex(zeros(overgrid_mat_size)); % oddly slow
sosComplexVol = complex(tmpVol,tmpVol);
tmpData = zeros([1 keysPerTimePoint]);
%tmpComplexData = complex(zeros([1 keysPerTimePoint]));% oddly slow
tmpComplexData = complex(tmpData,tmpData);
%dcf = complex(zeros([1 keysPerTimePoint]));
dcf = complex(tmpData,tmpData);

% startingTime = toc;
% parfor iWin = 1:nWindows
fprintf('begin grid operations by window\n');
for iWin = 1:nWindows
%     disp(['Reconstructing ' num2str(iWin) '/' num2str(nWindows) ' traj subsets']);
    if(cropVolume)
        tv=zeros(mat_size);
        cropSosVol = complex(tv,tv);clear tv;
    end
    
    %% Make a trajectory for this window
    windowKeys = windowStartIdxs(iWin)+[1:keysPerTimePoint]-1;
    if(nAcq > 1)
        windowKeysTraj = mod(windowStartIdxs(iWin)+[1:keysPerTimePoint]-1,samplesPerAcq);
        windowKeysTraj(windowKeysTraj==0)=samplesPerAcq;
        windowTraj = squeeze(traj(:,windowKeysTraj));
    else
        windowTraj = squeeze(traj(:,windowKeys));
    end
    
    %% Construct system model for this trajectory This is a BIG OPERATION.
    % It may have room for optimization under the hood.
    % Also trashes previous copies, so no point in keeping an old one in
    % memory. 
    clear SG;
    fprintf('Creating SparseGrid window %i, massive memory expansion now.\n',iWin);
    SG = SparseGridder(windowTraj,mat_size,sharpness,extent,overgridding,opt_struct.nThreads);% mem rise 34->77->75
       
       
    % Option 1 - dcf
    %% Computing density compensation weights
    fprintf('Computing density compenation\n');
    clear tmpVol;
    tmpVol = ones(SG.overgridSize);% mem rise 73->76->71
    tmpData = SG.ungrid(tmpVol);
    dcf = 1./tmpData;% mem rise 71->72
    for iDcfIter=1:nPipeIter
%         disp(['      DCF iter ' num2str(iDcfIter) '/' num2str(nPipeIter)])
        SG.grid(dcf,tmpVol);
        SG.ungrid(tmpVol,tmpData);
        %% display if we can
        if exist('disp_vol_center','file')
            fig_id=disp_vol_center(tmpVol,1,100+iDcfIter);
            if fig_id>0
                set(fig_id,'Name',sprintf('kspace_dcf_i%i',iDcfIter));
            end
        end
        dcf = dcf./tmpData;
    end
    
    %% Grid data from all coils and compute SOS,
    % Must set data to 0 as we accumulate it into sos. This can cause a
    % memory surge depending on how its done.
    % Both of these cause a surge.
    %sosComplexVol = complex(zeros(overgrid_mat_size)); % reset to zero
    %sosComplexVol = complex(zeros(overgrid_mat_size),zeros(overgrid_mat_size)); % reset to zero
    clear tmpVol tmpComplexVol sosComplexVol;
    tmpVol=zeros(SG.overgridSize,'single');sosComplexVol=complex(tmpVol,tmpVol);tmpVol=double(tmpVol);tmpComplexVol=complex(tmpVol,tmpVol);
    %%% these two may be more memory efficient.
    % clear sosComplexVol;sosComplexVol=zeros(overgrid_mat_size);sosComplexVol=complex(sosComplexVol,sosComplexVol);
    % sosComplexVol = sosComplexVol.*(0+0i); % mem rise 74->85, interestingly this just delays when memory is used, and slows down later operations.
    
    % Create a data matrix of all repetitions of this trajectory
    windowData = double(data(:,windowKeys));
    for iCoil = 1:nCoils
        % Apply dcf
        windowData(iCoil,:)  = windowData(iCoil,:).*dcf;
        % Calculate gridded kspace
        SG.grid(windowData(iCoil,:),tmpComplexVol);
        % Reconstruct image domain with IFFT
        % tmpComplexVol = ifftn(tmpComplexVol);
        tmpComplexVol = (ifft(ifft(ifft(tmpComplexVol,[],1),[],2),[],3)); % In testing this way is ever so slightly faster.
        % Accumulate SOS
        sosComplexVol = sosComplexVol + tmpComplexVol.^2;% mem rise 79->90->85
    end
    % Finish SOS
    sosComplexVol = sqrt(sosComplexVol);
    sosComplexVol = ifftshift(sosComplexVol);
    if exist('disp_vol_center','file')
        fig_id=disp_vol_center(sosComplexVol,0,310);
        if fig_id>0
            set(fig_id,'Name',sprintf('sosComplex'));
        end
    end
    % clear tmpComplexVol windowData;
    %% Compute deapodization volume for this traj
    if(deapodize)
%         disp('   Deapodizing...');
        tmpData = ~any(windowTraj,1).*dcf;
        SG.grid(tmpData,tmpVol);
        tmpVol = ifftshift(ifftn(tmpVol));% mem rise 85->96->85
        if exist('disp_vol_center','file')
            fig_id=disp_vol_center(tmpVol,0,311);
            if fig_id>0
                set(fig_id,'Name',sprintf('deapodize_filter'));
            end
        end
        if sum(nonzeros(tmpVol(:)))==0
            warning('Deapoidze failure, cannot deapodize');
        else
            sosComplexVol = sosComplexVol./tmpVol;
        end
    end
    % clear tmpData;
    if exist('disp_vol_center','file')
        fig_id=disp_vol_center(sosComplexVol,0,312);
        if fig_id>0
            set(fig_id,'Name',sprintf('Complete_image_post_deapodize'));
        end
    end
    % Crop volume
    if(cropVolume)
%         disp('   Cropping volume...');
        cropSosVol = subvolume(sosComplexVol,...
            [round([0.5*(SG.overgridSize-SG.gridSize)+1]); ...
            round([0.5*(SG.overgridSize+SG.gridSize)])]);
    end
    
    %% Save/store this time point
    if(~saveFullVol)
        disp('   Saving Data');
        if(cropVolume)
            nii = make_nii(abs(cropSosVol));
        else
            nii = make_nii(abs(sosComplexVol));
        end
        save_nii(nii,[save_location '/reconWin_' num2str(iWin) '.nii']);
    else
        if(cropVolume)
            fullVolume(:,:,:,iWin) = cropSosVol;
        else
            fullVolume(:,:,:,iWin) = sosComplexVol;
        end
    end
    
    fprintf('Window %i completed.\n',iWin);
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
