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
% retrievalable via "git clone https://git@github.com/jamesjcook/civm_matlab_comon_utils"
 
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
        sliding_window_recon(data_buffer,opt_struct,data_in,data_work,data_out);
        return;
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
saveFullVol = 1;

%% Read the header file to get scan info
% nCoils      = coil count in acquision
% nPts        = ray_lengh+ramp_points
% nRaysPerKey = rays per step of trajectory
% nKeys       = number of steps in full trajectory
% nAcq        = repetitions of trajectory
nKeys = 13;% in ergys terminology a fully sampled acq has this many steps.

%The number of repetitions is called 'ray_block':
%data_buffer.headfile.ray_blocks

if isfield(data_buffer.headfile,'B_vol_type') ... 
        && isempty(regexpi(data_buffer.headfile.B_vol_type,'.*keyhole.*'))
    nKeys=1;
end
% To keep this code as it was, for new acquisitions this must be the number
% of steps of the trajectory.
% In john's super acquisition 16 would work, but that produces many volumes(57*16).
% Recommended 2, 
%{
nKeys=2; % using  factor(51360) found the factors ( 2 2 2 2 2 3 5 107 )
nRaysPerKey=107*3*5*2*2*2*2;
nAcq=57;
%}

if ~exist ('data_in','var')
    nPts=64;
    nCoils=4;
    nRaysPerKey=1980;
    nAcq=11;
else
    % This somewhat nicely reads the header. 
    % For expanding the sliding window code this needs to be expanded/overridden.
    nPts = data_in.ray_length;%64;
    if isfield(data_in,'ramp_points')
        nPts=nPts+data_in.ramp_points;%64+ramp points;
    end
    nCoils = data_in.ds.Sub('c');%2;%4;          % coil count in acquision
    nRaysPerKey = data_in.rays_per_block;%1980;  % stepsize of trajectory
    nAcq = data_in.ray_blocks/nKeys;%4;%11;      % total number of output volumes AND steps in the data.
    
    %===================== ADDED BY JOHN =====================%
    
%     %Do this part without KFiltering
%     data_in.nRepetitions=nAcq; %number of times the same group of trajectories was repeated in the acquisition
%     data_in.nRaysPerRepetition=data_in.rays_per_block; %205634 for 256^3
%     data_in.FirstReconIndex=25704; %Quarter a Full Sample
%     data_in.ReconstructionIncrement=8000 %1/16th of a Full Sample        6425 %1/16th of a Full Sample  %25704;         % number of rays (+1) between subsequent reconstructions
%     data_in.WindowWidth=8000; %1/16th of a Full Sample      %Sliding window width (maximal width), MUST BE ODD NUMBER
%     data_in.nRaysFullWeightAtKZero=0; %not used when no k-filtration
%     data_in.nRaysFullWeightAtKMax=data_in.WindowWidth;
%     data_in.DoKFiltering=0;
    
    %Do this part without KFiltering
    data_in.nRepetitions=nAcq; %number of times the same group of trajectories was repeated in the acquisition
    data_in.nRaysPerRepetition=data_in.rays_per_block; %205634 for 256^3
    data_in.FirstReconIndex=25704; %Quarter a Full Sample
    data_in.ReconstructionIncrement=8000 %1/16th of a Full Sample        6425 %1/16th of a Full Sample  %25704;         % number of rays (+1) between subsequent reconstructions
    data_in.WindowWidth=51408-1; %1/16th of a Full Sample      %Sliding window width (maximal width), MUST BE ODD NUMBER
    data_in.nRaysFullWeightAtKZero=0; %not used when no k-filtration
    data_in.nRaysFullWeightAtKMax=data_in.WindowWidth;
    data_in.DoKFiltering=0;
    
%     %Do part with KFiltering
%     data_in.nRepetitions=nAcq; %number of times the same group of trajectories was repeated in the acquisition
%     data_in.nRaysPerRepetition=data_in.rays_per_block; %205634 for 256^3
%     data_in.FirstReconIndex=25704;
%     data_in.ReconstructionIncrement=8000 %12847 %1/16th of a Full Sample       %25704;         % number of rays (+1) between subsequent reconstructions
%     data_in.WindowWidth=51408-1; %1/4th of a Full Sample      %Sliding window width (maximal width), MUST BE ODD NUMBER
%     data_in.nRaysFullWeightAtKZero=8000; %1/16th of a Full Sample  %Must be odd
%     data_in.nRaysFullWeightAtKMax=data_in.WindowWidth;
%     data_in.DoKFiltering=1;

    if mod(data_in.WindowWidth,2)==0
        data_in.WindowWidth=data_in.WindowWidth-1;
    end;
    if mod(data_in.nRaysFullWeightAtKZero,2)==0
        data_in.nRaysFullWeightAtKZero=data_in.nRaysFullWeightAtKZero-1;
    end;
    
%     %Test
%     data_in.nRaysAcquired=18;
%     data_in.WindowWidth=9;
%     data_in.FirstUsed_RayIndex=3;
%     data_in.ReconstructionIncrement=3;
    
    data_in.nRaysAcquired=data_in.nRepetitions*data_in.nRaysPerRepetition;
    data_in.FirstUsed_RayIndex=data_in.FirstReconIndex-(data_in.WindowWidth-1)/2; %We could discard and leave unused a number of rays at the beginning of the data
    nRaysSurroundedByFullWindow=floor((data_in.nRaysAcquired-(data_in.WindowWidth-1)-data_in.FirstUsed_RayIndex)/...
        data_in.ReconstructionIncrement)*data_in.ReconstructionIncrement+1;
    data_in.nReconstructions=1+floor((nRaysSurroundedByFullWindow-1)/data_in.ReconstructionIncrement);  % number fo volumes reconstructed (or time points)
    data_in.FirstRecon_RayIndex=(data_in.WindowWidth-1)/2+data_in.FirstUsed_RayIndex; 
    data_in.LastRecon_RayIndex=(data_in.WindowWidth-1)/2+data_in.FirstUsed_RayIndex+(data_in.nReconstructions-1)*data_in.ReconstructionIncrement; 
    data_in.LastUsed_RayIndex=data_in.LastRecon_RayIndex+(data_in.WindowWidth-1)/2; 

end

%samplesPerAcq = nPts*data_in.ReconstructionIncrement*nKeys;
mat_size = 2*nPts*[1 1 1];
if isfield(data_in,'ramp_points')
    mat_size = 2*(nPts-data_in.ramp_points)*[1 1 1];
end

overgrid_mat_size = ceil(mat_size.*overgridding);
overgridding = overgrid_mat_size./mat_size;

% %% Sliding window parameters
% %     keysPerTimePoint = nRaysPerKey*nKeys*nPts; % 
% %     timePointStep_sampleUnits = round(keysPerTimePoint/nKeys); % percent of keysPerTimePoint
% 
%     %===================== ADDED BY JOHN =====================%
%     nPtsPerReconstruction=nPts*data_in.WindowWidth; % number of points included in a reconstruction increment
%     data_in.ReconstructionIncrement=10000;       % number of rays (+1) between subsequent reconstructions
%     data_in.nReconstructions=floor(data_in.ray_blocks*data_in.rays_per_block/data_in.ReconstructionIncrement)  % number fo volumes reconstructed (or time points)
%     
%     % transfer to legacy terminology
%     keysPerTimePoint = nPtsPerReconstructionIncrement;
%     timePointStep_sampleUnits = round(keysPerTimePoint/nKeys); % percent of keysPerTimePoint
%     
%     %=================== END ADDED BY JOHN ===================%
% 
%% Read/reshape fid data, 
expected_dims=[nPts nCoils data_in.WindowWidth data_in.nReconstructions];
%expected_dims=[nPts nCoils data_in.nRaysPerRepetition data_in.nRepetitions];
%expected_dims=[nPts nCoils nRaysPerKey nKeys nAcq];
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
        
        % KEYS MUST BE REMOVED FROM ABOVE
        
        
    else
        %===================== ADDED BY JOHN =====================%
        %the original size is Y x data_in.nRepetitions x data_in.nRaysPerRepetition
        %we need: Y x x data_in.nReconstructions x data_in.ReconstructionIncrement
        Y=length(data_buffer.data)/(data_in.nRepetitions*data_in.nRaysPerRepetition);
        
%         data_in.ReconstructionIncrement=3;
%         data_in.WindowWidth=9;
%         data_in.FirstUsed_RayIndex=2;
        data=zeros(1,(data_in.WindowWidth*data_in.nReconstructions)*Y);
        
        %We Create the K-Filter       
        KFilter=logical(ones(nPts,data_in.WindowWidth,1,1));
        for ThisKPoint=1:nPts
            if ThisKPoint<=data_in.ramp_points
                NumberOfRaysAtFullWeight=data_in.nRaysFullWeightAtKZero;
            else
                ThisKPointAfterRamp=ThisKPoint-data_in.ramp_points;
                NumberOfRaysAtFullWeight=...
                    (data_in.WindowWidth-data_in.nRaysFullWeightAtKZero)/data_in.ray_length...
                    *ThisKPointAfterRamp+data_in.nRaysFullWeightAtKZero;
            end;
            
            for ThisRayIndex=1:data_in.WindowWidth
                if or(ThisRayIndex<(data_in.WindowWidth-NumberOfRaysAtFullWeight)/2,ThisRayIndex>(data_in.WindowWidth-NumberOfRaysAtFullWeight)/2+NumberOfRaysAtFullWeight)
                    ToggleValue=0;
                else
                    ToggleValue=1;
                end;
                if ~data_in.DoKFiltering %FORCES 'NO FILTER'
                    ToggleValue=1;
                end;
                KFilter(ThisKPoint,ThisRayIndex,1,:)=ToggleValue;
            end;
        end;
        
        figure(3),subplot(2,2,1),imagesc(KFilter),colorbar,drawnow;            
        DensityCompensation=size(KFilter,2);
        KFilterCompensated=KFilter*DensityCompensation./repmat(sum(KFilter,2),1,size(KFilter,2));
        figure(3),subplot(2,2,2),imagesc(KFilterCompensated),colorbar,drawnow;            
        
        %Sorry for the loop
        for WhichReconIndex=1:data_in.nReconstructions
            ThisRecon_RayIndex=(data_in.WindowWidth-1)/2+data_in.FirstUsed_RayIndex+(WhichReconIndex-1)*data_in.ReconstructionIncrement;
            DataStartIndex=ThisRecon_RayIndex-(data_in.WindowWidth-1)/2;
            DataEndIndex=ThisRecon_RayIndex+(data_in.WindowWidth-1)/2;
            ThisData=data_buffer.data(1,Y*(DataStartIndex-1)+1:Y*DataEndIndex); %Each chunk of data = one Sliding Window Width
            data(1,Y*data_in.WindowWidth*(WhichReconIndex-1)+1:Y*data_in.WindowWidth*WhichReconIndex)=ThisData;
        end;
        clear ThisData;  
        %=================== END ADDED BY JOHN ===================%
        
        %data = reshape(data_buffer.data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
        %data = reshape(data,[nPts nCoils data_in.nRaysPerRepetition data_in.nRepetitions]); % put data into matrix form
        data = reshape(data,[nPts nCoils data_in.WindowWidth data_in.nReconstructions]); % put data into matrix form
        data_buffer.data=[];
        if numel(data) ~= prod(expected_dims) ...
                || sum(size(data) ~= expected_dims)>0
            db_inplace('scott_grid','data didnt shape up');
        end
        data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
        data = reshape(data,[nPts data_in.WindowWidth data_in.nReconstructions nCoils]); % Vectorize all but coil dimmension
        figure(3),subplot(2,2,3),imagesc(squeeze(abs(data(:,:,1,1)))),colorbar,drawnow;
        data = data.*repmat(KFilterCompensated,[1 1 data_in.nReconstructions nCoils]);
        figure(3),subplot(2,2,4),imagesc(squeeze(abs(data(:,:,1,1)))),colorbar,drawnow;            

    end
    data_buffer.addprop('radial');
    data_buffer.radial=data;
else
    data=data_buffer.radial;
end

%% Read/reshape in trajectory
%===================== ADDED BY JOHN =====================%
expected_dims=[3 nPts data_in.WindowWidth data_in.nReconstructions];
%expected_dims=[3 nPts nRaysPerKey nKeys];
traj_rep=repmat(data_buffer.trajectory,[1 1 data_in.nRepetitions]);
data_buffer.trajectory=[];

%Sorry for the loop
data_traj=zeros(3,nPts,data_in.WindowWidth*data_in.nReconstructions);
for WhichReconIndex=1:data_in.nReconstructions
    ThisRecon_RayIndex=(data_in.WindowWidth-1)/2+data_in.FirstUsed_RayIndex+(WhichReconIndex-1)*data_in.ReconstructionIncrement;
    DataStartIndex=ThisRecon_RayIndex-(data_in.WindowWidth-1)/2;
    DataEndIndex=ThisRecon_RayIndex+(data_in.WindowWidth-1)/2;
    ThisData=traj_rep(:,:,DataStartIndex:DataEndIndex);
    data_traj(:,:,(WhichReconIndex-1)*data_in.WindowWidth+1:WhichReconIndex*data_in.WindowWidth)=ThisData;
end;
clear traj_rep;

if ~isprop(data_buffer,'straj')
    if ~isprop(data_buffer,'trajectory')
        trajectory_name='traj';
        base_path=fileparts(opt_struct.dataFile);
        trajFile=[base_path '/' trajectory_name ];
        fid = fopen(trajFile);
        traj = fread(fid,inf,'double');
        fclose(fid);
        
        %===================== ADDED BY JOHN =====================%
        %we need the file to be longer (include all repetitions)
        %THIS PART OF CODE NOT FUNCTIONAL YET
        traj = repmat(traj,[1 data_in.nRepetitions]);
        traj = reshape(traj,[nPts nRaysPerKey*nKeys*data_in.nRepetitions 3]); % vectorize keys and Acq
        traj = reshape(traj,[3 nPts nRaysPerKey*data_in.nRepetitions nKeys]); % put trajectory into matrix form
        traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
        %=================== END ADDED BY JOHN ===================%
        
%         traj = reshape(traj,[nPts nRaysPerKey*nKeys 3]); % vectorize keys and Acq
%         traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
%         traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
        
        
    else
        if numel(data_traj) ~= prod(expected_dims) %...
        %if numel(data_buffer.trajectory) ~= prod(expected_dims) %...
               % || sum(size(data_buffer.trajectory) ~= expected_dims)>0
            db_inplace('sliding_window_recon','traj didnt shape up');
        end
        traj = permute(data_traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
        data_traj=[];
    end
    traj = reshape(traj,[nPts data_in.WindowWidth*data_in.nReconstructions 3]);
%    traj = reshape(traj,[nPts nRaysPerKey*nKeys 3]); % vectorize keys and Acq
    data_buffer.addprop('straj');
    data_buffer.straj=traj;
else
    traj=data_buffer.straj;
end

%% Vectorize data and traj
% traj is rl*rpk*k*3
traj = reshape(traj,[nPts*data_in.WindowWidth*data_in.nReconstructions 3])'; % vectorize all but [kx,ky,kz] dimmension
data = reshape(data,[nPts*data_in.WindowWidth*data_in.nReconstructions nCoils])'; % Vectorize all but coil dimmension
%traj = reshape(traj,[samplesPerAcq 3])'; % vectorize all but [kx,ky,kz] dimmension
%data = reshape(data,[samplesPerAcq*nAcq nCoils])'; % Vectorize all but coil dimmension

%% Perform sliding window Reconstruction
totalSamples = nPts*data_in.WindowWidth*data_in.nReconstructions; % all keys, acqs
%totalSamples = samplesPerAcq*nAcq; % all keys, acqs
% windowStartIdxs = 1:timePointStep_sampleUnits:(totalSamples-keysPerTimePoint+1);
% nWindows = length(windowStartIdxs);
windowStartIdxs=nPts*(data_in.FirstUsed_RayIndex-1)+1:nPts*data_in.WindowWidth:nPts*(data_in.LastRecon_RayIndex-(data_in.WindowWidth-1)/2);
%windowStartIdxs=data_in.FirstUsed_RayIndex:data_in.WindowWidth:data_in.LastRecon_RayIndex-(data_in.WindowWidth-1)/2;
nWindows=data_in.nReconstructions;
keysPerTimePoint=data_in.WindowWidth*nPts;
samplesPerAcq=data_in.WindowWidth*nPts;

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
tmpComplexVol = complex(tmpVol,tmpVol);
sosComplexVol = complex(tmpVol,tmpVol);
tmpData = zeros([1 keysPerTimePoint]);
tmpComplexData = complex(tmpData,tmpData);
dcf = complex(tmpData,tmpData);

% startingTime = toc;
% parfor iWin = 1:nWindows
fprintf('begin grid operations by window\n');
for iWin=1:nWindows
%     disp(['Reconstructing ' num2str(iWin) '/' num2str(nWindows) ' traj subsets']);
    if(cropVolume)
        tv=zeros(mat_size);
        cropSosVol = complex(tv,tv);clear tv;
    end
    
    %% Make a trajectory for this window
%     windowKeys = windowStartIdxs(iWin)+[1:keysPerTimePoint]-1;
    if(data_in.nReconstructions>1)
%         windowKeysTraj = mod(windowStartIdxs(iWin)+[1:keysPerTimePoint]-1,samplesPerAcq);
%         windowKeysTraj(windowKeysTraj==0)=samplesPerAcq;
        windowTraj=traj(:,(iWin-1)*data_in.WindowWidth*nPts+1:iWin*data_in.WindowWidth*nPts);
    else
        windowTraj = squeeze(traj(:,(iWin-1)*data_in.WindowWidth*nPts+1:iWin*data_in.WindowWidth*nPts));
        %THIS CASE HAS NOT BEEN RESOLVED
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
    windowData = double(data(:,(iWin-1)*data_in.WindowWidth*nPts+1:iWin*data_in.WindowWidth*nPts));
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
    if(saveFullVol)
        disp('   Saving Data');
        if(cropVolume)
            nii = make_nii(abs(cropSosVol));
        else
            nii = make_nii(abs(sosComplexVol));
        end
        save_nii(nii,[save_location '/reconWin_' num2str(iWin) '.nii']);
        fprintf('Saving Window %i completed.\n',iWin);

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
