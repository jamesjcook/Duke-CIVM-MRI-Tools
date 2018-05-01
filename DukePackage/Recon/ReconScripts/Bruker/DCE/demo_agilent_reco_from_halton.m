%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% A demonstration of basic radial reconstruction
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IAMSCOTT=false
%% sort out code paths and setup stuff.
if ~IAMSCOTT
    % code_path='/Users/james/Desktop/DCE_proto';
    WKS_HOME=getenv('WORKSTATION_HOME');
    if isempty(WKS_HOME)
        WKS_cf=fullfile('/Volumes/workstation_home/','software');
    else
        WKS_cf=WKS_HOME;
    end
    if isdir(WKS_cf)
        %code_path='/Volumes/workstation_home/software/recon/External/ScottHaileRobertson';
        code_path=fullfile(WKS_HOME,'recon','External','ScottHaileRobertson');
        packages=strsplit('Duke-CIVM-MRI-Tools GE-MRI-Tools Non-Cartesian-Reconstruction');
        % compiles fail on for civm cluster, beacuse mex compile has to be done out of band with gcc6.3 due to matlab errors.
        for p=1:numel(packages)
            run(fullfile(code_path,packages{p},'setup.m'));
        end
    end
end
try
    profile('-memory','on');
catch
    profile off
    profile('-memory','on');
end
%% Set acq directory containing all information
% or just load a .mat file
if ~IAMSCOTT
    if exist('reconDir','var')
        lastReco=reconDir;
    else
        lastReco='';
    end
    u_dir=fullfile(getenv('BIGGUS_DISKUS'),'S68003_10');
    %% Axial scan ordering.
    %% rad3d002
    % per pt loading implemented in 
    % np=128
    % pt_index=[0,1,2] nv=1024 (3072 total)
%     reconDir = [u_dir '/ser10.fid']; % bs=16
    % pt_index=[0,1,2] nv=4096 (12288 total)
%     reconDir = [u_dir '/ser11.fid']; % 
    % pt_index=[0,1,2,3]  nv=4096 (16384 total)
%     reconDir = [u_dir '/ser12.fid']; %
    % pt_index= pt_index=[0]  nv=20480 (20480 total)
%     reconDir = [u_dir '/ser13.fid']; % sl 163.1256
    % pt_index= pt_index=[0,1,2,3]  nv=16384 (65536 total)
%     reconDir = [u_dir '/ser14.fid']; % sl crashed on cluster, ~ 100-120 GiB used and 1110.655945 seconds for slow
    % np=256 pt_index= pt_index=[0,1,2,3]  nv=1024 (4096 total)
%     reconDir = [u_dir '/ser15.fid']; % sl() fa(68.2400) TR=12
%     reconDir = [u_dir '/ser16.fid']; % sl() fa(67.6260)  TR=24
    % pt_index= pt_index=[0,1,2,3]  nv=16384 (65536 total)
%      reconDir = [u_dir '/ser17.fid']; % sl() og 1.0 fa(137.3230) og 1.2 fa(276.7263)
     % set overgridding to 
%       reconDir = [u_dir '/ser18.fid']; % sl() fa(567.389) 
%       reconDir = [u_dir '/ser19.fid']; % sl() fa(215.8760)
%       reconDir = [u_dir '/ser20.fid']; % sl() fa(1.6719e+03) ) 
     %% rad3ddw014 adding diffusion stuff back, 
     % np=128, pt_index=0, pt_index_skip=1 nv=3072
%      reconDir = [u_dir '/ser24.fid']; % sl() fa(4.9464)  diffusion  panel
%      disabled using minTE and minTR
%      all 1 not sure).
%      reconDir = [u_dir '/ser25.fid']; % sl() fa(4.9592) diffusion off dro=0,dpe=0,dsl=1
%      reconDir = [u_dir '/ser26.fid']; % sl() fa(4.1793) diffusion off dro=0,dpe=1,dsl=0
%      reconDir = [u_dir '/ser27.fid']; % sl() fa(4.1824) diffusion off dro=1,dpe=0,dsl=0
%      reconDir = [u_dir '/ser28.fid']; % sl() fa(4.1868) diffusion off dro=1,dpe=1,dsl=1
%       reconDir = [u_dir '/ser29.fid']; % sl() fa(4.5535) diffusion off dro=0.001,dpe=0.001,dsl=0.001
%       reconDir = [u_dir '/ser30.fid']; % sl() fa(4.2394) diffusion panel disabled, kept TR and TE from the diffusion scans
   
else
end
% should  try to calc this, its based on dwell time. 
% It should be some constant amount of time. 
garbagePts=3;%garbagePts=3;
acqdelay=50e-6;% this is hard set in our radial code right now, in rad3ddw007 forward and rad3d002 forward
trajShift=1;% + numbers add more pts to contibute to k0 
% in testing, what we really want is a trajectory offset, our first pt is
% probably not exactly the 0th point. 

%% gross reco settings fast/slow
reco_speed='slow';
reco_speed='fast';

if strcmp(reco_speed,'fast')
    rs=0; 
    % reconTime =   2.3634 ( og 1.0 )
elseif strcmp(reco_speed,'slow')
    rs=1; 
    % reconTime = 405.7575 ( og 1.5 )
    % reconTime =  92.0684 ( og 1.0 )
else
    
end

%% fine Reconstruction parameters where we set what fast/slow mean
scale = 1;

% For oversampling use at least 2. To tune the value, turn "crop" off 
% on reconObject, then increase this until you dont see signal at the
% edge of FOV, then turn "crop" back on and use that oversampling
% THIS HAS VERY STRONG EFFECT ON MEMORY USAGE.
% at 1.75 with 128 nPts and 12800 rays used over 64GiB's of memory and was killed.
% at 1.5  with 128 nPts and 12800 rays used over 73GiB's memory(horray for mem compression) 
% at 1.5  with 128 nPts and 12800 rays used over 20GiB's memory
% oversampling = 1.75; 
% oversampling = 1.5;
oversampling = 1;
if ~rs
    oversampling = 1;
end
% oversampling = 1.2;
% Sharpness is a key parameter that tradesoff SNR and resolution.
% making sharpness smaller will blurr the object,
% but increase SNR and vice versa.
% sharpness = 0.3; % inital demo code
% sharpness = 0.5; % a supposidly reasonable number
sharpness = 0.7;
if ~rs
sharpness = 0.35; % 0.15 recommended by scott when told we're slow
sharpness = 0.5; % 0.15 recommended by scott when told we're slow
end

%extent = 1*sharpness; % 9 is a good value
extent = 9*sharpness; % 9 is a good value
if ~rs
extent = 6*sharpness; % recommended by scott when told we're slow
end

verbose = 1;
nPipeIter = 3; % Use 3-15 for option 1, lots more (~75) for iterative recon. This should scale recon time linearly(ish)
crop = 1;

%% Read in fid data, put in 
if ~IAMSCOTT
    if isempty(lastReco) ...
            || ~strcmp(lastReco,reconDir)
        system(sprintf('dumpHeader -d0 kamy %s',reconDir));
        ah=read_headfile(sprintf('%s/agilent.headfile',reconDir),1);
        dataFile = [ reconDir '/fid'];
        [RE,IM,nPts,NB,NT,HDR] = load_fid(reconDir);
    end
    %{
    save('agilent_test_S68003_06_ser09.mat','RE','IM','nPts','NB','NT','HDR','ah');
    %}
else
    % load('agilent_test_S68003_06_ser09.mat');
end
% gPts=ceil(acqdelay/(1/ah.z_Agilent_sw));% this appears right. 
% dwellTime=(nPts*2)/ah.z_Agilent_sw
% dwellTime=(2*nPts)/ah.z_Agilent_at
dwellTime=ah.z_Agilent_at/(2*nPts);% dwellTime is per 1/2 complex pt.
gPts=ceil(acqdelay/dwellTime/2);% the /2 *2 stuff is to keep an even multiple of pts.
% at looks to be in ms
% sw is in ?
if garbagePts~=gPts
    msg=sprintf('garbagePts%i proscribed doesnt match auto guess %f',garbagePts,gPts);
    warning(msg);
    %db_inplace(mfilename,msg);
    % garbagePts=round(gPts)
end
data=complex(RE,IM);
nKeys = NB*NT;%nv rad007 its just NB, rad008 is NB*NT
nRaysPerKey = 1;% ns, or maybe NT
%nPts = 128;%np
nCoils = 1;
nAcq = 1;
% data=load_fid(reconDir); this loading code doesnt work on agilent radial
% files.
% the data load was retrieved by using radmat, and debug_stop_load to get
% the data in memory, then saving all vars there.
if 0 % this is bruker load code
fid = fopen(dataFile);
data = fread(fid,inf,'single');
fclose(fid);
data = complex(data(1:2:end),data(2:2:end)); % Make the data actually complex
end
%% reshape data
%%% these vars are here for dividing up our echo train+diffusion to test
%%% getting a better pic using only one element of the echo train.
%{
nCoils=8; 
nAcq=4;
nKeys=nKeys/nAcq/nCoils;
%}
data = reshape(data,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
data = permute(data,[1 3 4 5 2]);% Make coil last dimmension
data = reshape(data,[nPts nRaysPerKey*nKeys nAcq nCoils]); % Vectorize all but coil dimmension
% data plot test
%{
 plot(t,g);hold on
c_i=1;ac_i=1;
for c_i=1:nCoils
figure(1);
d=real(data(:,c_i,1,:,ac_i));
d2=imag(data(:,c_i,1,:,ac_i));
d=squeeze(d);
d2=squeeze(d2);
subplot(2,1,1);
imagesc(log(abs(d)));
subplot(2,1,2);
imagesc(log(abs(d)));

figure(2);
t = 1:numel(d);
subplot(2,1,1);
plot(d)
xlabel('pt');
ylabel('real');
subplot(2,1,2);
plot(d2);
xlabel('pt');
ylabel('imag');


figure(3);
t = 1:numel(d);
subplot(2,1,1);
plot(abs(d))
xlabel('pt');
ylabel('abs real');
subplot(2,1,2);
plot(abs(d2));
xlabel('pt');
ylabel('abs imag');
pause(3);
end

d_max=round(size(data,2)/20);
d_off=0;
for d_off=128:128:3096
d_max=128;
d_idx=(1:d_max)+d_off;
figure();
subplot(3,1,1);
plot(abs(data(:,d_idx))); %plot all data in new fig
subplot(3,1,2);
plot(real(data(:,d_idx))); %plot all data in new fig
subplot(3,1,3);
plot(imag(data(:,d_idx))); %plot all data in new fig
end
figure();imagesc(log(abs(data(:,1:round(size(data,2)/20))))); %implot first 144

%}


%% set scan details from header
samplesPerAcq = nPts*nRaysPerKey*nKeys;
unscaled_size = 2*nPts*[1 1 1];
scaled_output_size = round(scale*unscaled_size);
scale = scaled_output_size(1)/unscaled_size(1);
overgrid_mat_size = 2*ceil(scaled_output_size*oversampling/2);
oversampling=overgrid_mat_size/scaled_output_size; % update over sampling 

if(crop)
    reconMatSize = scaled_output_size;
else
    reconMatSize = overgrid_mat_size;
end

%% Sliding window parameters
% keysPerWindow = nPts*nRaysPerKey*13; 
keysPerWindow = nPts*nRaysPerKey*nKeys; 
% windowStep = nPts*nRaysPerKey*1; % Step by one key of data
windowStep = nPts*nRaysPerKey*nKeys; % Step by one key of data

%% adjust data n points of garbage
if garbagePts>0
    data(1:garbagePts,:)=0;
    %%% initial idea was to move garbage to eprofilend of array where zeros are
    %%% kinda expected. Unclear which gives better results. 
    %data=circshift(data,[-garbagePts,0]);
end

%% display first 5%(up to 128) of data in plot and image. 
d_max=min(round(size(data,2)/20),128);
figure(10+garbagePts);plot(abs(data(:,1:d_max))); %plot all data in new fig
figure(30+garbagePts);imagesc(log(abs(data(:,1:d_max)))); %implot first 5%
%{
figure();imagesc(log(abs(data))); 
%}
%% flesh out the halton sequence
if ~IAMSCOTT
%     haltonFile = [ '/delosspace/'  '/Halton_trajpoints_1GiB.dat'];
    haltonFile=fullfile(getenv('WORKSTATION_DATA'),'trajectory','Halton_trajpoints_1GiB.dat');
    trajFile = [ reconDir  '/traj'];
else
    haltonFile='trajpoints.dat';
    trajFile= 'traj';
end
if ~exist(trajFile,'file')
    sim_traj(trajFile,haltonFile,nPts,nRaysPerKey*nKeys,-(garbagePts+trajShift),0.00*(1/nPts));%
end

%% Read in trajectory
fid = fopen(trajFile);
traj = fread(fid,inf,'double');
fclose(fid);
traj = reshape(traj,[3 nPts nRaysPerKey nKeys]); % put trajectory into matrix form
traj = permute(traj,[2 3 4 1]); % Put [kx,ky,kz] dimmension last
traj = reshape(traj,[nPts nRaysPerKey*nKeys 3])/scale; % vectorize keys and Acq

%{
% if coils are not coils, but rahter foreach point 
traj = reshape(traj,[nPts nCoils nRaysPerKey nKeys nAcq]); % put data into matrix form
traj = permute(traj,[1 3 4 5 2]);% Make coil last dimmension
traj = reshape(traj,[nPts nRaysPerKey*nKeys nAcq nCoils]); % Vectorize all but coil dimmension
%}
if trajShift>0
%     traj=circshift(traj,[-(trajShift+garbagePts),0,0]);
%     traj(end-(trajShift+garbagePts):end,:,:)=0;
end

% % % Throw away any aliased data
% aliased_Pts = any(any(abs(traj)>=0.5,3),2);
% nPts = sum(~aliased_Pts);
% traj = traj(1:nPts,:,:);
% data = data(1:nPts,:,:,:);
% keysPerWindow = nPts*nRaysPerKey*13; % 13 keys of data 
% windowStep = nPts*nRaysPerKey*1; % Step by one key of data
% samplesPerAcq = nPts*nRaysPerKey*nKeys;

% Vectorize data and traj
traj = reshape(traj,[nPts*nRaysPerKey*nKeys 3]); % vectorize all but [kx,ky,kz] dimmension
data = reshape(data,[nPts*nRaysPerKey*nKeys*nAcq nCoils]); % Vectorize all but coil dimmension

%% Create the recon objects that dont change 
fprintf('Recon(%s) ready to start at output_size^3(%i) using %i rays.\n',reco_speed,2*nPts, NB*NT);
startTime=tic;
% Construct Gridding kernel
kernelObj = Recon.SysModel.Kernel.Gaussian(sharpness, extent, verbose);

% Construct Proximity object
proxObj = Recon.SysModel.Proximity.L2Proximity(kernelObj, verbose);

%% Perform sliding window Reconstruction
ptsPerCoil = nPts*nRaysPerKey*nKeys*nAcq; % all keys, acqs
% windowStartIdxs = 1:windowStep:(ptsPerCoil-keysPerWindow);
% just doing one volume. 
windowStartIdxs=1;
nWindows = length(windowStartIdxs);
slidingWindowReconVol = zeros([reconMatSize nWindows]);

% Batch reconstructions of the same trajectories to save on DCF
% calculations and the computation of system matrices
startMod = mod(windowStartIdxs-1,samplesPerAcq)+1;
[uniqueStarts,ia,ic] = unique(startMod); % find all unique trajectories

% and yet he nevr uses either edges or counts.
try % new code 2014b and newer
[uniqueCounts, uniqueEdges] = histcounts(ic,length(uniqueStarts));
catch me % old codeb pre. 
[uniqueCounts, uniqueEdges] = histc(ic,length(uniqueStarts));
end


% Presort traj so parfor works...
shit = struct;
nSysMat = length(uniqueStarts);
% shit.deltaInput = ones([keysPerWindow 1]);
% deapVol = zeros(overgrid_mat_size);
shit.tmpVol = zeros(overgrid_mat_size);
assignVol='tmpVol';
if crop
    shit.cropTmpVol = complex(zeros(reconMatSize));
    assignVol='cropTmpVol';
end
disp(['Completed 0/' num2str(nSysMat) ' traj subsets']);
shit.windowTraj = zeros([keysPerWindow 3]);


for iSysMat = 1:nSysMat % This can be done in parallel, but takes LOTS of memory
    % Make a trajectory for each window
    shit.windowTraj = squeeze(traj(mod(uniqueStarts(iSysMat)+[0:(keysPerWindow-1)]-1,samplesPerAcq)+1,:));
    
    if max(shit.windowTraj(:))>0.5 || min(shit.windowTraj(:))<-0.5  
        warning('naughty trajectory, should be between -0.5 and 0.5');
        shit.windowTraj=shit.windowTraj*0.5 ;
    end
    % Construct system model for this trajectory
    disp('   Creating System model');
    shit.systemObj = Recon.SysModel.MatrixSystemModel(shit.windowTraj, oversampling, ...
        scaled_output_size, proxObj, verbose);   % This can be stored
    
    % Calculate (Pipe) Iterative Density Compensation weights
    disp('   Calculating DCF');
    shit.dcfObj = Recon.DCF.Iterative(shit.systemObj, nPipeIter, verbose); % This can be stored
     
%     % Compute deapodization volume for this traj
%     disp('   Calculating Deapodization');
%     deapVol = shit.systemObj'*(shit.deltaInput.*dcfObj.dcf);
%     deapVol = reshape(full(deapVol),overgrid_mat_size); % make unsparse;
%     deapVol = ifftshift(ifftn(deapVol));
    
    % Create a data matrix of all repetitions of this trajectory
    sameStartIdx = find(ic==iSysMat); 
    nSameStart = length(sameStartIdx);
    shit.dataIdxRep = repmat(windowStartIdxs(sameStartIdx),[keysPerWindow 1]) + repmat([0:(keysPerWindow-1)]',[1 nSameStart]);
    shit.windowData = reshape(data(shit.dataIdxRep(:),:),[keysPerWindow nSameStart*nCoils]);
    shit.dcfRep = repmat(shit.dcfObj.dcf,[1 nSameStart*nCoils]);
%     clear dcfObj;
    
    % Grid all data that share this trajectory
    disp(['   Gridding (' num2str(nSameStart) ' time points)x(' num2str(nCoils) ' coil channels) datasets...']);
    shit.windowData = shit.windowData.*shit.dcfRep;
    shit.ATrans = shit.systemObj.ATrans;
    shit.windowRecon = shit.ATrans*double(shit.windowData); % We will get a huge speek boost if you can get this to take more advantage of CPU
    
    % Perform SOS across coil channels;
    disp(['   Performing IFFT and SOS on ' num2str(nSameStart) ' time points and ' num2str(nCoils) ' coils...']);
    shit.windowRecon = reshape(shit.windowRecon, [size(shit.systemObj.A,2) nSameStart nCoils]);
%     clear systemObj;

    for iSameStart = 1:nSameStart 
        % Figure out this windows index
        iWindow = sameStartIdx(iSameStart);
        
        for iCoil = 1:nCoils
            % Reconstruct image domain with IFFT
            shit.tmpVol = reshape(full(shit.windowRecon(:,iSameStart, iCoil)),overgrid_mat_size); % make unsparse;
            %% display kspace 
            volfig_id=100+garbagePts+10*trajShift+10000*rs;
            if IAMSCOTT
                imslice(abs(slidingWindowReconVol),sprintf('kspace %i',volfig_id));
            else
                disp_vol_center(shit.tmpVol,1,volfig_id);
            end
            %% do ifft
            shit.tmpVol = ifftshift(ifftn(shit.tmpVol));
            if crop
                shit.cropTmpVol = subvolume(shit.tmpVol,...
                    [round([0.5*(size(shit.tmpVol)-size(slidingWindowReconVol))+1]); ...
                    round([0.5*(size(shit.tmpVol)+size(slidingWindowReconVol))])]);
            end
            % Accumulate SOS
            slidingWindowReconVol(:,:,:,iWindow) = slidingWindowReconVol(:,:,:,iWindow) + shit.(assignVol).^2;%(shit.tmpVol.*conj(shit.tmpVol));
            disp(['      Finished Coil ' num2str(iCoil) '/' num2str(nCoils)]);
        end
        % Finish SOS
        slidingWindowReconVol(:,:,:,iWindow) = sqrt(slidingWindowReconVol(:,:,:,iWindow));
        
        % Show some progress
        disp(['   Completed ' num2str(iSameStart) '/' num2str(nSameStart) ' time Points']); 
    end
     
    % Show some progress
    disp(['Completed ' num2str(iSysMat) '/' num2str(nSysMat) ' traj subsets']);
end
reconTime = toc(startTime)
%% Show the reconstruction and/or save it.
volfig_id=200+garbagePts+10*trajShift+10000*rs;
if IAMSCOTT
    imslice(abs(slidingWindowReconVol),sprintf('imgspace %i',volfig_id));
    profile viewer
else
    disp_vol_center(slidingWindowReconVol,0,volfig_id);
    outfile=sprintf('%s_%s%05.2f.nii',reconDir,reco_speed,oversampling);
    profout=sprintf('%s_%s%05.2f.profile',reconDir,reco_speed,oversampling);
    fprintf('saving %s\n',outfile);
    save_nii(make_nii(abs(slidingWindowReconVol)),outfile);
    profsave(profile('info'),profout);
end
