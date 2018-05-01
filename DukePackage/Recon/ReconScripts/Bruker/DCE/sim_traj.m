function sim_traj(trajFile,endpointsFile,nPts,rayCount,trajShift,trajMin)
% function sim_traj(trajFile,endpointsFile,nPts,rayCount,garbagePts)
% simulate trajectories given required variables. 
%
% trajFile      - destination traj file
% endpointsFile - input endpoints file
% nPts          - complex points per ray
% rayCount      - number of rays we'll simulate, input file must have at 
%                 least this N points. 
% trajShift     - shift resulting trajectory this nPts. Normally a negative
%                 number equivalent to the bad pts per line. eg, if we 
%                 have 2 garbage pts per line, this would be -2.
% trajMin       - value from 0-1 setting the lowest number for traj
%                 calc,(instead of staring at 0)
% 
% Saves output while it traverses the input file. 
%
% This stupidly loops when some smart matrix operations could have been
% done. 
% Oh well -\_(O_o)_/-
if ~exist('trajShift','var')
    % How late does acq window start after gradients. 
    % Ideally acq window would start first, but this accounts for that. 
    trajShift=0;
end
if ~exist('trajMin','var')
    trajMin=0;
elseif trajMin>0.25  || trajMin<0
    error('invalid trajMin given, should be between 0-0.25');
end
%%
% happy little gradient multiplier and plot code from our dear friend scott. 
figure(2000);
% Plot gradients (g) vs time (t)
t = 1:nPts;
% how many ramp points do we have in our acquisition window. 
% This was set in the acquisition pulse sequence.
ramp_pts=round(nPts/3);
% ramp points seemed valid from round(nPts/3) (+/-) 2
g = [linspace(0,1,ramp_pts) ones(1,nPts-ramp_pts)];
if trajMin>0
    g=g+trajMin;
    g(g>1)=1;
end
if trajShift>0
    g=circshift(g,[0,-trajShift]);
    g(end-trajShift:end)=1;
elseif trajShift<0
    g=circshift(g,[0,-trajShift]);
    g(1:-trajShift)=0;
end
plot(t,g);
hold on
ylim([-0.05 1.1])
xlabel('sample number');
ylabel('gradients');

% Plot trajectory (traj) vs time (t)
traj = cumtrapz(g);
traj = 0.499999*traj/max(traj);% normalize values to +/- 0.499999
plot(t,traj);
legend('Gradients','Trajectory');
%%
hfid=fopen(endpointsFile,'r');
tfid=fopen(trajFile,'w');
% read 3x rayCount doubles representing each end pt,
% then fill inthe blanks and save
traj_p=fread(hfid,3*rayCount,'double',0,'l');
traj_p=reshape(traj_p,[3,rayCount]);
g_delay=[0.0,0.0,0.0];% delays for each gradient. 
fprintf('Preparing full trajectory for %i rays with %i pts each\n',rayCount,nPts);
fprintf('This could be a while : ( \n');
for tn=1:rayCount
    dx=traj_p(:,tn);
    l=zeros(3,nPts);
    for gx=1:3
        if dx(gx)~=0
            % linear
            % l(gx,:)=(0:(dx(gx)/(nPts-1)):traj_p(gx,tn));
            % using grad calced version
            % G_delay doesnt work right becuase it increases the length
            % of vector. I should normalize the vector to its max to fix.
            % max(traj_p(gx,tn))
            l(gx,:)=traj*(traj_p(gx,tn)-g_delay(gx));
        end
    end
    fwrite(tfid,l,'double',0,'l');
end
fclose(hfid);fclose(tfid);