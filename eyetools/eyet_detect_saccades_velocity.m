function [movement] = eyet_detect_saccades_velocity(data,param)

% FROM FT_SACCADE_DETECTION performs micro/saccade detection on time series data
% over multiple trials
%
% The configuration should contain:
%  param.method   = different methods of detecting different movement types
%                , Micro/saccade detection based on Engbert R,
%                   Kliegl R (2003) Vision Res 43:1035-1045. The method
%                   computes thresholds based on velocity changes from
%                   eyetracker data (horizontal and vertical components).
%                'clustering', Micro/saccade detection based on
%                   Otero-Millan et al., (2014) J Vis 14 (not implemented
%                   yet)
%
% METHOD SPECIFIC OPTIONS AND DESCRIPTIONS
%
%  VELOCITY2D
%   VELOCITY2D detects micro/saccades using a two-dimensional (2D) velocity
%   space velocity. The vertical and the horizontal eyetracker time series
%   (one eye) are transformed into velocities and microsaccades are
%   indentified as "outlier" eye movements that exceed a given velocity and
%   duration threshold.
%     param.velocity2D.kernel   = vector 1 x nsamples, kernel to compute velocity (default = [1 1 0 -1 -1].*(data.fsample/6);
%     param.velocity2D.demean   = 'no' or 'yes', whether to apply centering correction (default = 'yes')
%     param.velocity2D.mindur   = minimum microsaccade durantion in samples (default = 3);
%     param.velocity2D.velthres = threshold for velocity outlier detection (default = 6);
%
% The output argument "movement" is a Nx3 matrix. The first and second
% columns specify the begining and end samples of a movement period
% (saccade, joystic...), and the third column contains the peak
% velocity/acceleration movement. This last thrid column will allow to
% convert movements into spike data representation, making the spike
% toolbox functions compatible (not implemented yet).
%

fsample = param.FS;

% set the defaults for the various microsaccade detection methods
% Engbert R, Kliegl R (2003) Microsaccades uncover the orientation of
% covert attention. Vision Res 43:1035-1045.
kernel = [1 1 0 -1 -1].*(fsample/6); % this is equivalent to Engbert et al (2003) Vis Res, eqn. (1)
if ~isfield(param.velocity2D, 'kernel'),   param.velocity2D.kernel  = kernel; end
if ~isfield(param.velocity2D, 'demean'),   param.velocity2D.demean  = 'yes';  end
if ~isfield(param.velocity2D, 'mindur'),   param.velocity2D.mindur  =  3;     end % minimum microsaccade duration in samples
if ~isfield(param.velocity2D, 'velthres'), param.velocity2D.velthres = 6;     end

% determine the size of the data
ntrial = size(data,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the actual computation is done in the middle part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
movement = [];
% do all the computations
for i=1:ntrial
    fprintf('finding microsaccades trial %d of %d\n', i, ntrial);
    
    dat = data(i,:);
    ndatsample = size(dat,2);
    
    % demean horizontal and vertical time courses
    if strcmp(param.velocity2D.demean, 'yes');
        dat = ft_preproc_polyremoval(dat, 0, 1, ndatsample);
    end
    
    %% eye velocity computation
    % deal with padding
    n = size(param.velocity2D.kernel,2);
    pad = ceil(n/2);
    dat = ft_preproc_padding(dat, 'localmean', pad);
    
    % convolution. See Engbert et al (2003) Vis Res, eqn. (1)
    if n<100
        % heuristic: for large kernel the convolution is faster when done along
        % the columns, weighing against the costs of doing the transposition.
        % the threshold of 100 is a bit ad hoc.
        vel = convn(dat,   param.velocity2D.kernel,   'same');
    else
        vel = convn(dat.', param.velocity2D.kernel.', 'same').';
    end
    % cut the eges
    vel = ft_preproc_padding(vel, 'remove', pad);
    
    %% microsaccade detection
    % compute velocity thresholds as in Engbert et al (2003) Vis Res, eqn. (2)
    medianstd = sqrt( nanmedian(vel.^2,2) - (nanmedian(vel,2)).^2 );
    
    % Engbert et al (2003) Vis Res, eqn. (3)
    radius = param.velocity2D.velthres*medianstd;
    
    % compute test criterion: ellipse equation
    test = sum((vel./radius(:,ones(1,ndatsample))).^2,1);
    sacsmp = find(test>1);% microsaccade's indexing
    
    %% determine microsaccades per trial
    % first find eye movements of n-consecutive time points
    j = find(diff(sacsmp)==1);
    j1 = [j; j+1];
    com = intersect(j,j+1);
    cut = ~ismember(j1,com);
    sacidx = reshape(j1(cut),2,[]);
    
    for k=1:size(sacidx,2);
        duration = sacidx(1,k):sacidx(2,k);
        if size(duration,2) >= param.velocity2D.mindur;
            % finding peak velocity by Pitagoras
            begtrl = sacsmp(duration(1,1));
            endtrl = sacsmp(duration(1,end));
            
            [peakvel smptrl] = max(sqrt(sum(vel(:,begtrl:endtrl).^2,1)));
            veltrl = sacsmp(duration(1,smptrl));% peak velocity microsaccade sample -> important for spike conversion
            
%             trlsmp = data.sampleinfo(i,1):data.sampleinfo(i,2);
%             begsample = trlsmp(1, begtrl); % begining microsaccade sample
%             endsample = trlsmp(1, endtrl); % end microsaccade sample
%             velsample = trlsmp(1, veltrl); % velocity peak microsaccade sample
            movement(end+1,:) = [begtrl endtrl veltrl];
        end
    end
    
end
