function [logname]=PD_Fit_saveParams(outDir,degrees,M0file,Qmapfile,BMfile,...
                                            SEGfile,outMm,boxSize,... 
                                            percent_overlap,Inclusion_Criteria,...
                                            Qmap_factor)
                                        
% This function creates a structure of information to fit the M0 boxes for
% PD.
%
%   ~INPUTS~
%             outDir:   The output directory. 
%                           The function also reads the file from there.
%            degrees:   Polynomial degrees for the coil estimation.
%             M0file:   The combined/aligned M0 data.
%             Qmapfile: The map.
%             BMfile:   The defined brain mask.
%             SEGfile:  Segmentation file, where 0 is not brain, and each
%                       tissue type has a different label (typically WM, GM,
%                       CTX and CSF).
%              outMm:   The resample (undersample) resolution of those 
%                           images. (Default is 2mm x 2mm x 2mm.) This can
%                           shorten the fits. The fit procedure was written
%                           to this resolution. This resolution is assumed
%                           to be high, given the low frequency change of
%                           the coils gain.
%            boxSize:   The box size that is used for the fit (in mm)
%                           Default is 14.
%    percent_overlap:   The overlap between the boxes. 
%                           (Default is 0.5 --> 50%)
% Inclusion_Criteria:   Default is [0.8 200].
%        Qmap_factor:   Default is 1.
%
%     ~OUTPUTS~
%            logname:   A structure that saves the fit parameters. It is 
%                           saved in the outDir, and in the tmp directory,
%                           with the name fitLog.mat. 
%
%
%
% AM (C) Stanford University, VISTA
%
%

%% I. Check inputs and set defaults
%  Saving parameters and relevant information for the Gain fit in the "opt"
%  structure. This allows us to send them to all the grid calls running 
%  in parallel

if (notDefined('outDir') || ~exist(outDir,'dir'))
    outDir = uigetDir(pwd,'Select outDir');
end
   opt.outDir  = outDir;
   opt.PDfit_Method=1;
    Coilsinfo.maxCoil=1; Coilsinfo.minCoil=1; Coilsinfo.useCoil=1; 
    opt.maxCoil=Coilsinfo.maxCoil;
    opt.minCoil=Coilsinfo.minCoil;
    opt.useCoil=Coilsinfo.useCoil;
    
    if (~exist('Qmap_factor','var'))
        Qmap_factor = 1;
    end
        opt.Qmap_factor = Qmap_factor;

if(~exist('degrees','var'))
    disp('Using the default polynomials: degrees = 3 for coil estimation');
    degrees = 3;
end
    opt.degrees = degrees;
    opt.M0file = M0file;
    opt.T1file=Qmapfile;

if notDefined('boxSize')
    boxSize =14;
end

if notDefined('percent_overlap')
    percent_overlap =0.5;
end

% In cases of high-resolution data, we can undersample the data to outMm
% resolution
if notDefined('outMm')
    outMm=[2 2 2];
end

 

if (notDefined('Inclusion_Criteria'))
    opt.Inclusion_Criteria=[0.8 200];
else
    opt.Inclusion_Criteria=Inclusion_Criteria;
    
end

if (exist(BMfile,'file'))
    disp(['Loading brain Mask data from ' BMfile '...']);
    brainMask = readFileNifti(BMfile);
    mmPerVox  = brainMask.pixdim;
    opt.BMfile = BMfile;
    
    % In cases of high resolution data, we can undersample the data when we
    % fit the coil gains. This may shorten the fit time. Note that the coil
    % gain was written for 2mm x 2mm x 2mm resolution. Images with a
    % different resolution may need some changes in the regularization
    % protocol.
    
    if outMm(1)~=0 % Unless we don't want to change the fitting resolution
       if (outMm(1)~=mmPerVox(1) || outMm(2)~=mmPerVox(2) || outMm(3)~=mmPerVox(3) )
          [opt]=mrQ_resamp4G_fit(opt,outMm);
           brainMask = readFileNifti(opt.BMfile);
           mmPerVox  = brainMask.pixdim;
       end
    end
    brainMask = logical(brainMask.data);
else
    error('Cannot find the file: %s', BMfile);
end

% segmenting R1 to tissue types
    opt.segfile=SEGfile;
 
%% II. Identify the boxes we want to fit

% Try to parcel the brain into boxes of roughly boxSize mm^3 and with an
% overlap of 2;
sz   = (size(brainMask));

boxS = round(boxSize./mmPerVox);
even = find(mod(boxS,2)==0);

boxS(even)  = boxS(even)+1;
opt.boxS = boxS;

% Determine the percentage of percent_overlap  (0.1, 0.5, 0.7)
overlap = round(boxS.*percent_overlap);

% Create a grid of the center of the boxes that will be used to fit
[opt.X,opt.Y,opt.Z] = meshgrid(round(boxS(1)./2):boxS(1)-overlap(1):sz(1),...
                        round(boxS(2)./2):boxS(2)-overlap(2):sz(2),...
                        round(boxS(3)./2):boxS(3)-overlap(3):sz(3));


%donemask is a volume of the center locations and is used for bookkeeping. 
%(i.e., which boxes are done and which need be done or skipped.)
donemask = zeros(size(opt.X));

opt.HboxS = (boxS-1)/2;

%% III. Loop over the identified boxes and check that the data is there

ii = 1;
opt.donemask = donemask;

for i=1:prod(size(opt.X))
    
    [fb(1) fb(2) fb(3)] = ind2sub(size(opt.X),i);
    [empty] = M0toPD_isDataBox(opt,brainMask,fb,opt.Inclusion_Criteria);
    if empty == 0 % this is a good box
        opt.wh(ii) = i;
        ii = ii+1;
        
    elseif empty == 1 %box we won't use
        donemask(fb(1),fb(2),fb(3)) = -1e3;
        
    elseif empty == -1 %box we won't use
        donemask(fb(1),fb(2),fb(3)) = -2e3;
        opt.wh(ii) = i;
        ii = ii+1;
    end
    
end

opt.donemask = donemask;

%% IV. Initiate other parameters for the fit and SGE call

opt.jumpindex = length(opt.wh);
opt.Kfold=3;   % the fold for cross validation (use split half)
opt.BasisFlag = 'qr'; %orthonormal basis for the coil polynomials

opt.smoothkernel=0;
dirname    = [outDir '/tmpSGM0' ];
dirDatname = [outDir '/tmpSGM0dat'];

opt.dirDatname = dirDatname;
opt.name = [dirname '/M0boxfit_iter'] ;
opt.date = date;
opt.dirname=dirname;

% Save out a logfile with all the options used during processing
logname = [outDir '/fitLog.mat'];
opt.logname=logname;

%saving an information file we can load afterwards if needed
save(opt.logname,'opt');

if ~isfolder(dirname)
    mkdir(dirname);
end
return
