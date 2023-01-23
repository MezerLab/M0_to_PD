function [opt]=M0toPD_CSFnorm(opt);

% This function allows for normalization by the median of the CSF
%
% Inputs:
%       - opt - the structure for all previous analyses
% Outputs:
%       - opt - the same structure
%

PD = readFileNifti(opt.PDfile);
PD = PD.data;

BM=readFileNifti(opt.BMfile);
xform=BM.qto_xyz;
BM=BM.data;

seg = readFileNifti(opt.segfile);
seg = seg.data;

CSF = median(PD(BM==1&seg==0));
PD = PD./CSF;

opt.WF_file = fullfile(opt.outDir, 'WF.nii.gz')
dtiWriteNiftiWrapper(PD,xform,opt.WF_file);
end

