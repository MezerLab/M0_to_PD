function [Boxes, PositiveBoxs, UnCorBoxs, UnSTDBoxs]=CalBoxPD_step1a(opt,BoxesToUse,CoilGains)

% This is step 1 of 6 for building the WF map.
%     Step 1: Get the x, y, z, and PD values of the boxes. 
%
% ~INPUTS~ 
%          opt:
%   BoxesToUse:
%    CoilGains:
%
% ~OUTPUTS~
%        Boxes:
% PositiveBoxs:
%    UnCorBoxs:
%    UnSTDBoxs:
%
% See also: mrQ_buildPD_ver2
%           Step_0: none
%           Step_2: mrQ_ScaleBoxes_step2
%           Step_3: mrQ_BoxJoinBox
%           Step_4: mrQ_smoothGain_step4b
%           Step_5: mrQ_PD2WF_step5
%
% AM Vistalab team 2013


%% I. Book keeping

% bookkeeping of the box fits that went wrong
PositiveBoxs=zeros(length(opt.wh),1);
UnCorBoxs=zeros(length(opt.wh),1);
UnSTDBoxs=zeros(length(opt.wh),1);

%% II. Get M0, T1 and Brain Mask 
% Get the M0 and T1 information
% multi coil M0
M0=readFileNifti(opt.M0file);
M0=M0.data;
SZ=size(M0);
%T1
%T1=readFileNifti(opt.T1file);
%T1=T1.data;

%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

smoothkernel=opt.smoothkernel;
pBasis = M0toPD_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);

nPolyCoef=size(pBasis,2);
nVoxels=size(pBasis,1);

for ii=BoxesToUse
    clear M01  t1  BM1  SZ M0_v R1basis PD
    
    [fb(1,1,1), fb(1,1,2), fb(1,1,3)]=ind2sub(size(opt.X),opt.wh(ii));
    
    % get all the relevant box data for the fit
    [M01, ~, ~, ~, ~, ~,~, XX, YY, ZZ ]= M0toPD_GetM0_boxData(opt,[],M0,BM,fb(1,1,:),smoothkernel);
    M0_v=M01(:);
        

%% III. G, PD, coil gain and coil PD
    % 1. Get G
    g=CoilGains(ii).g;
    G = pBasis*g;
    
    % 2. Solve PD
    Clist=CoilGains(ii).Clist;
    PD = zeros(nVoxels,1);
    for jj=1:nVoxels
        PD(jj) = G(jj,:)' \ M0_v(jj,Clist)';
    end
    
% If there are zeroes, we will get NaN. 
% If it less then zero, it is just wrong
    mask=PD>0; 
    
% Check that the between-coil error is minimal
    PDC=M0_v(:,Clist)./G;
    EstimateErr=std(PDC./repmat(mean(PDC,2),1,size(PDC,2)),[],2);
    mask=mask & EstimateErr<0.08;
    
           PDSTD=std(PD(mask));
           
% Check that there are no outlier voxels
    mask=mask & PD<PD+3*PDSTD & PD>PD-3*PDSTD ;
    
  
    % 3. Solve for all coils Gain
    G  = zeros(nVoxels,1);
    g0 = zeros(nPolyCoef,1);
    G(mask)  = M0_v(mask,:) ./ PD(mask);         % Raw estimate
    g0(:) = pBasis(mask,:) \ G(mask,:);  % Polynomial approximation
    G = pBasis*g0;
    
    % 4. Solve all coils PD
    PD = zeros(nVoxels,1);
    V=1:nVoxels & mask';
    for jj=find(V)
        PD(jj) = G(jj,:)' \ M0_v(jj,:)';
    end
    
    Boxes(ii).PD=PD;
    Boxes(ii).XX=XX;
    Boxes(ii).YY=YY;
    Boxes(ii).ZZ=ZZ;
 
%% IV. Various Checks    
% Check for correlation role
    Bad=[];
    
% Check for bad solution negative PD. 
    if     length(find(mask))/length(mask)<0.5 
%        If there are 50% negative values, it is clearly a wrong solution
         Boxes(ii).NegativeBad=1;
       
         
    else
         Boxes(ii).NegativeBad=0;
         PositiveBoxs(ii)=1;
         
         mask1=zeros(size(BM));
         mask1(XX(1):XX(2),YY(1):YY(2),ZZ(1):ZZ(2))=1;
         wh=find(mask1);
         Boxes(ii).loc=wh(mask);
          Boxes(ii).PD=PD(mask);
          
    end
    
% Check if the M0 vector is more correlated than the G. 
    % They must be or the solution is wrong.
    if isempty(Bad)
        Boxes(ii).Corgood=1;
        UnCorBoxs(ii)=1;
    else
        Boxes(ii).Corgood=0;
    end
    Boxes(ii).CorBad=Bad;
    
    
end
