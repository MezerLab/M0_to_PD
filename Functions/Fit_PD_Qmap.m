function Fit_PD_Qmap(opt)
%
% INPUTS:
%         opt:   This is the optimization structure that was passed from
%                   mrQ_fitPD_multicoil. It has all the needed information.
% OUTPUTS:
%                The function will save an output file with fitted
%                parameters in a tmp directory. This will be used later by
%                mrQfitPD_multiCoils_M0 to make the PD map.
%
% AM (C) Stanford University, VISTA
%
%

%% I. Initialization

%Find the box to work on
j=0;
st=1 ;
ed=length(opt.wh);

nIteration=ed-st+1;
%Initialize the parameters and saved outputs

% Get the M0 and Qmap information

% Multi coil M0
M0=readFileNifti(opt.M0file); % PDw
M0=M0.data;

%T1
Tmap=readFileNifti(opt.T1file); %T1w
Tmap=Tmap.data;

%Brain mask
BM=readFileNifti(opt.BMfile);
BM=BM.data;

BM(Tmap>3000)=0; % Clear areas that are not GM or WM.

%seg mask
seg=readFileNifti(opt.segfile);
seg=seg.data;

if size(seg)~=size(Tmap)
    disp('Please choose another segmentation file!')
    segfile = mrvSelectFile('r','*.nii.gz','Select segmentation file',opt.segfile);
    seg = readFileNifti(segfile);
    seg = seg.data;
end

smoothkernel=opt.smoothkernel;

% The poly basis to fit the coil gains
pBasis = mrQ_CreatePoly(opt.boxS,opt.degrees,3,opt.BasisFlag);

nVoxels=size(pBasis,1);
nPolyCoef=size(pBasis,2);

% Initiate the saved parameters
fb=zeros(nIteration,1,3);
gEst=zeros(nPolyCoef,nIteration);
resnorm=zeros(nIteration,1);
exitflag=zeros(nIteration,1);
skip=zeros(nIteration,1);

Iter=0;
Qmap_factor = opt.Qmap_factor;

%%  II. Go over it, box by box
kk=0;
for jj= st:ed
    %run over the box you like to fit
    clear M01  qmap  BM1  SZ M0_v Qmapbasis PDinit Segmask g0 G0 mask1
    Iter= Iter+1;
    tic
    %Find the x,y,z location of the box.
    % (This is not the x,y,z location in image space, but rather the grid
    % of boxes we made by meshgrid in  mrQ_PD_multicoil_RgXv_GridCall.m)
    [fb(Iter,1,1), fb(Iter,1,2), fb(Iter,1,3)]=ind2sub(size(opt.X),opt.wh(jj));
    
    % Get all the relevant box data for the fit
    [M01, tmap, BM1, SZ, skip(Iter)]= mrQ_GetM0_boxData(opt,Tmap,M0,BM,fb(Iter,1,:),smoothkernel,seg);
    tmap = tmap.*Qmap_factor;
    
    
    
    if  skip(Iter)==1
       % disp(['skipping box ' num2str(jj) ' bad data'])
        
    else
        
    %% Fit
        %% Initiate the multi-box fit parameters  1/PD=A +B/T1
        A=zeros(1,7);B=zeros(1,7);
        rmap = 1./tmap;
        
        [a,b] = find_A_B_coeffs(M01,rmap,BM1);
        A(1,1) = a(1); B(1,1) = b(1); 

        % a basis for estimating the new A and B.
        Qmapbasis(:,2)=rmap(BM1);
        Qmapbasis(:,1)=1;
        
        %iterate to converge A and B
        for ii=2:7
            PDp=1./(A(ii-1)+B(ii-1)./tmap);

                      %  PDp=PDp./median(PDp(BM1)); scale

            % the sensitivity the receive profile
            RPp=M01./PDp;
            
            % Raw estimate
            g = pBasis(BM1,:) \ RPp(BM1);  % Polynomial approximation
            RPi=pBasis*g;
            
            % calculate PD from M0 and RP
            PDi=M01(:)./RPi;
            
            % solve for A B given the new PD estimation
            % ( 1./PDi(BM1) )= A* R1basis(:,1) + B*R1basis(:,2);
            
            co     = Qmapbasis \ ( 1./PDi(BM1) );
            A(ii)=co(1);
            B(ii)=co(2);
            
        end
        gEst(:,Iter)=g;
        
    end
end

name=[ opt.name '_' num2str(st) '_' num2str(ed)];

save(name,'gEst','st','ed','skip','fb')

