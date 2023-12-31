set(groot,'defaulttextinterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultAxesFontsize',40);
if ~exist('OutputFigures')
    mkdir('OutputFigures')
end
OutputFolder='OutputFigures';
FigPreamble='SimulatedData_CommonPolodial';

%----data params----
DataParams=struct;
DataParams.NumberOfDatapoints=1e3;
DataParams.Layout='Torus_Torodial';
DataParams.R1=10;
DataParams.r1=3;
DataParams.R2=10;
DataParams.r2=5;

%----kernels construnction params----
NormFac=0.1;
NormFacVisualization=0.3;
% 0<NormFac<1 - global scale factor that eqauls to NormFac*(median of the pairwise distances)
% NormFac>1 - local scale factor thar equals to the median of the distances of the NormFac nearest neighbours
UseMutualScaleFac='None';
% UseMutualScaleFac='Mean' - use the same scale factor for each kernel, the mutual scale factor is set to be their mean
% UseMutualScaleFac='Min' - use the same scale factor for each kernel, the mutual scale factor is set to be the lower scale factor
% UseMutualScaleFac='None' - use different scale factor for each kernel
KernelsParams=v2struct(NormFac,NormFacVisualization,UseMutualScaleFac);

%----eigenvalues flow diagram params----
NumberOfEigenVals=20; %Number of eigenvalues calculated at each point on the geodesic
ntVec=200;% Number of points on the geodesic grid
TolFac=1e-5;% Tolerance factor for fixed rank approximation
Interpolator='Geodesic';
% Interpolator='Linear' - linear interpolation: (1-t)*K1+t*K2
% Interpolator='Geodesic' - geodesic interpolation: K1^(-1/2)*(K1^(1/2)*K2^(-1)*K1^(1/2))^t
% Interpolator='Harmonic' - Harmonic interpolation: K1^(1-t)*K2^t
EvfdParams=v2struct(NumberOfEigenVals,TolFac,ntVec);
