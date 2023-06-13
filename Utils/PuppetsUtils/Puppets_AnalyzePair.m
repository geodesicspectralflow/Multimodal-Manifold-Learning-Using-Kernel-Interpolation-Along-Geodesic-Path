disp(sprintf('\n--------- Generating figures for %s --------- ',FigPreamble))

%% Compute kernels
d1 = pdist2(s1,s1,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 1.
d2 = pdist2(s2,s2,'euclidean');             % The distance between samples, based only on snapshots taken by Camera 2.
D1=sqrt(d1);
D2=sqrt(d2);

UseMutualEps='none';
eps1=median(D1(:))*NormFac;
eps2=median(D2(:))*NormFac;

switch UseMutualEps
    case 'Min'
        eps1=min([eps1,eps2]);eps2=eps1;
    case 'Mean'
        eps1=(eps1+eps2)/2;eps2=eps1;
end

A1=exp(-D1.^2/(eps1^2));
[ColumnStochasticK1,K1] =SingleIterationNorm(A1,1);
A2=exp(-D2.^2/(eps2^2));
[ColumnStochasticK2,K2] =SingleIterationNorm(A2,1);

%% Compute EVFDs
GT=MapEmbdbulldog(:,2)';
cX=GT-mean(GT);

Lmax=100;
PolyNormFac=Inf;
Dim=GetEffectiveDim(K1,K2,TolFac);
tVec=linspace(0,1,ntVec);

Linear_KernelEigenValuesMat=[] ;Linear_KernelEigenValuesCorrMat= [];
Geodesic_KernelEigenValuesMat=[] ;Geodesic_KernelEigenValuesCorrMat=[];
LogEuc_KernelEigenValuesMat=[] ;LogEuc_KernelEigenValuesCorrMat=[];
BuresWass_KernelEigenValuesMat=[] ;BuresWass_KernelEigenValuesCorrMat=[];

tmpBW=(K1*K2)^(1/2)+(K2*K1)^(1/2);
logK1=logm(K1); logK2=logm(K2);


reverseStr = '';
tVec=linspace(0,1,ntVec);
disp('Calculating the eigenvalues flow diagram');
for tind=1:length(tVec)
    %     tind
    percentDone = 100 * tind / length(tVec);
    msg = sprintf('   percentage done: %3.1f', percentDone); %Don't forget this semicolon
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
    t=tVec(tind);
    K_gt=FixedGeodes( K1,K2,t,Dim );K_gt=real(K_gt);
    K_lt=(1-t)*K1+t*K2;
    K_le=expm((1-t)*logK1+t*logK2);
    K_bw=(1-t)^2*K1+t^2*K2+t*(1-t)*tmpBW;

    
    [Vs,EigenVals]=GetSortedEVs(K_lt,NumberOfEigenVals);
    Linear_KernelEigenValuesMat=[Linear_KernelEigenValuesMat;EigenVals];
    Linear_KernelEigenValuesCorrMat=[Linear_KernelEigenValuesCorrMat;abs(corr(GT',Vs(:,1:NumberOfEigenVals)))];
    
    [Vs,EigenVals]=GetSortedEVs(K_gt,NumberOfEigenVals);
    Geodesic_KernelEigenValuesMat=[Geodesic_KernelEigenValuesMat;EigenVals];
    Geodesic_KernelEigenValuesCorrMat=[Geodesic_KernelEigenValuesCorrMat;abs(corr(GT',Vs(:,1:NumberOfEigenVals)))];
    
    [Vs,EigenVals]=GetSortedEVs(K_bw,NumberOfEigenVals);
    BuresWass_KernelEigenValuesMat=[BuresWass_KernelEigenValuesMat;EigenVals];
    BuresWass_KernelEigenValuesCorrMat=[BuresWass_KernelEigenValuesCorrMat;abs(corr(GT',Vs(:,1:NumberOfEigenVals)))];


    [Vs,EigenVals]=GetSortedEVs(K_le,NumberOfEigenVals);
    LogEuc_KernelEigenValuesMat=[LogEuc_KernelEigenValuesMat;EigenVals];
    LogEuc_KernelEigenValuesCorrMat=[LogEuc_KernelEigenValuesCorrMat;abs(corr(GT',Vs(:,1:NumberOfEigenVals)))];
end
fprintf([reverseStr]);

%% AD
K_ad=0.5*( GetCS(K1)* GetCS(K2)+ GetCS(K2)* GetCS(K1));%K=K^2;

%% Show EVFDs side by side
[Linear_EigenValuesTupple,Linear_CorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,Linear_KernelEigenValuesMat,Linear_KernelEigenValuesCorrMat,tVec);
[Geodesic_EigenValuesTupple,Geodesic_CorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,Geodesic_KernelEigenValuesMat,Geodesic_KernelEigenValuesCorrMat,tVec);
[BuresWass_EigenValuesTupple,BuresWass_CorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,BuresWass_KernelEigenValuesMat,BuresWass_KernelEigenValuesCorrMat,tVec);
[LogEuc_EigenValuesTupple,LogEuc_CorrEigenValuesRawStack] = GetTupleForScatter(NumberOfEigenVals,LogEuc_KernelEigenValuesMat,LogEuc_KernelEigenValuesCorrMat,tVec);

figure();
subplot(1,4,1);
scatter((Geodesic_EigenValuesTupple(:,1)),Geodesic_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Geodesic($\gamma(t)$)','FontSize', 15);caxis([0 1]);colormap jet;
subplot(1,4,2);
scatter((Linear_EigenValuesTupple(:,1)),Linear_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Linear($\mathbf{L}(t))$','FontSize', 15);caxis([0 1]);colormap jet;

subplot(1,4,3);
scatter((BuresWass_EigenValuesTupple(:,1)),BuresWass_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Bures-Wasserstein($\gamma^{BW}(t)$)','FontSize', 15);caxis([0 1]);colormap jet;

subplot(1,4,4);
scatter((LogEuc_EigenValuesTupple(:,1)),LogEuc_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
title('Log-Euclidean($\gamma^{LE}(t)$)','FontSize', 15);caxis([0 1]);colormap jet;

sgtitle(sprintf('Eigenvalues flow diagram (EVFD)'),'FontSize', 30);
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFDsComparison_SideBySide");

%% Overlay EVFDs
figure();
subplot(1,3,1);
scatter((Geodesic_EigenValuesTupple(:,1)),Geodesic_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
hold on;
scatter((Linear_EigenValuesTupple(:,1)),Linear_EigenValuesTupple(:,2),50,'r','filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
legend('$\gamma(t)$','$\mathbf{L}(t))$');

subplot(1,3,2);
scatter((Geodesic_EigenValuesTupple(:,1)),Geodesic_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
hold on;
scatter((BuresWass_EigenValuesTupple(:,1)),BuresWass_EigenValuesTupple(:,2),50,'r','filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
legend('$\gamma(t)$','$\gamma^{BW}(t)$');

subplot(1,3,3);
scatter((Geodesic_EigenValuesTupple(:,1)),Geodesic_EigenValuesTupple(:,2),50,'b','filled','o');
ylabel('$t$');xlabel('$\log(\mu_{t}^{k})$');axis tight;
hold on;
scatter((LogEuc_EigenValuesTupple(:,1)),LogEuc_EigenValuesTupple(:,2),50,'r','filled','o','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);
legend('$\gamma(t)$','$\gamma^{LE}(t)$');

sgtitle(sprintf(['A comparison between the EVFD obtained using a geodesic interpolation ($\\gamma(t)$) \n' ...
    'and the EVFD obtained using linear, Bures-Wasserstein and Log-Euclidean interpolation ($\\mathbf{L}(t)),\\gamma^{BW}(t),\\gamma^{LE}(t)$)']),'FontSize', 15);
SaveFig(gcf,OutputFolder,FigPreamble+"_EVFDsComparison_Overlay");

%% Show insets
FontSize=10;xlFontSize=60;ylFontSize=60;
PointSize=100;

for t0=[0,0.5,1]
    K=FixedGeodes( K1,K2,t0,Dim );K=real(K);
    [MapEmbd, ~, ~,singvals, ~,~] =DiffusionMapsFromKer( K , 1);
    
    figure();
    subplot(3,3,1);
    scatter(AngleYoda,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none');
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,2);
    scatter(AngleBulldog,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$',t0));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,3);
    scatter(AngleBunny,MapEmbd(:,2),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{2}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    
    %
    subplot(3,3,4);
    scatter(AngleYoda,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{2}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,5);
    scatter(AngleBulldog,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{3}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,6);
    scatter(AngleBunny,MapEmbd(:,3),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{3}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{3}_{t^{*}}$'));
    end
    ax.YLabel.Visible = 'on';
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    %
    subplot(3,3,7);
    scatter(AngleYoda,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$'));
    end
    %     if t0==0
    xlabel('$\theta_{\mathrm{yoda}}$');
    %     end
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    axis off;
    
    subplot(3,3,8);
    scatter(AngleBulldog,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$',t0));
    end
    xlabel('$\theta_{\mathrm{bulldog}}$');
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    subplot(3,3,9);
    scatter(AngleBunny,MapEmbd(:,4),'filled');axis tight
    ax = gca; axis(ax,'off');
    ylabel(sprintf('$v^{4}_{%g}$',t0));
    if t0==0.5
        ylabel(sprintf('$v^{4}_{t^{*}}$'));
    end
    xlabel('$\theta_{\mathrm{bunny}}$');
    ax.YLabel.Visible = 'on';
    ax.XLabel.Visible = 'on';
    
    set(gca,'color','none')
    set(get(gca,'XAxis'),'FontSize', FontSize);    set(get(gca,'YAxis'),'FontSize', FontSize);
    set(get(gca,'XLabel'),'FontSize', xlFontSize);    set(get(gca,'YLabel'),'FontSize', ylFontSize);
    
    SaveFig(gcf,OutputFolder,[FigPreamble,sprintf('_Inset_%g',t0)]);
end