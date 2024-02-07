clear; close all; clc
addpath([pwd '\functions\']);

%% parameters
M=1; C=1; % univariate analysis
b=0.01:0.01:1; % internal dynamics parameter -- it weights the pole radius
l=length(b);
nfft=1000; % number of points on frequency axis (positive)
q=20;
fs=1; % sampling freq.
z = 1i*2*pi/fs; % complex variable normalized by fs
rangeLF=[0.04 0.15];
rangeHF=[0.15 0.4];
range_all=[rangeLF; rangeHF];
delta_LF=rangeLF(2)-rangeLF(1);
delta_HF=rangeHF(2)-rangeHF(1);
bandLF=round((nfft*2/fs)*rangeLF);
bandHF=round((nfft*2/fs)*rangeHF);

DimFont=24;
axislinewidth=1;
grey_line=[0.5 0.5 0.5];
% colormap for coefficient b
colrange=[0 0 1];
lim1=floor(l/2);
for i=1:lim1
    colrange=[colrange; 0 i/lim1 (lim1-i)/lim1]; %#ok
end 
for i=1:l-lim1-1
    colrange=[colrange; i/(l-lim1-1) (l-lim1-1-i)/(l-lim1-1) 0]; %#ok
end 

% init
Sm=cell(1,l); % spectrum
LFspect=zeros(nfft,l); % LF spectrum
HFspect=LFspect; % HF spectrum
LSP=zeros(1,l); % time domain LSP
LSPf_LF=LFspect; LSPf_HF=LFspect; % total term in a given band
H_HF=LFspect; H_LF=LFspect; % transfer functions relative to specific oscillations
H=LFspect; % total transfer function
lspf_HF=LFspect; lspf_LF=LFspect; % freq specific term of the spectral LSP integrated in band
LSPf_integral_HF=LSP;
LSPf_integral_LF=LSP;

%% theoretical process 
for il=1:l

    par.poles{1}=([b(il)*0.8 0.1; 0.9 0.3]); 
    par.Su=1; 
    par.coup=[];
    
    % theoretical VAR parameters
    [Am,Su,Ak]=theoreticalVAR(M,par); 
    Am=Am';

    % compute LSP
    out=lsp_computation(Am,Su,range_all,q,nfft,fs);

    LSP(il)=out.LSP; % time domain
    LSP_tmp=LSP(il).*ones(1,nfft);
    LSP_HF=mean(LSP_tmp(bandHF(1):bandHF(2)))*delta_HF;
    LSP_LF=mean(LSP_tmp(bandLF(1):bandLF(2)))*delta_LF;

    index_LF = out.index_band(1);
    index_HF = out.index_band(2);

    % freq specific term
    lspf_all=out.lspf; 
    lspf_LF(:,il)=lspf_all(:,index_LF+1); % LF oscillation
    lspf_HF(:,il)=lspf_all(:,index_HF+1); % HF oscillation

    % total term in time and frequency domain
    LSPf_LF(:,il)=out.LSPf(:,index_LF+1); % (sum: constant term + single contribute in LF)
    LSPf_HF(:,il)=out.LSPf(:,index_HF+1); % (sum: constant term + single contribute in HF)

    % integral of the LSP in the freq domain
    LSPf_integral_LF(il)=mean(LSPf_LF(bandLF(1):bandLF(2),il))*delta_LF;
    LSPf_integral_HF(il)=mean(LSPf_HF(bandHF(1):bandHF(2),il))*delta_HF;

    % TFs
    H(:,il)=out.H;
    H_LF(:,il)=out.Hk(:,index_LF); 
    H_HF(:,il)=out.Hk(:,index_HF);

    %% spectral decomposition
    A=[1 -Am];
    [Sk,Sm{il},freq,pot,potn,poli,f]=lsp_spect_dec(A,C,Su,fs,nfft);
    % Sk - spectral components
    % Sm - whole spectrum
    % freq - pole frequencies
    % pot, potn - pole powers
    % poli - poles

    %%%%%%%% total power
    Pow=sum(pot); 
    
    %%%%%%%% spectral components
    LFspect(:,il)=Sk(:,2);
    HFspect(:,il)=Sk(:,1);
    
end

%% plot 
h0=figure('Color','w','WindowState','maximized');
subplot(1,3,3); % spectral information storage - HF term
for il=1:l
    plot(f,LSPf_HF(:,il),'LineWidth',2,...
        'Color',colrange(il,:));
    hold on;
end
% ylabel('[nats]')
legend('$s_Y^{(HF)}(f)$','Interpreter','Latex','Location','northwest');
legend Box off
xlim([0 fs/2]);
xticks(0:0.1:fs/2);
% xlabel('f [Hz]');
ylim([-1 2.5])
yticks([ ]);
line([0 fs/2],[0 0],'LineStyle','--',...
        'Color',grey_line,...
        'LineWidth',axislinewidth,...
        'HandleVisibility','off')
hold on;
set(gca,'FontSize',DimFont,'FontName','Times', ...
    'Position',[0.7 0.3 0.2 0.6]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end
colormap(ax,colrange) 
cb1 = colorbar;
cb1.Position = [0.93 0.3 0.01 0.6];
cb1.Ticks=[0 0.5 1];
% cb1.TickLabels={horzcat('{\it b=0}'),horzcat('{\it b=0.5}'),horzcat('{\it b=1}')};

subplot(1,3,1); % spectral information storage - constant term
for il=1:l
    plot(b(il),LSP(il),'o',...
        'MarkerSize',5,...
        'MarkerEdgeColor',colrange(il,:),...
        'MarkerFaceColor',colrange(il,:));
    hold on;
end
legend('$S_Y(f)$','Interpreter','Latex','Location','southeast');
legend Box off
xticks([0 0.5 1])
% xticklabels({'b=0', 'b=0.5', 'b=1'});
ylim([0 0.6])
yticks([0 0.3 0.6]);
xlim([0 b(end)]);
hold on;
set(gca,'FontSize',DimFont,'FontName','Times', ...
    'Position',[0.1 0.3 0.2 0.6]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end

subplot(1,3,2); % spectral information storage - freq specific term
for il=1:l
    plot(f,LSPf_LF(:,il),'LineWidth',2,...
        'Color',colrange(il,:));
    hold on;
end
% ylabel('[nats]')
legend('$s_Y^{(LF)}(f)$','Interpreter','Latex','Location','northeast');
legend Box off
xlim([0 fs/2]);
xticks(0:0.1:fs/2);
% xlabel('f [Hz]');
ylim([-1 2.5])
yticks([-1 0 1 2]);
line([0 fs/2],[0 0],'LineStyle','--',...
        'Color',grey_line,...
        'LineWidth',axislinewidth,...
        'HandleVisibility','off')
hold on;
set(gca,'FontSize',DimFont,'FontName','Times', ...
    'Position',[0.4 0.3 0.2 0.6]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end

%% plot 
h1=figure('Color','w','WindowState','maximized');
subplot(2,3,[1 4]) % spectra for b==1
plot(f,Sm{1,l},'Color',[0.85 0.325 0.098],'LineWidth',2.5); % spettro
hold on;
plot(f,LFspect(:,l),'LineStyle','--',...
     'Color',[0.4660 0.6740 0.1880],...
     'LineWidth',2.2);
hold on;
plot(f,HFspect(:,l),'LineStyle','--',...
     'Color',[0.4940 0.1840 0.5560],...
     'LineWidth',2.2);
legend({'$P(f)$','$P^{(LF)}(f)$','$P^{(HF)}(f)$'},...
    'Interpreter','Latex','Location','northwest');
legend('boxoff')
ylim([0 11]);
xlim([0 fs/2]);
xticks(0:0.1:fs/2);
% xlabel('f [Hz]');
ylabel('PSD');
set(gca,'FontSize',DimFont,'FontName','Times',...
    'Position',[0.08 0.15 0.25 0.8]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end

subplot(2,3,[2 5]); % spectral information storage - total term
for il=1:l
    plot(f,H(:,il),'LineWidth',2,...
        'Color',colrange(il,:));
    hold on;
end
ylabel('TF')
legend('$H(f)$','Interpreter','Latex','Location','northwest');
legend Box off
xticks(0:0.1:fs/2);
xlim([0 fs/2]);
% xlabel('f [Hz]');
ylim([0 32])
yticks([0 10 20 30]);
line([0 fs/2],[0 0],'LineStyle','--',...
        'Color',grey_line,...
        'LineWidth',axislinewidth,...
        'HandleVisibility','off')
hold on;
set(gca,'FontSize',DimFont,'FontName','Times',...
    'Position',[0.4 0.15 0.25 0.8]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end

subplot(2,3,3); % spectral information storage - freq specific term
for il=1:l
    plot(f,H_HF(:,il),'LineWidth',2,...
        'Color',colrange(il,:));
    hold on;
end
% ylabel('[nats]')
legend('$H^{(HF)}(f)$','Interpreter','Latex','Location','northwest');
legend Box off
xlim([0 fs/2]);
xticks([]);
% xlabel('f [Hz]');
ylim([0 32])
yticks([0 10 20 30]);
line([0 fs/2],[0 0],'LineStyle','--',...
        'Color',grey_line,...
        'LineWidth',axislinewidth,...
        'HandleVisibility','off')
hold on;
set(gca,'FontSize',DimFont,'FontName','Times',...
    'Position',[0.7 0.58 0.25 0.37]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end

subplot(2,3,6); % spectral information storage - constant term
for il=1:l
    plot(f,H_LF(:,il),'LineWidth',2,...
        'Color',colrange(il,:));
    hold on;
end
% ylabel('[nats]')
legend('$H^{(LF)}(f)$','Interpreter','Latex','Location','northeast');
legend Box off
xticks(0:0.1:fs/2);
% xlabel('f [Hz]');
ylim([0 25])
yticks([0 10 20]);
xlim([0 fs/2]);
line([0 fs/2],[0 0],'LineStyle','--',...
        'Color',grey_line,...
        'LineWidth',axislinewidth,...
        'HandleVisibility','off')
hold on;
set(gca,'FontSize',DimFont,'FontName','Times',...
    'Position',[0.7 0.15 0.25 0.37]);
ax=gca;
ax.LineWidth=axislinewidth;
box off;
isholdonque = ishold;
ax1=axis;
hold on
plot(ax1(2)*[1,1],ax1(3:4),'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
hold on
plot(ax1(1:2),ax1(4)*[1,1],'k', ...
    'linewidth',axislinewidth,...
    'HandleVisibility','off')
if isholdonque == 0
hold off
end
