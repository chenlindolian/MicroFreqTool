% give the example of signal evoultion of HCFM with two fiber population
% Results are shown in Figure 1 of manuscript
% H0 = [1 0 0]' and H0 = [0 1 0]' for parallel fibers, two fibers crossing
% at 45° and 90°.
addpath('~/Desktop/WM/code/SimulateSignal/')
clc,clear all

theta = 90;
phi = 90;
Params.B0 = 7;              % Tesla
Params.gamma = 42.577e6;      % Hz/T
H0 = [sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)]';

% myelin (required: T2, xi, xa)
Params.myelin.Mag= 0.7; 
Params.myelin.T2 = 8*1e-3;
Params.myelin.T1 = 500*1e-3;
Params.myelin.proton_density = 1; 

Params.myelin.chii = -0.10;  % myelin anisotropic susceptibility (ppm)
Params.myelin.chia = -0.10;  % myelin isotropic susceptibility (ppm)

% intra axonal (required: T2)
Params.intra_axonal.Mag = 1;
Params.intra_axonal.T2 = 36*1e-3;
Params.intra_axonal.T1 = 1.5;
Params.intra_axonal.proton_density= 1; 
Params.intra_axonal.chii= 0; 

% extra axonal (required: T2)
Params.extra_axonal.Mag = 1;
Params.extra_axonal.T2 = 36*1e-3;
Params.extra_axonal.T1 = 1.5;
Params.extra_axonal.proton_density = 1; 
Params.extra_axonal.chii = 0; 
Params.E = 0.02;

load('~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus/317Fiber2BundleCrossAngle0.mat');
pad_size = [50 50 50];
total_chi = padarray(total_chi,pad_size,0,'both');
Phantom.model = padarray(Phantom.model,pad_size,1,'both');

deltaB0 = SimulateFieldFrom3DTensor(total_chi, Phantom, Params, H0);

deltaB0 = deltaB0(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
Phantom.model = Phantom.model(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
options.mask = Phantom.mask;
figure(11),createHistogramField(Phantom.model, deltaB0, options)
title('Parallel fiber along Z axis', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')

Params.TEs = (0:1:30)*1e-3;
signal0 = simulateSignalFrom3DField(Phantom, deltaB0, Params);
magn_signal0 = abs(signal0.total_normalized);
phase_signal0 = phase(signal0.total_normalized);
freq_signal0 = phase_signal0./(2*pi*Params.TEs);

figure(12)
imshow(squeeze(deltaB0(:,:,100)),[-30,30]),...
hold on
rectangle('position',[21,21,160,160],'edgecolor','b','LineWidth',1)
set(gca, 'YDir','Normal'),
% title('Zslice:100', 'FontWeight','bold','Fontsize',30,'fontname','Times New Roman'),...
h = colorbar;
set(get(h,'Title'),'string','Hz')
set(h,'Ticks',[-30,-20,-10,0,10,20,30]);
set(h,'TickLabels',{'-30','-20','-10','0','10','20','30'});
set(h,'LineWidth',2);
set(h,'FontSize',25);
set(h,'fontname','Times New Roman');
set(gca,'LineWidth',1,'FontSize',25);

load('~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus/317Fiber2BundleCrossAngle45.mat');
pad_size = [50 50 50];
total_chi = padarray(total_chi,pad_size,0,'both');
Phantom.model = padarray(Phantom.model,pad_size,1,'both');

deltaB45 = SimulateFieldFrom3DTensor(total_chi, Phantom, Params, H0);

deltaB45 = deltaB45(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
Phantom.model = Phantom.model(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
options.mask = Phantom.mask;
figure(21),createHistogramField(Phantom.model, deltaB45, options)
title('Two fibers crossing at 45\circ', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')

Params.TEs = (0:1:30)*1e-3;
signal45 = simulateSignalFrom3DField(Phantom, deltaB45, Params);
magn_signal45 = abs(signal45.total_normalized);
phase_signal45 = phase(signal45.total_normalized);
freq_signal45 = phase_signal45./(2*pi*Params.TEs);

figure(22)
imshow(squeeze(deltaB45(:,:,100)),[-30,30]),...
hold on
rectangle('position',[21,21,160,160],'edgecolor','b','LineWidth',1)
set(gca, 'YDir','Normal'),
% title('Zslice:100', 'FontWeight','bold','Fontsize',30,'fontname','Times New Roman'),...
h = colorbar;
set(get(h,'Title'),'string','Hz')
set(h,'Ticks',[-30,-20,-10,0,10,20,30]);
set(h,'TickLabels',{'-30','-20','-10','0','10','20','30'});
set(h,'LineWidth',1);
set(h,'FontSize',25);
set(h,'fontname','Times New Roman');
set(gca,'LineWidth',1,'FontSize',25);

load('~/Desktop/WM/code/Exampledata/HCFM/HCFM2FiberSus/317Fiber2BundleCrossAngle90.mat');
pad_size = [50 50 50];
total_chi = padarray(total_chi,pad_size,0,'both');
Phantom.model = padarray(Phantom.model,pad_size,1,'both');

deltaB90 = SimulateFieldFrom3DTensor(total_chi, Phantom, Params, H0);

deltaB90 = deltaB90(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
Phantom.model = Phantom.model(pad_size(1)+1:end-pad_size(1),pad_size(2)+1:end-pad_size(2),pad_size(3)+1:end-pad_size(3));
options.mask = Phantom.mask;
figure(31),createHistogramField(Phantom.model, deltaB90, options)
title('Two fibers crossing at 90\circ', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')

Params.TEs = (0:1:30)*1e-3;
signal90 = simulateSignalFrom3DField(Phantom, deltaB90, Params);
magn_signal90 = abs(signal90.total_normalized);
phase_signal90 = phase(signal90.total_normalized);
freq_signal90 = phase_signal90./(2*pi*Params.TEs);

figure(32)
imshow(squeeze(deltaB90(:,:,100)),[-30,30]),...
hold on
rectangle('position',[21,21,160,160],'edgecolor','b','LineWidth',1)
set(gca, 'YDir','Normal'),
% title('Zslice:100', 'FontWeight','bold','Fontsize',30,'fontname','Times New Roman'),...
h = colorbar;
set(get(h,'Title'),'string','Hz')
set(h,'Ticks',[-30,-20,-10,0,10,20,30]);
set(h,'TickLabels',{'-30','-20','-10','0','10','20','30'});
set(h,'LineWidth',2);
set(h,'FontSize',25);
set(h,'fontname','Times New Roman');
set(gca,'LineWidth',1,'FontSize',25);

 figure(41)
 plot(Params.TEs*1e3, magn_signal0, 'r', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, magn_signal45, 'b', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, magn_signal90, 'k', 'LineWidth',1)
 xlabel('Echo time (ms)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 ylabel('Magnitude (a.u.)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 title('Signal magnitude', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')
 axis([0 30 0 1]);
 set(gca,'XTick',[0:10:30],'LineWidth',2,'YTick',[0:0.2:1],'LineWidth',1,'FontSize',25);
 leg = legend('Parallel fibers', 'Two fibers crossing at 45\circ', 'Two fibers crossing at 90\circ');
 set(leg,'FontSize',20);
 set(leg,'fontname','Times New Roman');

 figure(42)
 plot(Params.TEs*1e3, freq_signal0, 'r', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, freq_signal45, 'b', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, freq_signal90, 'k', 'LineWidth', 1)
 xlabel('Echo time (ms)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 ylabel('{\itf} (Hz)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 %title('Signal frequency', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')
 axis([1 30 -8 0]);
 set(gca,'XTick',[0:10:30],'LineWidth',1,'YTick',[-8:2:0],'LineWidth',1,'FontSize',25,'fontname','Times New Roman');
leg = legend('Parallel fibers', 'Two fibers crossing at 45\circ', 'Two fibers crossing at 90\circ');
 set(leg,'FontSize',20);
 set(leg,'fontname','Times New Roman');
% set(leg,'Location','southwest');

 figure(43)
 plot(Params.TEs*1e3, phase_signal0, 'r', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, phase_signal45, 'b', 'LineWidth', 1)
 hold on
 plot(Params.TEs*1e3, phase_signal90, 'k', 'LineWidth', 1)
 xlabel('Echo time (ms)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 ylabel('Phase (rad)','FontSize',25,'FontWeight','Bold','fontname','Times New Roman')
 title('Signal phase', 'FontWeight', 'bold','Fontsize',30,'fontname','Times New Roman')
 axis([0 30 -1.4 0]);
 set(gca,'XTick',[0:10:30],'LineWidth',2,'YTick',[-1.4:0.2:0],'LineWidth',1,'FontSize',25,'fontname','Times New Roman');
 leg = legend('Parallel fibers', 'Two fibers crossing at 45\circ', 'Two fibers crossing at 90\circ');
 set(leg,'FontSize',20);
 set(leg,'fontname','Times New Roman');
