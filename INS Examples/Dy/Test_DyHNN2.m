clear all;
clc;

Bkq=importdata('DyO8H8_Bkq.txt');
k=num2str(Bkq(:,1));
q=num2str(Bkq(:,2));
e=Bkq(:,4)';

STR=cell(length(k),1);

for i=1:length(k)
STR{i}= strcat('B',k(i,:),q(i,:)); %String with space
STR{i}=STR{i}(~isspace(STR{i})); %String without space
end


% conversions:
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz


%Setup of J=15/2 for Dy3+ and S=1/2 for each of the two HNN radical ligands
Sys.S =[15/2;1/2;1/2];

Sys.B2 = [flip(e(1:5)); zeros(1,5);zeros(1,5)].*rcm; % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [flip(e(6:14)); zeros(1,9);zeros(1,9)].*rcm; % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
Sys.B6 = [flip(e(15:27));zeros(1,13);zeros(1,13)].*rcm; %[0 0 B64 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13); zeros(1,13)].*meV;  % B(k=6,q) with q = +6,+4,+3,+2,+1,0,-1,-2,-3,-4,-6
Sys.B8 = [flip(e(28:44));zeros(1,17);zeros(1,17)].*rcm;
Sys.B10 = [flip(e(45:65));zeros(1,21);zeros(1,21)].*rcm;
Sys.B12 = [flip(e(66:90));zeros(1,25);zeros(1,25)].*rcm;


gDy = 1.53; % gDy=4/3
ge=2.0023;
Sys.g = [gDy;ge;ge]; 

%Electronic Interaction

Jex1 = -1.06;% Dy-R(1) z-component %-1.76;%-0.878;%-1.76;%-0.869;%-0.4; %Unit in kelvin %The difference between J1 and J2 defines the separation between two cold peaks. 
Jex2 = -0.95;% Dy-R(2) z-component %-1.13;%-1.13;%-1.15;%-0.3; %Unit in kelvin 
Jex3 = -6.4; % R-R z-component
Jex4 = 0;    % Dy-R(1) xy-component
Jex5 = 0;    % Dy-R(2) xy-component
Jex6 = -4.5; % R-R  xy-component

%By including xx and yy components, no effects are observed.

% Exchange coupling
%Between Dy (1) and radical (2)
% Jxx12 = Jex4;
% Jyy12 = Jex4;
% Jzz12 = Jex1;
% %Between Dy(1) and radical (3)
% Jxx13 = Jex5;
% Jyy13 = Jex5;
% Jzz13 = Jex2;
% %Between radical (2) and radical (3)
% Jxx23 = Jex6;
% Jyy23 = Jex6;
% Jzz23 = Jex3;
% 
% Sys.ee =[Jxx12 Jyy12 Jzz12; Jxx13 Jyy13 Jzz13; Jxx23 Jyy23 Jzz23].*kelvin;
% 
N_electrons = 3;
%A0 = zfield(Sys,1:3,'sparse')
Sys.J = zeros(1,(N_electrons-1)*(N_electrons)/2);
H = sham(Sys,[0,0,0]); %magnetic field (mT) [x,y,z]
[Vecs,E]=eig(H);
EE = diag(E);
Eigs = (EE-EE(1))./meV; %convert to meV with the lowest eig set to zero.
A0 = H;

%% calculate magnetisation versus magnetic field
Exp.Temperature = 2; % Kelvin
Exp.Field = [0.1:100:10000]; %mT
Opt.Output = 'MvsBCGS'; % specify output format
mu = curry(Sys,Exp);

load('Dy_Magnetisation.mat');
% make figure
figure
plot(Dy_Magnetisation(:,1),Dy_Magnetisation(:,2),'LineWidth',2);
hold on;
plot(Exp.Field./1000,mu,'LineWidth',2);
xlabel('H (Tesla)'); ylabel('M (\mu_B)'); legend('Experimental','Simulation');
ylim([0,10]);
set(gca  ,'FontSize', 16);
hold off

% calculate magnetic susceptibility versus temperature in 0.T magnetic
% field
Exp.Temperature = 1:300; % Kelvin
Exp.Field = [100]; % 100 mT = 0.1 T
Opt.Output = 'ChiCGS'; % specify output format
chi = curry(Sys,Exp,Opt);
load('Dy_Susceptibility.mat');
%make figure
figure
plot(Dy_Susceptibility(:,1),Dy_Susceptibility(:,2),'LineWidth',2);
hold on
plot(Exp.Temperature,chi.*Exp.Temperature,'LineWidth',2);
xlabel('Temperature (K)'); ylabel('$\chi$T (cm$^3$ K mol$^{-1}$)','interpreter','latex');
%ylim([0,16]);
legend('Experimental','Simulation');
set(gca  ,'FontSize'   , 16);
hold off

%% INS Simulation

%6A Exp data
folder = '/Volumes/mlbakerlab/neutron_data/FRMII/6396_TOFTOF/analysis/DyHNN2_analysis/2021/data/feb2022/for_Alan/';
filename = 'DyHNN2_IofE_6A_1K_2theta_gt_30_bkgdsub=0p1.inx';
Dy_6A1K = importdata([folder,filename]);
filename = 'DyHNN2_IofE_6A_9K_2theta_gt_30_bkgdsub_0p1.inx';
Dy_6A9K = importdata([folder,filename]);
A= {Dy_6A1K,Dy_6A9K};

% 8A Exp data
filename = 'DyHNN2_IofE_8A_3p5K_2theta_gt_30_bkgdsub_0p1.inx';
Dy_8A3p5K = importdata([folder,filename]);
Dy_8A9K = load('IofE_Dy8A9K.mat');
B= {Dy_8A3p5K,Dy_8A9K.IofE_Dy8A9K};


% SE calculation
Exp.Energy = [0:0.01:5]; % spectrum energy range in meV
Exp.lwfwhm = 0.1; % gaussin broadening FWHM in meV
Exp.Q = [0.1,0.5,1.0,1.5,2.0]; % Q range in A^-1
Exp.SpectrumType = 'SE' ; % S(E) only
Exp.Temperature = [1,3.5,9]; % Kelvin
Sys.FormFactor = {'Dy3';'Cu2';'Cu2'}; %elements included i.e. Dy(3+) and neutral oxygen radical
Sys.Coords = [0 0 0; 0.84070 -1.85950 -1.18881; 0.76832   1.93275   1.04084]; % distance between Dy and each of HNN radical

[cross_sect,Eigs]=mint(Sys,Exp); %Output

Multiple_Plotter(A); %This function contains a code for "figure"
hold on;
colour={'b','m','r'};
for ii= 1:size(cross_sect,1)

plot(Exp.Energy,cross_sect(ii,:),'color',colour{ii},'LineWidth',2);
end
xlim([0,1.6]); %1K
ylim([0,80]);% 1K
legend('Exp 2K', 'Exp 9K ','Sim 1K','Sim 3.5K','Sim 9K' ,'Fontsize',16);
set(gca,'FontSize',16);
hold off


Multiple_Plotter(B);
hold on
for ii= 1:size(cross_sect,1)

plot(Exp.Energy,cross_sect(ii,:).*40,'color',colour{ii},'LineWidth',2);
end
xlim([0,0.8]); %1K
ylim([0,800]);% 1K
legend('Exp 3.5K', 'Exp 9K ','Sim 1K','Sim 3.5K','Sim 9K' ,'Fontsize',16);
set(gca,'FontSize',16);
hold off

%% Eigenvectors

figure;
plot_eigenvectors(abs(Vecs(1:20,:)),Eigs);
title('DyHNN_{2}','FontSize',16);
