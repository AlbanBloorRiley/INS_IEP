clear Sys

% conversions:
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz

% assume two magnetic centres have isotropic interaction
B20 = -0.173585036568371;%-0.1736;%-0.197;
B22 = 0.443886669455654;%0.4439;%0.066;
B40 = -0.000767621238475;%3.137e-4;%1.5e-4;
B44 = 0;%3.6e-3;%5e-4;
B60 = -0.000125378265669;%-1.8351e-6;%-1.25e-5;
B64 = 0;%-1.7855e-5;%1.17e-4;
B66 = 0;%3.3974e-5;


% exchange coupling in -2JS1S2 formalism
Jex1 = 0.076588615657158;%0.9; %Unit in meV %-0.54;
Jex2 = 0.048947841762338;%0.6; %Unit in meV %-0.298;
Jex1 = 1;%0.9; %Unit in meV %-0.54;
Jex2 = 1;%0.6; %Unit in meV %-0.298;
Jex3=0;
%Jex4 = 1.7;
%Jex5 = 0.5;


%Setup of J=15/2 for Dy3+ and S=1/2 for each of the two HNN radical ligands
Sys.S = [15/2 1/2 1/2];

% Exchange coupling
%Between Dy (1) and radical (2)
Jxx12 = Jex1;
Jyy12 = Jex1;
Jzz12 = Jex1;
%Between Dy(1) and radical (3)
Jxx13 = Jex2;
Jyy13 = Jex2;
Jzz13 = Jex2;
%Between radical (2) and radical (3)
Jxx23 = Jex3;
Jyy23 = Jex3;
Jzz23 = Jex3;

% Build up Hamiltonian
Sys.B2 = [B22 0 B20 0 0 ; zeros(1,5); zeros(1,5)].*meV;        % B(k=2,q) with q = +2,+1,0,-1,-2
Sys.B4 = [B44 0 0 0 B40 0 0 0 0 ; zeros(1,9); zeros(1,9)].*meV;  % B(k=4,q) with q = +4,+3,+2,+1,0,-1,-2,-3,-4
Sys.B6 = [B66 0 B64 0 0 0 B60 0 0 0 0 0 0 ; zeros(1,13); zeros(1,13)].*meV;  % B(k=6,q) with q = +6,+4,+3,+2,+1,0,-1,-2,-3,-4,-6


Sys.ee = -2*[Jxx12 Jyy12 Jzz12; Jxx13 Jyy13 Jzz13; Jxx23 Jyy23 Jzz23].*meV;

% gDy = 4/3;
% g = 2;
% Sys.g = [gDy g g]; %mean of g1 for radical 1 and g2 for radical 2

H = sham(Sys,[0,0,0]); %magnetic field (mT) [x,y,z]
[Vecs,E]=eig(H);
EE = diag(E);
Eigs = (EE-EE(1))./meV; %convert to meV with the lowest eig set to zero.
figure;
plot_eigenvectors(Vecs(1:12,1:16));
% str = sprintf('B20 = %1.3f J1= %1.2f J2=%1.2f',B20,Jex1,Jex2);
% title(str,'FontSize',16);


d=1;
mat = zeros(64,3); %Build a matrix to record eigenvectors
for i=-15/2:15/2
    for j=-1/2:1/2
        for z=-1/2:1/2
            mat(d,:)=[i,j,z];
            d = d+1;
        end
    end
end


%% Calculate magnetic properties

% calculate magnetisation versus magnetic field
Exp.Temperature = 2; % Kelvin
Exp.Field = [0.1:100:10000]; %mT
Opt.Output = 'MvsBCGS'; % specify output format
mu = curry(Sys,Exp,Opt);

load('Dy_Magnetisation.mat');
% make figure
figure
plot(Dy_Magnetisation(:,1),Dy_Magnetisation(:,2));
hold on;
plot(Exp.Field./1000,mu);
xlabel('H (Tesla)'); ylabel('M (\mu_B)'); legend('Experimental','Simulation');
ylim([0,10]);
set(gca  ,'FontSize', 16);
hold off
% calculate magnetic susceptibility versus temperature in 0.T magnetic
% field
Exp.Temperature = [1:300]; % Kelvin
Exp.Field = [100]; % 100 mT = 0.1 T
Opt.Output = 'ChiTCGS'; % specify output format
chiT = curry(Sys,Exp,Opt);
load('Dy_Susceptibility.mat');

%make figure
figure
plot(Dy_Susceptibility(:,1),Dy_Susceptibility(:,2));
hold on
plot(Exp.Temperature,chiT);
xlabel('Temperature (K)'); ylabel('$\chi$T (cm$^3$ K mol$^{-1}$)','interpreter','latex');
%ylim([0,16]);
legend('Experimental','Simulation');
set(gca  ,'FontSize'   , 16);
hold off

%% Calculation of INS spectra

% SE calculation
Exp.Energy = [0:0.01:5]; % spectrum energy range in meV
Exp.lwfwhm = .1; % gaussin broadening FWHM in meV
Exp.Q = [0.1,0.5,1.0,1.5,2.0]; % Q range in A^-1
Exp.SpectrumType = 'SE' ; % S(E) only
Exp.Temperature = [2,9]; % Kelvin
Sys.FormFactor = {'Dy3','Cu2','Cu2'}; %elements included i.e. Dy(3+) and neutral oxygen radical
Sys.Coords = [0 0 0 ; 0 0 2.362; 0 0 -2.326]; % distance between Dy and each of HNN radical

[cross_sect,Eigs]=mint(Sys,Exp);




%make figure
folder = '/Volumes/mlbakerlab/neutron_data/FRMII/6396_TOFTOF/analysis/DyHNN2_analysis/2021/data/feb2022/for_Alan/';
filename = 'DyHNN2_IofE_6A_1K_2theta_gt_30_bkgdsub=0p1.inx';
Dy_6A1K = importdata([folder,filename]);
filename = 'DyHNN2_IofE_6A_9K_2theta_gt_30_bkgdsub_0p1.inx';
Dy_6A9K = importdata([folder,filename]);
A= {Dy_6A1K,Dy_6A9K};

Multiple_Plotter(A); %This function contains a code for "figure"
hold on;
plot(Exp.Energy,cross_sect);
xlim([0,2.5]); %1K
ylim([0,80]);% 1K
legend('Exp 2K', 'Exp 9K ','Sim 2K','Sim 9K' ,'Fontsize',16);
set(gca,'FontSize',16);


% SQE map calculation
Exp.Energy = [0:0.05:15]; % spectrum energy range in meV
Exp.Q = [0.1:0.05:3]; %Q range in A^-1
Exp.SpectrumType = 'SQE' ; % SQE map calculation
[cross_sect,Eigs]=mint(Sys,Exp);

%make figure
figure
surf(Exp.Q,Exp.Energy,cross_sect)
view(0,90); shading interp; colorbar; axis tight
caxis([0,1]);
xlabel('Q ($\AA^{-1}$)','interpreter','latex')
ylabel('Energy (meV)')
ylim([0,1.5]);
