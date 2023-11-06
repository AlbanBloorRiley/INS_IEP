%% Simulation and Comparison  6A energy
%Only need to run this once!

%    folder = '/Volumes/mlbakerlab/neutron_data/FRMII/6396_TOFTOF/analysis/DyHNN2_analysis/2021/data/feb2022/for_Alan/';
   folder = '/Users/user/Dropbox (The University of Manchester)/alban.bloorriley@manchester.ac.uk’s files/MatLab/mint/Dy_data/';
   filename = 'DyHNN2_IofE_6A_1K_2theta_gt_30_bkgdsub=0p1.inx';
   Dy_6A1K = importdata([folder,filename]);
   filename = 'DyHNN2_IofE_6A_9K_2theta_gt_30_bkgdsub_0p1.inx';
   Dy_6A9K = importdata([folder,filename]);
   A= {Dy_6A1K,Dy_6A9K};
   %for i=1:length(A)
   n=2;
       x_avg_1 = arrayfun(@(j) mean(A{1}(j:j+n-1,1)),1:n:length(A{1}(:,1))-n+1)';
    y_avg_1 = arrayfun(@(j) mean(A{1}(j:j+n-1,2)),1:n:length(A{1}(:,1))-n+1)';
    err_avg_1 = arrayfun(@(j) sqrt(mean(A{1}(j:j+n-1,3).^2)),1:n:length(A{1}(:,1))-n+1)'; 
    x_avg_2 = arrayfun(@(j) mean(A{2}(j:j+n-1,1)),1:n:length(A{2}(:,1))-n+1)';
    y_avg_2 = arrayfun(@(j) mean(A{2}(j:j+n-1,2)),1:n:length(A{2}(:,1))-n+1)';
    err_avg_2 = arrayfun(@(j) sqrt(mean(A{2}(j:j+n-1,3).^2)),1:n:length(A{2}(:,1))-n+1)'; 
   %end

   
   Dy_coords = Dy_Coords;
   rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
%Tesla = meV*0.116;   % Tesla         to MHz


MintOpt.NumEigs = 10;
MintExp.SpectrumType = 'SE'; 
MintExp.Temperature = [1, 3.5, 9]; 
MintExp.Energy = -1.5:0.001:1.5; 
MintExp.lwfwhm = 0.05;
MintExp.Q = [0.1,0.5,1,1.5,2]; 
%%
clear Sys
Sys.S = [1/2 15/2 1/2];
B20 = -0.173585036568371*meV;%-0.1736;%-0.197;
B22 = 0.443886669455654*meV;%0.4439;%0.066;
B40 = -0.000767621238475*meV;%3.137e-4;%1.5e-4;
B60 = -0.000125378265669*meV;%-1.8351e-6;%-1.25e-5;

Jex1 = 0.076588615657158*meV;%0.9; %Unit in meV %-0.54;
Jex2 = 0.048947841762338*meV;%0.6; %Unit in meV %-0.298;
% Sys.dint = [B20;B22;B40;B60;Jex1;Jex2];

Sys.B2 = [0 0 0 0 0 ; B22 0 B20 0 0 ; 0 0 0 0 0 ];
% Sys.B2 = [0 0 0 0 0 ; 0 0 B20 0 0 ; 0 0 0 0 0 ];
Sys.B4 = [0 0 0 0 0 0 0 0 0; 0 0 0 0 B40 0 0 0 0; 0 0 0 0 0 0 0 0 0];
Sys.B6 = [0 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 B60 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 0];
Sys.J = [Jex1 0 Jex2];
Vary.J = [1 0 1];
% Sys.ee = [diag([Jex1,Jex1,Jex1]);zeros(3);diag([Jex2,Jex2,Jex2])];
% 
% Sys.ee = [diag([Jex1,Jex1,Jex1+1]);zeros(3);diag([Jex2,Jex2,Jex2+1])];
% Vary.ee = [[1,0,1;0,1,0;1,0,1];zeros(3);[1,0,1;0,1,0;1,0,1]];

Vary.B2 = [zeros(1,5) ; 1 0 1 0 0 ; zeros(1,5)];
Vary.B4 = [zeros(1,9); 0 0 0 0 1 0 0 0 0; zeros(1,9)];
Vary.B6 = [zeros(1,13); 0 0 0 0 0 0 1 0 0 0 0 0 0; zeros(1,13)];

%dint = round(dint,2);
EE = [0;0;0.72;0.72;1.12;1.12;1.28;1.28].*meV;
Exp.ev=EE;
%% Run whatever code you want to set up your Sys structure 
% - inputing a nonzero value for all the parameters you wish to use.
% [Sys,Vary,Exp]=Dy_Sys;

Opt = struct('NMinima',2,'Method','Gauss-NewtonT2','Linesearch','Basic',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-2,'ObjectiveTolerance',1e-2,...
    'GradientTolerance',1e-1,'Minalpha',1e-18,'Scaled',0,...
    'deflatelinesearch',1,'IEPType','Classic','Verbose',0,'tau',0.8);

[SysOut, NIter, Flags, Iters, FinalError] = INS_IEP(Sys,Vary,Exp,Opt)

% adding (...,'SysFound',SysFound) to the end will start the loop as if
% SysFound has already been found.
%% Run this to simulate the INS Spectrum
startminima= 1; %change this to start simulating at a differnt minima
for k = startminima:size(SysOut,2)
%     SysOut(k).B2 = Sysfixed.B2;
%     SysOut(k).B4 = Sysfixed.B4;
%     SysOut(k).B6 = Sysfixed.B6;
%         if ~(SysOut.Flags{k} == "error<tol")&&~(Flags{k} == "Input Minima")
    if (Flags{k}  == "Max Iterations reached")
        continue
    end
%     if SysOut(k).B2(2,3)>0
%         continue
%     end
    
    Dy_Mint;
    SysOut(k)
   
    k

    %50?
    %plot Experimental data
    hold on
    errorbar(x_avg_1, y_avg_1, err_avg_1,'o-','LineWidth',1.5);
    errorbar(x_avg_2, y_avg_2, err_avg_2,'o-','LineWidth',1.5);
    %
    xlim([0,1.5]); %1K
    ylim([0,80]);% 1K
    legend('Sim 1 K','Sim 3.5 K','Sim 9 K','Exp 1K','Exp 9K')
    hold off
    pause

end
%% Good fits?
goodfits = [1 13 24 36 37  52 54 65 97 99 154 155 156 164 169 192 200];
for k = goodfits
%     k = goodfits()
    MintSys = SysFound_NBC(k);
    MintSys.FormFactor = {'Cu2','Dy3','Cu2'}; %elements included i.e. Dy(3+) and neutral oxygen radical
    MintSys.Coords = Dy_Coords;
    [cross,Eigs] = mint(MintSys,Exp);
    plot(Exp.Energy,cross);
    xlim([0;Exp.Energy(end)])
    xlabel('Energy (meV)')
    xticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5])
    ylabel('Intensity (arb. units)')
    legend('1 K','3.5 K','9 K')

    %Plot experimental data
    hold on
    errorbar(x_avg_1, y_avg_1, err_avg_1,'o-','LineWidth',1.5);
    errorbar(x_avg_2, y_avg_2, err_avg_2,'o-','LineWidth',1.5);
    ylim([0,80]);% 1K
    legend('Sim 1 K','Sim 3.5 K','Sim 9 K','Exp 1K','Exp 9K')
    hold off
    SysOut(k)
    k
    pause
end




%% New easyspin
Exp.Field = 100;
Exp.Temperature = Dy_Susceptibility(:,1);
[~,chizz] = curry(SysOut(1),Exp);

%one of these
plot(Exp.Temperature, chizz)
figure
plot(Exp.Temperature, Exp.Temperature'.*chizz./(4*pi*1e-6))

%plot Dy_Susceptibility for the exp data

