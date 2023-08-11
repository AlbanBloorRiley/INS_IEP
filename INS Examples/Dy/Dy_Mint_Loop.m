%% Simulation and Comparison  6A energy
%Only need to run this once!

%    folder = '/Volumes/mlbakerlab/neutron_data/FRMII/6396_TOFTOF/analysis/DyHNN2_analysis/2021/data/feb2022/for_Alan/';
   folder = '/Users/user/Dropbox (The University of Manchester)/alban.bloorriley@manchester.ac.uk’s files/MatLab/mint/for_Alan/';
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
%% Run whatever code you want to set up your Sys structure 
% - inputing a nonzero value for all the parameters you wish to use.
[Sys,Vary,Exp]=Dy_Sys;

Opt = struct('NMinima',5,'Method','Newton','Linesearch','Basic',...
    'MaxIter',1000,'p',2,'StepTolerance',1e-6,'ObjectiveTolerance',1e-1,'GradientTolerance',1e-1,'Minalpha',1e-18,'Scaled',1,...
    'deflatelinesearch',0,'IEPType','Difference','Verbose',0,'tau',0.5);

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
%     if ~(SysOut(k).Flags  == "error <tol")
%         continue
%     end
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

%% New easyspin
Exp.Field = 100;
Exp.Temperature = Dy_Susceptibility(:,1);
[~,chizz] = curry(SysOut(4),Exp);

%one of these
% plot(Exp.Temperature, chizz)
%plot(Exp.Temperature, Exp.Temperature'.*chizz./(4*pi*1e-6))

%plot Dy_Susceptibility for the exp data

