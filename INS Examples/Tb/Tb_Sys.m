% Simulation reproducing Figure 4b of % Chem. Commun., 2016,52, 2091-2094 clear all
% some constants
function Sys1 = Tb_Sys(Tb_Type)
close all
rcm = 29979.2458;   %rcm
meV = rcm*8.065;    %meV
clear Sys1
Sys1.S = 6;
Exp.SpectrumType = 'SE';
Exp.Temperature = 30;
Exp.Energy = 0:0.05:11; Exp.lwfwhm = 0.87; Exp.Q = 1.15:0.05:1.65;
Sys1.FormFactor = 'Tb3';
Sys1.Coords = [0 0 0];
%Tb_Type = 'a';
 %Tb_Type = 'b';
if Tb_Type =='a'
%    Sys1.B2 = [0 0 0.705 0 0.0250].*meV; 
    Sys1.B2 = [0. 0 0.705 0 00250].*meV; 
    Sys1.B4 = [-8.5E-4 0 0 0 0 0 0 0 0].*meV;
%     [cross,Sys1.ev] = mint(Sys1,Exp);
    EE = eig(sham(Sys1,[0,0,0]));
    Sys1.EE = EE-EE(1);
elseif Tb_Type =='b'
%    Sys1.B2 = [0 0 0.712 0 0.0729].*meV; 
    Sys1.B2 = [0.0729 0 0.712 0 0].*meV; 
    Sys1.B4 = [0 0 0 0 0 0 0 0 0].*meV;
%      [cross,Sys1.ev] = mint(Sys1,Exp);
         EE = eig(sham(Sys1,[0,0,0]));
    Sys1.EE = EE-EE(1);
end
 %[cross,Sys1.ev] = mint(Sys1,Exp); figure
 %plot(Exp.Energy,cross.*3); hold on

%%
% % run simulation and plot over previous simulation in a 3 to 7 scaled ratio to match figure 4b of Chem. Commun., 2016,52, 2091-2094 
% clear Sys
% close all
%  [A,dint,ev,Ops] = SysInput(Sys1);
%  Y=General_IEP(A,ev,'dint',round(dint,1,'significant'),'NMinima',2,'method','N','linesearch','basic');
% Sys_Out = Sys_Output(Y,Ops,Sys1.S);
% Sys = Sys_Out(1);
% Sys.FormFactor = 'Tb3';
% Sys.Coords = [0 0 0];
% 
% 
% % Tb3+ J=6
% [cross,Eigs] = mint(Sys,Exp);
% plot(Exp.Energy,cross.*7);
% % format figure axes
% set(gca,'xtick',[0:1:11]); xlabel('Energy (meV)'); ylabel('Intensity (arb. units)');