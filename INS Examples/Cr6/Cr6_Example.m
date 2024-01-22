%Cr6 - Example
clear Sys Exp
%Simulate experimental data
rcm = 29979.2458;   
meV = rcm*8.065;  
N_electrons = 6;    %Number of Cr atoms in chain
N_eigenvalues = 21; %Number of Eigenvalues known experimentaly
[Exp] = Cr_solution(N_electrons,N_eigenvalues); %Calculates simulated eigenvalues


%Set up model to be optimised
r=1;    %Number of significant figures initial guiess uses.
%Set up initial guesses for parameters
B20 = round((-0.041/3)*meV,r,'significant'); %B20 = 1000;
B22 = round(0.007*meV,r,'significant');    %B22=-1000;
Jval = round(1.46*meV,r,'significant');     %Jval = 10000
S=3/2;  Sys.S = [S];    %set up Spins
B2 = [B22 0 B20 0 0];   Sys.B2 = [B2];  %set up Stevens parameter
Sys.J = [];    %set up exchange terms
for i = 2:N_electrons   %Loop over all electrons
    Sys.S = [Sys.S,S];
    Sys.B2 = [Sys.B2;B2];
    %     Sys.J = [Jval,zeros(1,i-2),Sys.J];  %same value of J for all nearest neighbours
    Sys.J = [Jval+i,zeros(1,i-2),Sys.J];   %different value of J for nearest neighbours
end
% Sys.J = Sys.J+[1 0 0 0 0 2 0 0 0 3 0 0 2 0 1];
Vary = Sys; %This will vary all non-zero parameters


%Optimse over parameters given
Opt = struct('NMinima',1,'Method','Newton','Linesearch','Basic',...
    'MaxIter',1000,'theta',2,'StepTolerance',1e-6,'GradientTolerance',1e-1,...
    'ObjectiveTolerance',1e-1,'Minalpha',1e-10,'Scaled',1,...
    'deflatelinesearch',1,'IEPType','Difference','Verbose',1,'tau',0.5);

[SysOut, NIter, Flags, Iters, FinalError]= INS_IEP(Sys,Vary,Exp,Opt);

% SysOut.B2

%%
%Set up Mint variables - only works for Nelectrons = 6.
% Mint still runs with less than 6 electrons? Just ignores unncessary
% coordinates?
a = [0;-3.3807;-5.7290;-5.5960;-3.1148;0.2482];
b = [0;0.0840;-2.3170;-5.6470;-7.8270;-7.7410];
c = [0; -0.2765 ; 0.1103 ; -0.2956 ; 0.2920 ; -0.1803];
SysMint.Coords = [a b c];
SysMint.FormFactor = 'Cr3';
MintOpt.NumEigs = 24;
ExpMint.SpectrumType = 'SE';
ExpMint.Temperature = [1.5, 6, 15];
ExpMint.Energy = -1:0.001:4;
ExpMint.lwfwhm = 0.2;
ExpMint.Q = [0.1,0.5,1,1.5,2]; 
%have to 'press any key' in the command window to move between minima
for i = 1:length(SysOut)
    SysMint = SysOut(i);
    SysMint.Coords = [a b c]; SysMint.FormFactor = 'Cr3';
    [cross,Eigs] = mint(SysMint,ExpMint,MintOpt);
    plot(ExpMint.Energy,cross)
    xlim([0 3.25])
    xlabel('Energy (meV)')
    ylabel('Intensity (arb. units)')
    legend('1.5 K','6 K','15 K')
    if i <length(SysOut)
        pause
    end
end



%function to simulate the experimental eigenvalues
function [Exp,Sys] = Cr_solution(N_electrons,Neigenvalues)
rcm = 29979.2458;   meV = rcm*8.065;    %conversions
B20 = (-0.041/3).*meV; B22 = 0.007.*meV; % Known stevens parameters
S=3/2;  Sys.S = [S];    %set up Spins
B2 = [B22 0 B20 0 0];   Sys.B2 = [B2];  %set up Stevens parameters
Jval =1.46*meV; Sys.J = [];    %set up exchange terms
for i = 2:N_electrons   %Loop over all electrons
    Sys.S = [Sys.S,S];
    Sys.B2 = [Sys.B2;B2];
    Sys.J = [Jval,zeros(1,i-2),Sys.J];
end
%Calculate simulated eigenvalues
[~,EE] = eigs(ham(Sys,[0,0,0],'sparse'),Neigenvalues+10,'smallestreal');
EE=diag(EE);
Exp.ev=EE(1:Neigenvalues) - EE(1);
end
