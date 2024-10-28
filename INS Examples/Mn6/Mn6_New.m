
%Define unit conversions
rcm = 29979.2458;    % reciprocal cm to MHz
kelvin = rcm*0.695;  % kelvin        to MHz
meV = rcm*8.065;     % meV           to MHz
Tesla = meV*0.116;   % Tesla         to MHz
% The measured eigenvalues are
% exp_eigs = [0.14153, 0.59095, 0.59095, 0.9428, 0.9428, 0.9428, 1.2285+0.59095,1.2285+0.59095, 1.5128+0.9428, 2.603, 2.86];
% exp_errs = [0.00016, 0.00011, 0.0003, sqrt(0.0016^2+0.00011^2), sqrt(0.0016^2+0.0003^2), 0.012, 0.02];
% exp_eigs = [0.14153, 0.59095, 0.9428, 1.2285+0.59095, 1.5128+0.9428, 2.86];
% 0 meV, 1-fold
% 0.1414(3) meV, 1-fold
% 0.59070(16) meV, 2-fold
% 1.0841(5) meV, 3-fold (experimentally). In my simulation it is actully one singlet and one doublet which are pseudo-degenerate.
% 1.4134(3) meV, 2-fold
% 2.316(3) meV, 3- or 5-fold
% 2.5218(16) meV, 2-fold if the above is 5-fold. 4-fold is the above is 3-fold
% 2.603(12) meV, 1-fold
% 2.86(2) meV, 2-fold
  exp_eigs = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.5218, 2.5218, 2.5218, 2.5218]';   %or
% Difference: 1.4276e+05; Classic: 4.662164e+05
  % exp_eigs = [ 0, 0.1414, 0.59070, 0.59070, 1.0841, 1.0841 , 1.0841, 1.4134, 1.4134, 2.316, 2.316, 2.316, 2.316, 2.316,  2.5218, 2.5218]';


% Here are some somewhat reasonable starting parameters
clear Sys Exp Vary
Sys.S = [2 2 5/2 5/2 5/2 5/2]; % MnIII has S = 2, MnII has S = 5/2

J_S4_S4   = -7*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
J_S4_S5_1 = 0.41*meV;  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
J_S4_S5_2 = -0.4*meV; % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
J_S5_S5   = -0.1*meV;  % MnII  - MnII.                        Rodoplhe: Very weak and AFM. Free value



J_AB = J_S4_S4;
J_A1 = J_S4_S5_1; J_A3 = J_S4_S5_1; J_B2 = J_S4_S5_1; J_B4 = J_S4_S5_1; %For these, the MnIII JT axis is involved. Can be FM or AFM
J_A2 = J_S4_S5_2; J_A4 = J_S4_S5_2; J_B1 = J_S4_S5_2; J_B3 = J_S4_S5_2; %For these, the MnIII JT axis is NOT involved. Can only be AFM
J_12 = J_S5_S5;   J_34 = J_S5_S5; 
J_14 = -0.00.*meV;         J_23 = J_14; %Assumed zero. Only interacts via a pivalate bridge.
J_13 = -0.00.*meV;         J_24 = J_13; %Assumed zero. No interaction pathway

% %Trial
% J_S4_S4   = -5*meV;    % MnIII - MnIII.                       Rodolphe: Strong and AFM.    Keep fixed.
% J_S4_S5_1 = 0.1*meV;  % MnIII - MnII, MnIII JT involved.     Rodolphe: Weak, FM or AFM.   Free value
% J_S4_S5_2 = -0.3*meV; % MnIII - MnII, MnIII JT not involved. Rodolphe: Weak and AFM.      Free value
% J_S5_S5   = -0.05*meV;  % MnII  - MnII.  
% J_14 = -0.01;         J_23 = -0.001; %Assumed zero. Only interacts via a pivalate bridge.



Sys.J = [J_AB J_A1 J_A2 J_A3 J_A4 J_B1 J_B2 J_B3 J_B4 J_12 J_13 J_14 J_23 J_24 J_34].*(-2); %-2JS.S formalism

% Stevens Operators for MnIII ions. As long as |J_S4_S4| >> |J_S4_S5| or 
% |J_S5_S5|, the lowest eigenvalues are indepedent of these values
% For scenario 1, these stay constant. Middle/green ones

DIII = -0.029*meV; % MnIII anisotropy. Free Value 
EIII = 0*DIII;   % MnIII rhombicity. For INS, assume 0.
B20III = 3*DIII; B22III = EIII; %Converting to Stevens Operator formalism

DII = -0.00*meV; % MnII anisotropy. For INS, assume 0. 
EII = DII*0;      % MnII rhombicity. For INS, assume 0.
B20II = 3*DII; B22II = EII; %Converting to Stevens Operator formalism

Sys.B2 = [B22III 0 B20III 0 0;
          B22III 0 B20III 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0;
          B22II 0 B20II 0 0];
% Sys.B4=  [0 0 0 0 -0.01*meV 0 0 0 0;
%           0 0 0 0 -0.01*meV 0 0 0 0;
%           0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0;
%           0 0 0 0 0 0 0 0 0];

%set the parameters to be varied
Vary.B2= [0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0;
          0 0 0 0 0];
% Vary.B4= Sys.B4;
% Vary.J = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]; 
Vary.J = Sys.J;   Vary.J(1) = 0;
% Vary.J([3,5,6,8])=0;
Vary.B2 = Sys.B2;
 % NumEigs = 12;    %Number of experimental eigenvalues to be fit
 % NumEigs = length(exp_eigs); 
% Exp.ev = (exp_eigs(1:NumEigs)-exp_eigs(1)).*meV;%Add eigenvalues to be used to Sys
Exp.ev = exp_eigs.*meV;

%Mint Setup 
FormFactor = {'Mn3', 'Mn3', 'Mn2', 'Mn2', 'Mn2', 'Mn2'};
% Metal ion coordinates for INS. Only relative positions are important.
Coords = [24.60550  13.06674   7.07523; % A
              22.27025  11.49336   6.95626; % B
              21.99544  13.57222   9.30268; % 1
              24.62949  10.95050   9.41655; % 2
              24.24486  10.55088   4.68547; % 3
              22.84040  14.03029   4.66904; % 4
              ];

% Sim INS powder spectrum
MintExp.SpectrumType = 'SE'; %INS intensity vs. energy

Ei = 4; %Incident neutron energy in meV
MintExp.lwfwhm = 0.02*Ei/2.355; %calculated line width full-width-at-half-max, 2 % of incident energy
% MintExp.lwfwhm = 0.05;
MintExp.Energy = linspace(-Ei*0.8,Ei*0.8,1000); %calculating the spectrum in the interval -0.8*Ei to 0.8*Ei
MintExp.Q = 0.1:0.01:2.5; %Q-range which the simulation integrates over.
MintExp.Temperature = [1.5 5 10 30];

MintOpt.NumEigs = 100; %100 eigenvalues gives a good INS sim


b =[0 0.4470 0.7410];
r=[0.8500 0.3250 0.0980];
y=[0.9290 0.6940 0.1250];
g=[0.4660 0.6740 0.1880];
colours = [b;y;g;r];


%%
[A,A0,x0,Ops,SysFixed] = Sys_Input(Sys,Vary);
   INSOperators.A =A;   INSOperators.A0 =A0;    INSOperators.x0 =x0;    
   INSOperators.Ops =Ops; INSOperators.SysFixed =SysFixed;
   INSOperators.ED = 'eigs';
%%

regulariser = @(NIter) 1e-5/NIter;
regulariser = @(NIter,X)max(max(X.J))/max(max(diag(diag(X.J'*X.J))))/sqrt(NIter);
Scaled= false;
StepTol = 1e-2;    if Scaled;StepTol = StepTol/1e3;end

% GroundStateFound = SysOut.GroundStateFound;
GroundStateFound = [];
Opt = struct('NDeflations',1,'Method','Good_GN','Linesearch','Quadratic','INSOperators',INSOperators,...
    'MaxIter',500,'StepTolerance',StepTol,'GradientTolerance',1e-0,...
    'ObjectiveTolerance',1e-1,'Scaled' ,Scaled,'epsilon',0.1,'MinAlpha',1e-1,...
    'IEPType','Difference','Verbose',true,'c1',1e-8,'SysFound',[],...
    'Regularisation',regulariser,'LinearSolver','lsqminnorm','Alpha0',1e0,...
    'GroundStateFound',GroundStateFound);

[SysOut, options] = INS_IEP(Sys,Vary,Exp,Opt);
% SysOut
% [SysOut1, options1] = INS_IEP(SysOut,Vary,Exp,Opt);




%%       
% INS-specific information. Magnetic form factors.
for i = 1:length(SysOut)
    MintSys = SysOut(i);
    % MintSys = Sys;


    MintSys.FormFactor = FormFactor;
    MintSys.Coords = Coords;
    if ~exist('MintSysPrev','var') || ~isequal(MintSys,MintSysPrev)
        clear MintOpt Eigs Vecs
    else
        if exist('Eigs','var');MintOpt.Eigs = Eigs;end
        if exist('Vecs','var');MintOpt.Vecs = Vecs;end
    end

    MintOpt.NumEigs = length(Exp.ev);


    [cross_sect,Eigs,Vecs,I_nm] = mint(MintSys,MintExp,MintOpt);
    MintSysPrev = MintSys;
    figure(i+1)
    plotmint(cross_sect, MintExp, colours)
end
%%
Kb=8.6170e-2;
T = MintExp.Temperature;
% T = 1.5;
MintExp.InitialStates = 1:MintOpt.NumEigs;
MintExp.FinalStates = 1:MintOpt.NumEigs;


IT=zeros(MintOpt.NumEigs,MintOpt.NumEigs,numel(T));
for temp_index = 1:numel(T)
    part=sum(exp(-Eigs./T(temp_index)/Kb));  % (partition function)
    for bra=MintExp.InitialStates
        for ket=MintExp.FinalStates
            IT(bra,ket,temp_index)=I_nm(bra,ket) * exp(-Eigs(bra)./T(temp_index)/Kb)/part;        
        end             
    end

end
IT = real(IT);
IT(IT<1e-16)=0;


transitionI =    @(I) I(1,2);
transitionII =  @(I) sum(I(1,[3,4]));
transitioni =   @(I) sum(I(2,[5,6,7]));
transitionii =  @(I) sum(sum(I([3,4],[5,6,7])));
transitioniii = @(I) sum(sum(I([3,4],[8,9])));


%T = 5
temp_index = 2;
ITtemp = IT(:,:,temp_index)*(108.1/IT(1,2,temp_index));

ExperimentalInt =[108.1;36.1;86.8;38.2;59.7];
SimulatedInt = [transitionI(ITtemp);transitionii(ITtemp);transitionII(ITtemp) ;transitioniii(ITtemp); transitioni(ITtemp)];
Transition = ["I";"ii";"II";"iii";"i"];
Temp5K = table(Transition,SimulatedInt,ExperimentalInt)

%%
%T = 1.5(6)
temp_index = 1;
ITtemp = IT(:,:,temp_index)*(196.2/IT(1,2,temp_index));
ExperimentalInt = [196.2;167.0;4.5;55.6];
SimulatedInt = [transitionI(ITtemp);transitionII(ITtemp) ;transitioniii(ITtemp); transitioni(ITtemp)];
Transition = ["I";"II";"iii";"i"];
Temp1point6K = table(Transition,SimulatedInt,ExperimentalInt)





% transitionI =    @(I,temp_index) I(1,2,temp_index);
% transitionTII =  @(I,temp_index) sum(I(1,[3,4],temp_index));
% transitionTi =   @(I,temp_index) sum(I(2,[5,6,7],temp_index));
% transitionTii =  @(I,temp_index) sum(sum(I([3,4],[5,6,7],temp_index)));
% transitionTiii = @(I,temp_index) sum(sum(I([3,4],[8,9],temp_index)));

%%
cross_sect=zeros(numel(T),length(E_transfer));
for temp_index = 1:numel(T)
    part=sum(exp(-Eigs./T(temp_index)/Kb));  % (partition function)
    for bra=MintExp.InitialStates
        for ket=MintExp.FinalStates
            cross_sect(temp_index,:)=cross_sect(temp_index,:)+I_nm(bra,ket) .* exp(-Eigs(bra)./T(temp_index)/Kb)/part ...
                .* exp(-((E_transfer-Eigs(ket)+Eigs(bra)).^2 )./(2*sigma^2))/sqrt(2*pi*sigma);   %Boltzman factor ('normalising' to thermal poopulation)
        end              % this bit interpolates, to give width (Gaussian) 
    end
end



%% plotting INS sim
function plotmint(cross_sect, MintExp, colours)
subplot(2,1,1)
plot(MintExp.Energy,cross_sect*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('New Sim')
ax = gca;;ax.ColorOrder = colours;
%
 subplot(2,1,2)
load("Sys0_Sim.mat")
plot(MintExp.Energy,cross_sect_Sys0*1e-4,'linewidth',1.2)
legend('1.5 K', '5 K', '10 K', '30 K')
xlim([0 1.5]);xlabel('E [meV]');ylabel('Signal');title('Old Sim')
ax = gca;;ax.ColorOrder = colours;
f.Units = 'centimeters';
f.Position = [10 10 20 25];
end
%% Save a Sys struct
DoSave = true;
SysFound = SysOut;
OptionsFound = options;
[name,val] = nextname(strcat("INS_IEP/INS Examples/Mn6/Sys Found Values/",sprintf('%0.4e',SysFound.ErrorAtDeflatedPoint)),'HZ_1','.mat');
if val>1
    for i = 1:val-1
        PreviousSaves = load(strcat("INS_IEP/INS Examples/Mn6/Sys Found Values/",sprintf('%0.4e',SysFound.ErrorAtDeflatedPoint),'HZ_',num2str(i),'.mat'));
        if isequal(SysFound,PreviousSaves.SysFound)&&isequal(OptionsFound,PreviousSaves.OptionsFound)
            warning("Sys already found and saved");
            DoSave = false;
            break
        end
    end
end
if DoSave
    save(strcat("INS_IEP/INS Examples/Mn6/Sys Found Values/",name),'SysFound','OptionsFound')
end





%%

% %%
% % [A,A0,x0,Ops,SysFixed] = Sys_Input(Sys,Vary);
% %    INSOperators.A =A;   INSOperators.A0 =A0;    INSOperators.x0 =x0;    
% %    INSOperators.Ops =Ops; INSOperators.SysFixed =SysFixed;
% %    INSOperators.ev = Exp.ev'; INSOperators.ED = 'eigs';
% 
% [~,~,x,~,~] = Sys_Input(SysOut,Vary);
% % [~,~,x,~,~] = Sys_Input(MintSys,Vary);
% %%
% EvaluateINSOperators = INSOperators;
% EvaluateINSOperators.ev = Exp.ev;
% EvaluateINSOperators.A =options.constants.A;
% if options.IEPType=="Difference"
%     IEP_Evaluate_diff(x,EvaluateINSOperators)
% elseif options.IEPType=="Classic"
%     xC = x;
%     xC(end+1) = SysOut.GroundStateFound;
%     EvaluateINSOperators.A = EvaluateINSOperators.A;
%     EvaluateINSOperators.A{end+1} = speye(size(EvaluateINSOperators.A{1}));
%     IEP_Evaluate_full(xC,EvaluateINSOperators)
% end
% 
% 
% 
% 
% 
% 
% %%
% tic %to time Hamiltonian diagonalisation
% H = ham(MintSys, [0 0 0],'sparse'); %setting up Hamiltonian
% [Vecs,E]=eigs(H,MintOpt.NumEigs,'smallestreal'); %Diagonalising Hamiltonian
% toc %to time Hamiltonian diagonalisation
% EE = diag(E);
% %Vecs = Vecs(:,idx);
% Eigs = (EE-min(EE))./meV; %convert to meV and set lowest eigenvalue to 0.
% % Eigs % print lowest eigenvalues to command window
% [Eigs(1:length(Exp.ev))*meV,Exp.ev']
% norm([Eigs(1:length(Exp.ev))*meV,Exp.ev'])
% [Eigs(1:length(Exp.ev)),Exp.ev'./meV]
% norm([Eigs(1:length(Exp.ev)),Exp.ev'./meV])
% 



% %% Calculate susceptibility using curry
% 
% %NB! This code in theory calculates the susceptibility. However, the
% %Hamiltonian is too large and a normal computer cannot run this snippet
% clear MintExp
% clear MintOpt
% MintOpt.Units = 'CGS';
% MintOpt.Output = 'chimol';
% MintExp.Temperature = 1:1:300;
% MintExp.Field = 100; %mT
% chi = curry(SysOut,MintExp,MintOpt);
% 
% figure()
% plot(Exp.Temperature, chi.*Exp.Temperature,'o')
% xlabel('T [K]')
% ylabel('\chi T [cm^3 K/mol]')
% 
% figure()
% plot(Exp.Temperature, chi,'o')
% xlabel('T [K]')
% ylabel('\chi [cm^3/mol]')



% 0.01*meV is acceptable error
%
%Don't use any more parameters
%Fit intensity  - meet mike
%Range on each of the parameters
%
% 1% of teh magnitude of the value


% [0,0.8]