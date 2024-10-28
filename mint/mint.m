function [cross_sect,Eigs,Vecs,I_nm] = mint(Sys,Exp,Opt)
% MINT - Simulate the inelastic neutron scattering of magnetic ions and molecules.
%
% [cross_sect,Eigs,Vecs] = mint(Sys,Exp,Opt)
% cross_sect = mint(Sys,Exp)
%
% Standard <a href="matlab: 
% web('https://easyspin.org/')">easyspin</a> (Sys, Exp, Opt) type input structures are required
% with neutron scattering specific fields:
%
% Exp
%-----
% Exp.SpectrumType, this can be either:
% 'SE'   (energy spectrum with integrated Q for powder samples)
% 'SQ'   (Q dependence with integrated energy for powder samples)
% 'SQE'  (2D energy versus Q spectrum for powder samples)
% 'SQxyz (3D Q dependence in x, y z directions for single crystal samples).
%
% Exp.Energy is an array, the INS energy range in meV (for SE, SQE).
%
% Exp.lwfwhm is gaussian broadening FWHM of INS peaks in meV (for SE, SQE).
%
% Exp.Q is an array, the INS momentum transfer in A^-1 (for SQ, SE, SQE).
%
% Exp.Qx, Exp.Qy, Exp.Qz are arrays for INS momentum transfer in A^-1 (for SQxyz). 
%
% Exp.InitialStates is an array, initial states included in the calculation.
% Exp.FinalStates is an array, final states included in the calculation.
% If InitialStates or FinalStates are not given they are calulated based on 
% Exp.Temperature and the Exp.Energy range provided.
%
% Sys
%-----
% Sys.FormFactor is the magnetic formfactor, given by the element symbol
% followed by a number representing the oxidation state for heterometallic 
% compounds a cell array is given, i.e. Sys.FormFactor = {'Tb3','Cu2'}
% Transition metal and lanthanide ions included for typical oxidation
% states. If the required element / oxidation state is not available, give a 1x7 array
% for <jo> of form  [A a B b C c D].
%
% Sys.Coords = [site1x site1y site1z ; site2x site2y site2z; etc];
% to simulate the INS of a metal ion cluster you need to know relative the
% atomic cartesian coordinates of each metal sites in Angstroms.
% The atomic coordinates are in the same reference frame as the spin system.
% (Note - that Sys.J is +J S1S2)
%
% Opt
%-----
% Opt.NumEigs is an integer that caps the number of eigenvalues considered 
% in the INS simulation, suitable when only the lowest eigenvalues are accessed 
% experimentally. This is relevant when computing large spin systems with 
% Hamiltonian matrix dimensions that are prohibitively large to evaluate in
% their entirety. 
% 
% Dr Michael L. Baker
% Department of Chemistry, The University of Manchester, UK.
% michael.baker@manchester.ac.uk
% For more information, see <a href="matlab: 
% web('https://www.mlbakerlab.co.uk')">the M. L. Baker lab Web site</a>.
% 
tic

% some useful constants and conversions
rcm = 29979.2458;   %rcm    to MHz
meV = rcm*8.065;    %meV    to MHz
Kb  = 8.6170e-2;    % in meV/K
% *********
% Albans options

if exist('Opt','var') && isfield(Opt,'ShowMessages') && Opt.ShowMessages 
    verbose = true;
else 
    verbose = false;
end


if exist('Opt','var') && isfield(Opt,'Eigs') && isfield(Opt,'Vecs') && length(Opt.Eigs)>=Opt.NumEigs
    Eigs = Opt.Eigs; Vecs = Opt.Vecs;
    if length(Eigs)~=size(Vecs,2) 
        error('Please input the corresponding eigenvectors to the given eigenvalues')
    end
    disp('Using supplied eigenvalues and eigenvectors')
else
    disp('Hamiltonian Diagonalization in progress...')
if exist('Opt','var') && isfield(Opt,'NumEigs')
    H = ham(Sys,[0,0,0],'sparse');
    
    [Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal');
    [~,lw] = lastwarn;
    if lw=='MATLAB:eigs:NotAllEigsConverged'
     [Vecs,E]=eigs(H,Opt.NumEigs,'smallestreal','subspacedimension',3*(Opt.NumEigs+20));
    end
else
    H = ham(Sys,[0,0,0]);
    [Vecs,E]=eig(H);
end
EE = diag(E);
Eigs = (EE-EE(1))./meV; %convert to meV
disp('Diagonalization done.')
end



if verbose
disp('**************mint: simulate INS spectra*****************')
disp('    Michael L. Baker (michael.baker@manchester.ac.uk)')
disp('*********************************************************')
end
% some user input error checks and corrections
if ~isfield(Exp,'SpectrumType')
    % assume energy spectrum
    Exp.SpectrumType = 'SE';
end

% for calculations other than SQxyz and SQ check energy array versus peak width are resonable
if ~strcmp(Exp.SpectrumType,'SQxyz') && ~strcmp(Exp.SpectrumType,'SQ')
    if Exp.lwfwhm < min(diff(Exp.Energy))
        error('Energy resolution (Exp.lwfwhm) must be greater than energy step in Exp.Energy')
    end
end

if ~isfield(Sys,'Coords')
    if numel(Sys.S) == 1
        Sys.Coords = [0 0 0];
    else
        error('Coordinates for metal ion x y z positions must be given in Angstroms (Sys.Coords)')
    end
    
end

% check user input and get formfactor for INS calc 
if isfield(Sys,'FormFactor')
    if iscell(Sys.FormFactor) ~=1 % if not a cell make it a cell
        Sys.FormFactor = {Sys.FormFactor};
    end
    if numel(Sys.FormFactor) == numel(Sys.S)
        Fj0= formfactor(Sys.FormFactor);
    elseif numel(Sys.FormFactor) == 1
        for ii=1:1:numel(Sys.S); ion{ii} = Sys.FormFactor{1}; end
        Fj0= formfactor(ion);
    else
        error('The formfactor for each spin centre must be included')
    end
else
    error('No formfactor given, please define with Sys.FormFactor')
end

%coordinates
R_relative=zeros(numel(Sys.S),numel(Sys.S),3);
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        R_relative(ii,jj,:)=Sys.Coords(ii,:)-Sys.Coords(jj,:);
    end
end
abs_R=zeros(numel(Sys.S),numel(Sys.S));
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        abs_R(ii,jj)=norm([R_relative(ii,jj,1) R_relative(ii,jj,2) R_relative(ii,jj,3)]);
    end
end



% if Inital and Final states are not given estimate them based on
% population and energy spectrum range

if ~isfield(Exp,'InitialStates')
    Exp.InitialStates = 1:numel(Eigs);
    % reduce initial states to those that are significantly populated.
    part=sum(exp(-Eigs./max(Exp.Temperature)/Kb));
    pop = exp(-Eigs./max(Exp.Temperature)/Kb)/part;
    populated_indexes = find(pop > 0.005);
    
    if populated_indexes(end) < max(Exp.InitialStates)
        if verbose
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        disp(['Only the first ',num2str(populated_indexes(end)),' states have a population '])
        disp(['greater than 0.5% at ',num2str(max(Exp.Temperature)),' K.'])
        disp(['If there are more than ', num2str(populated_indexes(end)),' states the calculation'])
        disp(['will be sped up by reducing to Exp.InitialStates = 1:',num2str(populated_indexes(end))])
        disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        end
        Exp.InitialStates = 1:populated_indexes(end);
    elseif numel(Exp.InitialStates) > numel(Eigs)
        disp(['total number of states = ',num2str(numel(Eigs))])
        error('inital states given (Exp.InitialStates) exeeds total number of states:')
    end
end

if ~isfield(Exp,'FinalStates')
    Exp.FinalStates = 2:numel(Eigs);
    % reduce final states to be within the energy spectrum range
    % maximum transition energy is max(abs(Exp.Energy))
    [~,final_state_index2]=min( abs( Eigs-(max(abs(Exp.Energy))+Eigs(Exp.InitialStates(end))) ) );
    % first final state is set to be greater than zero energy transfer
    final_state_index1 = find(Eigs ~= 0,1,'first');
    Exp.FinalStates = final_state_index1:final_state_index2;
    if verbose
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp('Transition final states evaluated are')
    disp(['Exp.FinalStates = ',num2str(final_state_index1),':',num2str(final_state_index2)])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    end
elseif exist('Opt','var') && isfield(Opt,'NumEigs')
    % cut down final states to match reduced NumEigs if required
    if max(Exp.FinalStates) > Opt.NumEigs
        Exp.FinalStates = Exp.FinalStates(find(Exp.FinalStates <= Opt.NumEigs));
    end
end
if verbose; disp('Computing <n|s|m> for all spin operators...');end
[Sx_rot,Sy_rot,Sz_rot] = bra_n_s_m_ket(Sys,Vecs);
if verbose; disp('<n|s|m> done.');end

switch Exp.SpectrumType
    case 'SE'
        if verbose;disp('Calculating S(Energy)...');end
        [cross_sect, I_nm] = SE(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot);
        if verbose; disp('S(Energy) done.'); end
    case 'SQE'
        disp('Calculating S(Q,Energy)...')
        [cross_sect, I_nm] = SQE(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot);
        disp('S(Q,Energy) done.')
    case 'SQ'
        disp('Calculating S(Q)...')
        [cross_sect, I_nm] = SQ(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot);
        disp('S(Q) done.')
    case 'SQxyz'
        disp('Calculating S(Qxyz)...')
        [cross_sect, I_nm] = SQxyz(Sys,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot);
        disp('S(Qxyz) done.')
    otherwise
        disp('requested spectrum type is not recognised')
end

% rounding erorr can generate small imagninary component to intensity.
% check for this, report the magnetude, and then remove it.
if max(max(imag(cross_sect))) > 0
    disp(['Imaginary component to intensity found of maximum magnetude ',num2str(max(max(imag(cross_sect)))),'. This warning can be ignored if this value is small'])
    cross_sect = real(cross_sect);
else
    cross_sect = real(cross_sect);
end
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cross_sect, I_nm] = SQ(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot)

Kb=8.6170e-2;              % meV per K
T=Exp.Temperature;

% reduced version using selected spin operator elements
bra=Exp.InitialStates;
ket=Exp.FinalStates;
I_nm=zeros(numel(bra),numel(ket),numel(Exp.Q));


for Qindexer=1:numel(Exp.Q)
    Q = Exp.Q(Qindexer);
    s=Q/4/pi;

    
    for ii=1:numel(R_relative(:,1,1))
        for jj=1:numel(R_relative(:,1,1))
            x=Q*abs_R(ii,jj);
            if ii==jj
                C20=0;
                C22=0;
                C2_2=0;
                C21=0;
                C2_1=0;
                j_0=1;
                j_2=0;
            else
                C20=(3*(R_relative(ii,jj,3)/abs_R(ii,jj))^2 - 1)/2;
                C22=(R_relative(ii,jj,1)^2 - R_relative(ii,jj,2)^2)/abs_R(ii,jj)^2;
                C2_2=2*R_relative(ii,jj,1)*R_relative(ii,jj,2)/abs_R(ii,jj)^2;
                C21=R_relative(ii,jj,1)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                C2_1=R_relative(ii,jj,2)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                j_0=sin(x)/x;
                j_2=(3/x^2 - 1)*sin(x)/x - 3*cos(x)/x^2;
            end
            
            Fa=Fj0(ii,1)*exp(-Fj0(ii,2)*s^2)+Fj0(ii,3)*exp(-Fj0(ii,4)*s^2)+Fj0(ii,5)*exp(-Fj0(ii,6)*s^2)+Fj0(ii,7);
            Fb=Fj0(jj,1)*exp(-Fj0(jj,2)*s^2)+Fj0(jj,3)*exp(-Fj0(jj,4)*s^2)+Fj0(jj,5)*exp(-Fj0(jj,6)*s^2)+Fj0(jj,7);
            
            I_nm(1:numel(bra),1:numel(ket),Qindexer) = I_nm(1:numel(bra),1:numel(ket),Qindexer) + abs(Fa*Fb) *( (j_0+C20*j_2) * conj(Sz_rot(bra,ket,ii)).*Sz_rot(bra,ket,jj)*2/3 ...
                + (j_0-C20*j_2/2) * ( conj(Sx_rot(bra,ket,ii)).*Sx_rot(bra,ket,jj) + conj(Sy_rot(bra,ket,ii)).*Sy_rot(bra,ket,jj) )*2/3 ...
                + j_2 * ( C22*( conj(Sx_rot(bra,ket,ii)).*Sx_rot(bra,ket,jj) - conj(Sy_rot(bra,ket,ii)).*Sy_rot(bra,ket,jj) ) + C2_2*( conj(Sx_rot(bra,ket,ii)).*Sy_rot(bra,ket,jj) - conj(Sy_rot(bra,ket,ii)).*Sx_rot(bra,ket,jj) ) )/2 ...
                + j_2 * ( C21*( conj(Sz_rot(bra,ket,ii)).*Sx_rot(bra,ket,jj) - conj(Sx_rot(bra,ket,ii)).*Sz_rot(bra,ket,jj) ) + C2_1*( conj(Sz_rot(bra,ket,ii)).*Sy_rot(bra,ket,jj) - conj(Sy_rot(bra,ket,ii)).*Sz_rot(bra,ket,jj) ) ) );
            
        end
    end
end

cross_sect=zeros(numel(T),length(Exp.Q));
for temp_index = 1:numel(T)
    for Q_index = 1:numel(Exp.Q)
        part=sum(exp(-Eigs./T(temp_index)/Kb));
        for bra_index=1:numel(bra)
            for ket_index=1:numel(ket)
                cross_sect(temp_index,Q_index)=cross_sect(temp_index,Q_index)+I_nm(bra_index,ket_index,Q_index) .* exp(-Eigs(bra(bra_index))./T(temp_index)/Kb)/part;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cross_sect, I_nm] = SE(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot,varargin)
if isempty(varargin)
    verbose = false;
else
    verbose = varargin{1};
end
Kb=8.6170e-2; % meV per K
T=Exp.Temperature;
E_transfer = Exp.Energy;
sigma=Exp.lwfwhm/2.35;
I_nm=zeros(numel(Eigs),numel(Eigs));
if verbose; disp(' Calculating transition matrix elements'); end
for Q=Exp.Q
    
    s=Q/4/pi;
        
    for ii=1:numel(R_relative(:,1,1))
        for jj=1:numel(R_relative(:,1,1))
            % PHYSICAL REVIEW B 71, 174407 2005
            x=Q*abs_R(ii,jj);
            if ii==jj
                C20=0;
                C22=0;
                C2_2=0;
                C21=0;
                C2_1=0;
                j_0=1;
                j_2=0;
            else
                C20=(3*(R_relative(ii,jj,3)/abs_R(ii,jj))^2 - 1)/2;
                C22=(R_relative(ii,jj,1)^2 - R_relative(ii,jj,2)^2)/abs_R(ii,jj)^2;
                C2_2=2*R_relative(ii,jj,1)*R_relative(ii,jj,2)/abs_R(ii,jj)^2;
                C21=R_relative(ii,jj,1)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                C2_1=R_relative(ii,jj,2)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                j_0=sin(x)/x;
                j_2=(3/x^2 - 1)*sin(x)/x - 3*cos(x)/x^2;
            end
            
            Fa=Fj0(ii,1)*exp(-Fj0(ii,2)*s^2)+Fj0(ii,3)*exp(-Fj0(ii,4)*s^2)+Fj0(ii,5)*exp(-Fj0(ii,6)*s^2)+Fj0(ii,7);
            Fb=Fj0(jj,1)*exp(-Fj0(jj,2)*s^2)+Fj0(jj,3)*exp(-Fj0(jj,4)*s^2)+Fj0(jj,5)*exp(-Fj0(jj,6)*s^2)+Fj0(jj,7);
            I_nm = I_nm + abs(Fa*Fb) *( (j_0+C20*j_2) * conj(Sz_rot(:,:,ii)).*Sz_rot(:,:,jj)*2/3 ...
                + (j_0-C20*j_2/2) * ( conj(Sx_rot(:,:,ii)).*Sx_rot(:,:,jj) + conj(Sy_rot(:,:,ii)).*Sy_rot(:,:,jj) )*2/3 ...
                + j_2 * ( C22*( conj(Sx_rot(:,:,ii)).*Sx_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sy_rot(:,:,jj) ) + C2_2*( conj(Sx_rot(:,:,ii)).*Sy_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sx_rot(:,:,jj) ) )/2 ...
                + j_2 * ( C21*( conj(Sz_rot(:,:,ii)).*Sx_rot(:,:,jj) - conj(Sx_rot(:,:,ii)).*Sz_rot(:,:,jj) ) + C2_1*( conj(Sz_rot(:,:,ii)).*Sy_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sz_rot(:,:,jj) ) ) );
            
        end
    end
end
if verbose;disp(' Calculating INS cross section');end
cross_sect=zeros(numel(T),length(E_transfer));
for temp_index = 1:numel(T)
    part=sum(exp(-Eigs./T(temp_index)/Kb));  % (partition function)
    for bra=Exp.InitialStates
        for ket=Exp.FinalStates
            cross_sect(temp_index,:)=cross_sect(temp_index,:)+I_nm(bra,ket) .* exp(-Eigs(bra)./T(temp_index)/Kb)/part ...
                .* exp(-( (E_transfer-Eigs(ket)+Eigs(bra)).^2 )./(2*sigma^2))/sqrt(2*pi*sigma);   %Boltzman factor ('normalising' to thermal poopulation)
        end              % this bit interpolates, to give width (Gaussian) 
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cross_sect, I_nm] = SQE(R_relative,abs_R,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot)

Kb=8.6170e-2; % meV per K
T=Exp.Temperature(1);
E_transfer = Exp.Energy;
sigma=Exp.lwfwhm/2.355;

I_nm=zeros(numel(Eigs),numel(Eigs),numel(Exp.Q));

for Qindexer=1:numel(Exp.Q)
    Q = Exp.Q(Qindexer);
    s=Q/4/pi;
    
    for ii=1:numel(R_relative(:,1,1))
        for jj=1:numel(R_relative(:,1,1))
            x=Q*abs_R(ii,jj);
            if ii==jj
                C20=0;
                C22=0;
                C2_2=0;
                C21=0;
                C2_1=0;
                j_0=1;
                j_2=0;
            else
                C20=(3*(R_relative(ii,jj,3)/abs_R(ii,jj))^2 - 1)/2;
                C22=(R_relative(ii,jj,1)^2 - R_relative(ii,jj,2)^2)/abs_R(ii,jj)^2;
                C2_2=2*R_relative(ii,jj,1)*R_relative(ii,jj,2)/abs_R(ii,jj)^2;
                C21=R_relative(ii,jj,1)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                C2_1=R_relative(ii,jj,2)*R_relative(ii,jj,3)/abs_R(ii,jj)^2;
                j_0=sin(x)/x;
                j_2=(3/x^2 - 1)*sin(x)/x - 3*cos(x)/x^2;
            end
            
            Fa=Fj0(ii,1)*exp(-Fj0(ii,2)*s^2)+Fj0(ii,3)*exp(-Fj0(ii,4)*s^2)+Fj0(ii,5)*exp(-Fj0(ii,6)*s^2)+Fj0(ii,7);
            Fb=Fj0(jj,1)*exp(-Fj0(jj,2)*s^2)+Fj0(jj,3)*exp(-Fj0(jj,4)*s^2)+Fj0(jj,5)*exp(-Fj0(jj,6)*s^2)+Fj0(jj,7);
                   
            I_nm(:,:,Qindexer) = I_nm(:,:,Qindexer) + abs(Fa*Fb) *( (j_0+C20*j_2) * conj(Sz_rot(:,:,ii)).*Sz_rot(:,:,jj)*2/3 ...
                + (j_0-C20*j_2/2) * ( conj(Sx_rot(:,:,ii)).*Sx_rot(:,:,jj) + conj(Sy_rot(:,:,ii)).*Sy_rot(:,:,jj) )*2/3 ...
                + j_2 * ( C22*( conj(Sx_rot(:,:,ii)).*Sx_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sy_rot(:,:,jj) ) + C2_2*( conj(Sx_rot(:,:,ii)).*Sy_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sx_rot(:,:,jj) ) )/2 ...
                + j_2 * ( C21*( conj(Sz_rot(:,:,ii)).*Sx_rot(:,:,jj) - conj(Sx_rot(:,:,ii)).*Sz_rot(:,:,jj) ) + C2_1*( conj(Sz_rot(:,:,ii)).*Sy_rot(:,:,jj) - conj(Sy_rot(:,:,ii)).*Sz_rot(:,:,jj) ) ) );
            
        end
    end
end

cross_sect=zeros(length(E_transfer),numel(Exp.Q));
for Q_index = 1:numel(Exp.Q)
    part=sum(exp(-Eigs./T/Kb));
    for bra=Exp.InitialStates
        for ket=Exp.FinalStates
            cross_sect(:,Q_index)=cross_sect(:,Q_index)+I_nm(bra,ket,Q_index) .* exp(-Eigs(bra)./T/Kb)/part ...
                * exp(-( (E_transfer'-Eigs(ket)+Eigs(bra)).^2 )./(2*sigma^2))/sqrt(2*pi*sigma);
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cross_sect, I_nm] = SE_ave(Sys,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot)
error('the numerical intergration method is inaccurate')
% Q vectors
phi_steps = 5;
theta_steps = 10;
phi = linspace(0,pi-pi/phi_steps,phi_steps);
theta = linspace(0,2*pi-2*pi/theta_steps,theta_steps);
Qvecs_index =1;
for ii=1:theta_steps
    for jj=1:phi_steps
        Qvecs(Qvecs_index,1:3) = ang2vec(phi(jj),theta(ii))';
        Qvecs_index = Qvecs_index + 1;
    end
end
Qvecs = sphgrid('Ci',91)';

index_Qxyz = 1;
for ii=1:numel(Exp.Q)
    Qxyz(index_Qxyz:index_Qxyz-1+numel(Qvecs(:,1)),1:3) = Exp.Q(ii).*Qvecs;
    index_Qxyz = index_Qxyz + numel(Qvecs(:,1));
end
modQ=sqrt(Qxyz(:,1).^2+Qxyz(:,2).^2+Qxyz(:,3).^2);

%relative coordinates
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        R_relative_xyz(1:3,ii,jj) = Sys.Coords(ii,:)-Sys.Coords(jj,:);
    end
end

%formfactor coeffs
s=modQ/4/pi;

% I_nm
I_nm=zeros(numel(Eigs),numel(Eigs));   % Export this
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        for kk=1:numel(Qxyz(:,1))
            
            s=modQ(kk)/4/pi;
            Fa=Fj0(ii,1)*exp(-Fj0(ii,2)*s^2)+Fj0(ii,3)*exp(-Fj0(ii,4)*s^2)+Fj0(ii,5)*exp(-Fj0(ii,6)*s^2)+Fj0(ii,7);
            Fb=Fj0(jj,1)*exp(-Fj0(jj,2)*s^2)+Fj0(jj,3)*exp(-Fj0(jj,4)*s^2)+Fj0(jj,5)*exp(-Fj0(jj,6)*s^2)+Fj0(jj,7);

            I_nm = I_nm + (1-Qxyz(kk,1).^2./(modQ(kk).^2)) * abs(Fa*Fb) * exp(sqrt(-1)*Qxyz(kk,1:3)*R_relative_xyz(:,ii,jj)) * conj(Sx_rot(:,:,ii)).*Sx_rot(:,:,jj) + ...
                (1-Qxyz(kk,2).^2./(modQ(kk).^2)) * abs(F)^2 * exp(sqrt(-1)*Qxyz(kk,1:3)*R_relative_xyz(:,ii,jj)) * conj(Sy_rot(:,:,ii)).*Sy_rot(:,:,jj) + ...
                (1-Qxyz(kk,3).^2./(modQ(kk).^2)) * abs(F)^2 * exp(sqrt(-1)*Qxyz(kk,1:3)*R_relative_xyz(:,ii,jj)) * conj(Sz_rot(:,:,ii)).*Sz_rot(:,:,jj);
        end
    end
end

Kb=8.6170e-2;              % meV per K
T=Exp.Temperature;
E_transfer = Exp.Energy;
sigma=Exp.lwfwhm/2.35;
Exp.InitialStates = 1:numel(Eigs);
Exp.FinalStates = 1:numel(Eigs);

cross_sect=zeros(numel(T),length(E_transfer));
for temp_index = 1:numel(T)
    part=sum(exp(-Eigs./T(temp_index)/Kb));    % (partition function)
    for bra=Exp.InitialStates
        for ket=Exp.FinalStates
            cross_sect(temp_index,:)=cross_sect(temp_index,:)+I_nm(bra,ket) .* exp(-Eigs(bra)./T(temp_index)/Kb)/part ...
                * exp(-( (E_transfer-Eigs(ket)+Eigs(bra)).^2 )./(2*sigma^2))/sqrt(2*pi*sigma);  %Boltzman factor ('normalising' to thermal poopulation)
        end              % this bit interpolates, to give width (Gaussian) 
    end
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cross_sect, I_nm] = SQxyz(Sys,Exp,Fj0,Eigs,Sx_rot,Sy_rot,Sz_rot)

bra=Exp.InitialStates;
ket=Exp.FinalStates;

%relative coordinates
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        R_relative_xyz(1:3,ii,jj) = Sys.Coords(ii,:)-Sys.Coords(jj,:);
    end
end

% I_nm
I_nm=zeros(numel(bra),numel(ket),numel(Exp.Qx),numel(Exp.Qy),numel(Exp.Qz));
for ii=1:numel(Sys.S)
    for jj=1:numel(Sys.S)
        for xx=1:numel(Exp.Qx)
            for yy=1:numel(Exp.Qy)
                for zz=1:numel(Exp.Qz)
            modQ = sqrt(Exp.Qx(xx).^2+Exp.Qy(yy).^2+Exp.Qz(zz).^2);
            s=modQ/4/pi;
            Fa=Fj0(ii,1)*exp(-Fj0(ii,2)*s^2)+Fj0(ii,3)*exp(-Fj0(ii,4)*s^2)+Fj0(ii,5)*exp(-Fj0(ii,6)*s^2)+Fj0(ii,7);
            Fb=Fj0(jj,1)*exp(-Fj0(jj,2)*s^2)+Fj0(jj,3)*exp(-Fj0(jj,4)*s^2)+Fj0(jj,5)*exp(-Fj0(jj,6)*s^2)+Fj0(jj,7);

            I_nm(1:numel(bra),1:numel(ket),xx,yy,zz) = I_nm(:,:,xx,yy,zz) + ...
                    (1-Exp.Qx(xx).^2./(modQ.^2)) * abs(Fa*Fb) * exp(sqrt(-1)*[Exp.Qx(xx),Exp.Qy(yy),Exp.Qz(zz)]*R_relative_xyz(:,ii,jj)) * conj(Sx_rot(bra,ket,ii)).*Sx_rot(bra,ket,jj) + ...
                    (1-Exp.Qy(yy).^2./(modQ.^2)) * abs(Fa*Fb) * exp(sqrt(-1)*[Exp.Qx(xx),Exp.Qy(yy),Exp.Qz(zz)]*R_relative_xyz(:,ii,jj)) * conj(Sy_rot(bra,ket,ii)).*Sy_rot(bra,ket,jj) + ...
                    (1-Exp.Qz(zz).^2./(modQ.^2)) * abs(Fa*Fb) * exp(sqrt(-1)*[Exp.Qx(xx),Exp.Qy(yy),Exp.Qz(zz)]*R_relative_xyz(:,ii,jj)) * conj(Sz_rot(bra,ket,ii)).*Sz_rot(bra,ket,jj);
                end
            end
        end
    end
end

Kb=8.6170e-2;              % meV per K
T=Exp.Temperature(1);

cross_sect=zeros(numel(Exp.Qx),numel(Exp.Qy),numel(Exp.Qz));
bra=Exp.InitialStates;
ket=Exp.FinalStates;
part=sum(exp(-Eigs./T/Kb));
for Qx_index = 1:numel(Exp.Qx)
    for Qy_index = 1:numel(Exp.Qy)
        for Qz_index = 1:numel(Exp.Qz)
            for bra_index=1:numel(bra)
                for ket_index=1:numel(ket)
                    cross_sect(Qx_index,Qy_index,Qz_index)=cross_sect(Qx_index,Qy_index,Qz_index)+I_nm(bra_index,ket_index,Qx_index,Qy_index,Qz_index) .* exp(-Eigs(bra(bra_index))./T/Kb)/part;
                end
            end
        end
    end
end

% reduce and rotate for surf(Qdim1,Qdim2,cross_sect) format
cross_sect=rot90(squeeze(cross_sect));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sx_rot,Sy_rot,Sz_rot] = bra_n_s_m_ket(Sys,Vecs)


for ii=1:numel(Sys.S)
 LHScell{ii}=['S',num2str(ii),'x,' 'S',num2str(ii),'y,' 'S',num2str(ii),'z,'];
 RHScell{ii}=['''','x',num2str(ii),''',''y',num2str(ii),''',''z',num2str(ii),'''',','];
end
LHS = [LHScell{:}]; LHS(end)=[];
RHS = [RHScell{:}];RHS = RHS(1:end-1);

eval(['[',LHS,']=sop(Sys,',RHS,',''sparse''',');'])

% (' = Complex conjugate transpose)
for ii=1:numel(Sys.S)
    Sx_rot(:,:,ii)=Vecs'*eval(['S',num2str(ii),'x'])*Vecs;
    Sy_rot(:,:,ii)=Vecs'*eval(['S',num2str(ii),'y'])*Vecs;
    Sz_rot(:,:,ii)=Vecs'*eval(['S',num2str(ii),'z'])*Vecs;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  F = formfactor(ion)

for ii=1:numel(ion)
    
    if isnumeric(ion{ii})
        F(ii,1:7)=ion{ii};
    else
        
        switch ion{ii}
            case 'Sc0'
                F(ii,1:7) = [0.2512 	90.0296 	0.3290 	39.4021 	0.4235 	14.3222 	-0.0043];
            case 'Sc1'
                F(ii,1:7) = [0.4889 	51.1603 	0.5203 	14.0764 	-0.0286 	0.1792 	0.0185];
            case 'Sc2'
                F(ii,1:7) = [0.5048 	31.4035 	0.5186 	10.9897 	-0.0241 	1.1831 	0.0000];
            case 'Ti0'
                F(ii,1:7) = [0.4657 	33.5898 	0.5490 	9.8791 	-0.0291 	0.3232 	0.0123];
            case 'Ti1'
                F(ii,1:7) = [0.5093 	36.7033 	0.5032 	10.3713 	-0.0263 	0.3106 	0.0116];
            case 'Ti2'
                F(ii,1:7) = [0.5091 	24.9763 	0.5162 	8.7569 	-0.0281 	0.9160 	0.0015];
            case 'Ti3'
                F(ii,1:7) = [0.3571 	22.8413 	0.6688 	8.9306 	-0.0354 	0.4833 	0.0099];
            case 'V0'
                F(ii,1:7) = [0.4086 	28.8109 	0.6077 	8.5437 	-0.0295 	0.2768 	0.0123];
            case 'V1'
                F(ii,1:7) = [0.4444 	32.6479 	0.5683 	9.0971 	-0.2285 	0.0218 	0.2150];
            case 'V2'
                F(ii,1:7) = [0.4085 	23.8526 	0.6091 	8.2456 	-0.1676 	0.0415 	0.1496];
            case 'V3'
                F(ii,1:7) = [0.3598 	19.3364 	0.6632 	7.6172 	-0.3064 	0.0296 	0.2835];
            case 'V4'
                F(ii,1:7) = [0.3106 	16.8160 	0.7198 	7.0487 	-0.0521 	0.3020 	0.0221];
            case 'Cr0'
                F(ii,1:7) = [0.1135 	45.1990 	0.3481 	19.4931 	0.5477 	7.3542 	-0.0092];
            case 'Cr1'
                F(ii,1:7) = [-0.0977 	0.0470 	0.4544 	26.0054 	0.5579 	7.4892 	0.0831];
            case 'Cr2'
                F(ii,1:7) = [1.2024 	-0.0055 	0.4158 	20.5475 	0.6032 	6.9560 	-1.2218];
            case 'Cr3'
                F(ii,1:7) = [-0.3094 	0.0274 	0.3680 	17.0355 	0.6559 	6.5236 	0.2856];
            case 'Cr4'
                F(ii,1:7) = [-0.2320 	0.0433 	0.3101 	14.9518 	0.7182 	6.1726 	0.2042];
            case 'Mn0'
                F(ii,1:7) = [0.2438 	24.9629 	0.1472 	15.6728 	0.6189 	6.5403 	-0.0105];
            case 'Mn1'
                F(ii,1:7) = [-0.0138 	0.4213 	0.4231 	24.6680 	0.5905 	6.6545 	-0.0010];
            case 'Mn2'
                F(ii,1:7) = [0.4220 	17.6840 	0.5948 	6.0050 	0.0043 	-0.6090 	-0.0219];
            case 'Mn3'
                F(ii,1:7) = [0.4198 	14.2829 	0.6054 	5.4689 	0.9241 	-0.0088 	-0.9498];
            case 'Mn4'
                F(ii,1:7) = [0.3760 	12.5661 	0.6602 	5.1329 	-0.0372 	0.5630 	0.0011];
            case 'Fe0'
                F(ii,1:7) = [0.0706 	35.0085 	0.3589 	15.3583 	0.5819 	5.5606 	-0.0114];
            case 'Fe1'
                F(ii,1:7) = [0.1251 	34.9633 	0.3629 	15.5144 	0.5223 	5.5914 	-0.0105];
            case 'Fe2'
                F(ii,1:7) = [0.0263 	34.9597 	0.3668 	15.9435 	0.6188 	5.5935 	-0.0119];
            case 'Fe3'
                F(ii,1:7) = [0.3972 	13.2442 	0.6295 	4.9034 	-0.0314 	0.3496 	0.0044 ];
            case 'Fe4'
                F(ii,1:7) = [0.3782 	11.3800 	0.6556 	4.5920 	-0.0346 	0.4833 	0.0005];
            case 'Co0'
                F(ii,1:7) = [0.4139 	16.1616 	0.6013 	4.7805 	-0.1518 	0.0210 	0.1345];
            case 'Co1'
                F(ii,1:7) = [0.0990 	33.1252 	0.3645 	15.1768 	0.5470 	5.0081 	-0.0109];
            case 'Co2'
                F(ii,1:7) = [0.4332 	14.3553 	0.5857 	4.6077 	-0.0382 	0.1338 	0.0179];
            case 'Co3'
                F(ii,1:7) = [0.3902 	12.5078 	0.6324 	4.4574 	-0.1500 	0.0343 	0.1272];
            case 'Co4'
                F(ii,1:7) = [0.3515 	10.7785 	0.6778 	4.2343 	-0.0389 	0.2409 	0.0098  ];
            case 'Ni0'
                F(ii,1:7) = [-0.0172 	35.7392 	0.3174 	14.2689 	0.7136 	4.5661 	-0.0143];
            case 'Ni1'
                F(ii,1:7) = [0.0705 	35.8561 	0.3984 	13.8042 	0.5427 	4.3965 	-0.0118];
            case 'Ni2'
                F(ii,1:7) = [0.0163 	35.8826 	0.3916 	13.2233 	0.6052 	4.3388 	-0.0133];
            case 'Ni3'
                F(ii,1:7) = [-0.0134 	35.8677 	0.2678 	12.3326 	0.7614 	4.2369 	-0.0162];
            case 'Ni4'
                F(ii,1:7) = [-0.0090 	35.8614 	0.2776 	11.7904 	0.7474 	4.2011 	-0.0163];
            case 'Cu0'
                F(ii,1:7) = [0.0909 	34.9838 	0.4088 	11.4432 	0.5128 	3.8248 	-0.0124];
            case 'Cu1'
                F(ii,1:7) = [0.0749 	34.9656 	0.4147 	11.7642 	0.5238 	3.8497 	-0.0127];
            case 'Cu2'
                F(ii,1:7) = [0.0232 	34.9686 	0.4023 	11.5640 	0.5882 	3.8428 	-0.0137];
            case 'Cu3'
                F(ii,1:7) = [0.0031 	34.9074 	0.3582 	10.9138 	0.6531 	3.8279 	-0.0147];
            case 'Cu4'
                F(ii,1:7) = [-0.0132 	30.6817 	0.2801 	11.1626 	0.7490 	3.8172 	-0.0165];
            case 'Ce2'
                F(ii,1:7) = [0.2953 	17.6846 	0.2923 	6.7329 	0.4313 	5.3827 	-0.0194];
            case 'Nd2'
                F(ii,1:7) = [0.1645 	25.0453 	0.2522 	11.9782 	0.6012 	4.9461 	-0.0180];
            case 'Nd3'
                F(ii,1:7) = [0.0540 	25.0293 	0.3101 	12.1020 	0.6575 	4.7223 	-0.0216];
            case 'Sm2'
                F(ii,1:7) = [0.0909 	25.2032 	0.3037 	11.8562 	0.6250 	4.2366 	-0.0200];
            case 'Sm3'
                F(ii,1:7) = [0.0288 	25.2068 	0.2973 	11.8311 	0.6954 	4.2117 	-0.0213];
            case 'Eu2'
                F(ii,1:7) = [0.0755 	25.2960 	0.3001 	11.5993 	0.6438 	4.0252 	-0.0196];
            case 'Eu3'
                F(ii,1:7) = [0.0204 	25.3078 	0.3010 	11.4744 	0.7005 	3.9420 	-0.0220];
            case 'Gd2'
                F(ii,1:7) = [0.0636 	25.3823 	0.3033 	11.2125 	0.6528 	3.7877 	-0.0199];
            case 'Gd3'
                F(ii,1:7) = [0.0186 	25.3867 	0.2895 	11.1421 	0.7135 	3.7520 	-0.0217];
            case 'Tb2'
                F(ii,1:7) = [0.0547 	25.5086 	0.3171 	10.5911 	0.6490 	3.5171 	-0.0212];
            case 'Tb3'
                F(ii,1:7) = [0.0177 	25.5095 	0.2921 	10.5769 	0.7133 	3.5122 	-0.0231];
            case 'Dy2'
                F(ii,1:7) = [0.1308 	18.3155 	0.3118 	7.6645 	0.5795 	3.1469 	-0.0226];
            case 'Dy3'
                F(ii,1:7) = [0.1157 	15.0732 	0.3270 	6.7991 	0.5821 	3.0202 	-0.0249];
            case 'Ho2'
                F(ii,1:7) = [0.0995 	18.1761 	0.3305 	7.8556 	0.5921 	2.9799 	-0.0230  ];
            case 'Ho3'
                F(ii,1:7) = [0.0566 	18.3176 	0.3365 	7.6880 	0.6317 	2.9427 	-0.0248];
            case 'Er2'
                F(ii,1:7) = [0.1122 	18.1223 	0.3462 	6.9106 	0.5649 	2.7614 	-0.0235 ];
            case 'Er3'
                F(ii,1:7) = [0.0586 	17.9802 	0.3540 	7.0964 	0.6126 	2.7482 	-0.0251];
            case 'Tm2'
                F(ii,1:7) = [0.0983 	18.3236 	0.3380 	6.9178 	0.5875 	2.6622 	-0.0241];
            case 'Tm3'
                F(ii,1:7) = [0.0581 	15.0922 	0.2787 	7.8015 	0.6854 	2.7931 	-0.0224];
            case 'Yb2'
                F(ii,1:7) = [0.0855 	18.5123 	0.2943 	7.3734 	0.6412 	2.6777 	-0.0213];
            case 'Yb3'
                F(ii,1:7) = [0.0416 	16.0949 	0.2849 	7.8341 	0.6961 	2.6725 	-0.0229];
            case 'Pr3'
                F(ii,1:7) = [0.0504 	24.9989 	0.2572 	12.0377 	0.7142 	5.0039 	-0.0219];
            otherwise
                warning('Selected Sys.Formfactor not recognised')
                disp('Format should be element then oxidation state, e.g.')
                disp('Sys.Formfactor = {Cr3}; %for a Cr3+ ion')
                disp('Currently only transition metals and lanthanides are recognised')
                disp('For other elements <j0> coefficients can be input directly as a 1x7 array')
                error('stopping here')
        end
        
    end
    
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%