% Exp.SpectrumType = 'SE';Exp.Temperature = 24;
% Exp.Energy = 0:0.05:4; Exp.lwfwhm = 0.4; Exp.Q = 0.2:0.05:1;



% 



SysMint.FormFactor = 'Cr3';
a = [0;-3.3807;-5.7290;-5.5960;-3.1148;0.2482];    
b = [0;0.0840;-2.3170;-5.6470;-7.8270;-7.7410];
c = [0; -0.2765 ; 0.1103 ; -0.2956 ; 0.2920 ; -0.1803];
SysMint.Coords = [a b c];
Opt.NumEigs = 24;
ExpMint.SpectrumType = 'SE'; 
ExpMint.Temperature = [1.5, 6, 15]; 
ExpMint.Energy = -1:0.001:4; 
ExpMint.lwfwhm = 0.05;


ExpMint.Q = [0.1,0.5,1,1.5,2]; [cross,Eigs] = mint(SysMint,ExpMint,Opt);
figure
plot(ExpMint.Energy,cross)
xlim([0 3.25])
xlabel('Energy (meV)') 
ylabel('Intensity (arb. units)') 
legend('1.5 K','6 K','15 K')


% 
% 
% [cross,Eigs] = mint(Sys1,Exp);
% 
% figure
% plot(Exp.Energy,cross)