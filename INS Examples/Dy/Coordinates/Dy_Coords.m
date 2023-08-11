function Dy_coords = Dy_Coords
%plot3(Dy(:,1),Dy(:,2),Dy(:,3),'.')
 %hold on

theta = deg2rad( 67 )  ;
phi = deg2rad(  3.5   );
psi = deg2rad(  19.4   );
Rx = [1 ,      0       ,           0;...
     0  ,  cos(theta)  , -sin(theta);...
     0  ,  sin(theta)  ,  cos(theta)];

Ry = [cos(phi)  ,  0  ,  sin(phi);...
      0         ,  1  ,  0       ;...
     -sin(phi)  ,  0  ,  cos(phi)];

Rz = [cos(psi)  ,  -sin(psi)  ,  0;...
      sin(psi)  ,  cos(psi)   ,  0;...
      0         ,     0       ,  1];

Dy_Coords_unaligned=load('Dy_Coords_unaligned.mat');

 Dy_coords = Dy_Coords_unaligned.Dy;
 Dy_coords = Dy_coords - Dy_coords(1,:);
 Dy_coords = Dy_coords*Rx*Rz*Ry;

% %  plot3(Dy1(:,1),Dy1(:,2),Dy1(:,3),'.')
% 
%%
% % 
%   top=[3,4,5,9];
% bottom=[2,6,7,8];
% close all
% hold on
% plot3(Dy_coords(1,1),Dy_coords(1,2),Dy_coords(1,3),'bs','MarkerSize',10,'MarkerFaceColor','b')
%   plot3(Dy_coords(top,1),Dy_coords(top,2),Dy_coords(top,3),'r.','MarkerSize',10)
%   plot3(Dy_coords(bottom,1),Dy_coords(bottom,2),Dy_coords(bottom,3),'r.','MarkerSize',10)
%   plot3(Dy_coords(2,1),Dy_coords(2,2),Dy_coords(2,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
%   plot3(Dy_coords(3,1),Dy_coords(3,2),Dy_coords(3,3),'ro','MarkerSize',10,'MarkerFaceColor','r')
% %title('Coordinates of Dy')
% %legend('Dy atom','','','Bottom O-N','Top O-N','location','E')
% lg{1} = plot(nan,'bs','MarkerSize',10,'MarkerFaceColor','b');
% lg{2} = plot(nan,'r.','MarkerSize',10);
% lg{3} = plot(nan,'ro','MarkerSize',10,'MarkerFaceColor','r');
% legend([lg{:}],{'Dy','O','O-N'},'location','E')
%   view(90,0)
%   hold off
% % 
% % 
% 
% %-79 - 13
% %-38 -30