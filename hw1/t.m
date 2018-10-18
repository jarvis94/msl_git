%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written for Course :- Computational Electromagnetics, Fall 2011
%                       Department of Electrical Engineering
%                       Indian Institute of Technology Madras
%                       Chennai - 600036, India
%
% Authors            :- Sathya Swaroop Ganta, B.Tech., M.Tech. Electrical Engg.
%                       Kayatri, M.S. Engineering Design
%                       Pankaj, M.S. Electrical Engg.
%                       Sumanthra Chaudhuri, M.S. Electrical Engg.
%                       Projesh Basu, M.S. Electrical Engg.
%                       Nikhil Kumar CS, M.S. Electrical Engg.
%
% Instructor :- Ananth Krishnan
%               Assistant Professor
%               Department of Electrical Engineering
%               Indian Institute of Technology Madras
%
% Any correspondance regarding this program may be addressed to
% Prof. Ananth Krishnan at 'computational.em.at.iit.madras@gmail.com'
%
% Copyright/Licensing :- For educational and research purposes only. No
% part of this program may be used for any financial benefit of any kind
% without the consent of the instructor of the course.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Finite Difference Method (FDM) solution to Laplacian"
% 
% Objective of the program is to solve for the steady state voltage
% distribution in a region 0<x<30, 0<y<30, given that one of the sides of
% square is excited with a voltage of 45*(x)*(1-x) Volts and all other
% sides are maintained at 0 Volts. This voltage at the boundary is 
% symmetrical with its maximum value at centre of the boundary namely x=15.
% At any iteration, the value of voltage is updated as average of voltages
% of 4 nearest naighbors, until between consecutive iterations, the error
% is less than 0.01 V.
%
% The tolerance in error between iterations is kept at 0.01 V. This may be
% tweaked to a higher or lower value for lower or higher accuracy 
% respectively. Imagesc command by default uses image axis settings, which
% are different from normal plot command and hence x and y axis may look 
% flipped. Read Matlab documentation on imagesc for more details.
%
% Program stops when iteration number in plot does not change or can be
% closed anytime by just closing the plot window. On normal completion, the
% program plots the electric field in a quiver plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clearing variables in memory and Matlab command screen
clear all;
clc;

%Dimensions of the simulation grid in x (xdim) and y (ydim) directions
xdim=30;
ydim=15;

%Initializing previous (V_prev) and present (V_now) voltage matrices
V_now=zeros(xdim+1,ydim+1);
V_prev=zeros(xdim+1,ydim+1);

%Initializing boundary conditions only for V_now
i=1:1:xdim+1;%x-co-ordinates for boundary at y=ydim*grid_size 

%A voltage of 45*x*(1-x) is applied on one boundary, the remaining
%boundaries are going to remain at zero volts
V_now(i,ydim+1)=45*((i-1)/xdim).*(1-((i-1)/xdim));

%Iteration counter
iter=0;

%Calculation of maximum error between V_now and V_prev at all points
%By setting the applied voltage for only V_now, we have made V_no and
%V_prev different, hence error will be greater than zero and the program
%will enter the while loop following this command.
error=max(max(abs(V_now-V_prev)));

%Iteration loop
while(error>0.01)%Run this until convergence
    
    iter=iter+1; % Iteration counter increment
    
    % Updating present iteration using 4 point Central diffrence form
    % of Laplace equation obtained using Finite Difference method
    for i=2:1:xdim
        for j=2:1:ydim
            V_now(i,j)=(V_now(i-1,j)+V_now(i+1,j)+V_now(i,j-1)+V_now(i,j+1))/4;
        end
    end
    error=max(max(abs(V_now-V_prev))); % Calculate the maximum error between previous and current iteration at all points
    V_prev=V_now; % Updating previous iteration matrix to the last iteration performed
    
    %Movie type colour scaled image plot to see how solution progresses
    imagesc(V_now);colorbar;
    title(['Voltage distribution on a ',int2str(xdim),' x ',int2str(ydim),' grid at iteration no ',int2str(iter)],'Color','k'); 
    getframe;
end

%Plot the electric field distribution
figure;
[ex,ey]=gradient(V_now);
quiver(-ex,-ey); %Quiver command creates a plot, E=-grad(V), hence the negative sign

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF PROGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%