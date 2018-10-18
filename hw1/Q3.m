clc;
clear all;
%grid spacing
    h=1;
    %box dimensions in nm
    %choose in mutiples of cap dimensions
    a=500;
    b=500;

    nx=(a/h + 1);
    ny=(b/h + 1);
    v = zeros(nx,ny);
%cap dimensions in nm
    l1 = 20;
    l2 = 10;
    d = 10;
    %boundary conditions
    V2 = 1;
    V1 = 0;
    xedge_1b = a/2;
    nxedge_1b= xedge_1b/h+1;
    nxedge_1t= nxedge_1b-1;
    nxedge_21=  nxedge_1b + d/h;
    nxedge_22=  nxedge_21 + l2/h;
    yedge_11 = b/2;
    nyedge_11= yedge_11/h+1;
    nyedge_12 = nyedge_11+l1/h;
    nyedge_2b= nyedge_11;
    nyedge_2t = nyedge_2b+1;
    v(nxedge_1t:nxedge_1b,nyedge_11:nyedge_12)=V2;% postive plate
    v(nxedge_21:nxedge_22,nyedge_2b:nyedge_2t)=V1;%negative plate
    vnew = v;
    emax = 0;
    it=1;
    
    %solving using gauss-seidal iterative method 
    while true
        for i=2:nx-1
            for j=2:ny-1
                if ~(((i == nxedge_1b) || (i == nxedge_1t)) && (( j >= nyedge_11) ...
                        && ( j <= nyedge_12)) || ((j == nyedge_2b) || (j == nyedge_2t)) ...
                        && (( i >= nxedge_21)  && ( i <= nxedge_22)))
                    vnew(i,j) = (v(i-1,j)+v(i+1,j)+v(i,j-1)+v(i,j+1))/4;
                    e = abs((vnew(i,j) - v(i,j))/vnew(i,j));
                    if e > emax
                        emax = e;
                    end
                    v(i,j) = vnew(i,j);
                end
            end
        end
        disp(it);
        disp(emax);
        it=it+1;
        vt = vnew';
        if emax <= 0.1; break; end
        emax = 0;

    end
    
    %gradient
    for i=2:nx-1
        for j=2:ny-1
            if ( j == nyedge_2t)
                Ey(i,j) = -1*(vnew(i,j+1)-vnew(i,j));
            elseif (j == nyedge_2b)
                Ey(i,j) = -1*(vnew(i,j)-vnew(i,j-1));
            else
                Ey(i,j) = -1*(vnew(i,j+1)-vnew(i,j-1))/2;
            end
        end
    end
    for i=2:nx-1
        for j=2:ny-1
            if ( i == nxedge_1b)
                Ex(i,j) = -1*(vnew(i+1,j)-vnew(i,j));
            elseif (i == nxedge_1t)
                Ex(i,j) = -1*(vnew(i,j)-vnew(i-1,j));
            else
                Ex(i,j) = -1*(vnew(i+1,j)-vnew(i-1,j))/2;
            end
        end
    end
    
    %magnitude
    E = sqrt(Ex.^2 + Ey.^2);
    
    %applying Gauss law 
    %
    %find the net flux near the plate
    esum1 = sum(Ex(nxedge_1b,nyedge_11:nyedge_12));
    esum2 = -1*sum(Ex(nxedge_1t,nyedge_11:nyedge_12));
    
    %y components in the edge
    esum2e = -1*Ey(nxedge_1b,nyedge_11);
    esum1e = -1*Ey(nxedge_1t,nyedge_11);
    %etotal = esum1+esum2+esum3+esum4+esum5+esum6;
    etotal = esum1+esum2+esum1e+esum1e;
    emax = 0;
    %find emax
    for i=1:nx-1
        for j=1:ny-1
            if emax < E(i,j)
                emax = E(i,j);
                imax = i;
                jmax = j;
            end
        end
    end
    angle = (180/pi)*atan(Ey(imax,jmax)/Ex(imax,jmax));
    er = 1; % relative permitivity
    e0 = 8.85e-12; %permitivity of free space
    Dtotal = e0*etotal*er;
    Qtotal = Dtotal;
    cap = Qtotal/(V2-V1);
    [X,Y]=meshgrid(1:nx,1:ny);
    contour(X,Y,vnew',100)
    xlabel('x')
    ylabel('y')
    title('Equipotential lines')
    colorbar
    % hold on 
    figure
    [x,y]=meshgrid(1:nx-1,1:ny-1);
    quiver(x,y,Ex',Ey',2);
    xlabel('x')
    ylabel('y')
    title('Electric field profile')
    s1=sprintf('Cap : %d', cap);
    s2=sprintf('Emax : %f', emax);
    s3=sprintf('postition of Emax : %f %f', imax, jmax);
    s4=sprintf('angle w.r.t to xaxis in degree: %f', angle);
    disp(s1)
    disp(s2)
    disp(s3)
    disp(s4)

    