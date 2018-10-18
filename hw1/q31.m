 clc;
 clear all;
%grid spacing
    h=1;
    %box dimensions in nm
    %choose in mutiples of cap dimensions
    a=2000;
    b=500;

    nx=(a/h + 1);
    ny=(b/h + 1);
    v = zeros(nx,ny);
%     cap dimensions in nm
    l = 1000;
    d = 10;
    %boundary conditions
    xedge_1 = (a-l)/2;
    nxedge_1= xedge_1/h+1;
    xedge_2 = (a+l)/2;
    nxedge_2= xedge_2/h+1;
    yedge_1 = (b-d)/2;
    nyedge_1t= yedge_1/h+1;
    nyedge_1b = nyedge_1t-1;
    yedge_2 = (b+d)/2;
    nyedge_2b= yedge_2/h+1;
    nyedge_2t = nyedge_2b+1;
    V2 = 1;
    V1 = 0;
    v(nxedge_1:nxedge_2,nyedge_2t)=V2;% postive plate
    v(nxedge_1:nxedge_2,nyedge_2b)=V2;
    v(nxedge_1:nxedge_2,nyedge_1t)=V1;%negative plate
    v(nxedge_1:nxedge_2,nyedge_1b)=V1;
    vnew = v;
    emax = 0;
    it=1;
    
    %solving using gauss-seidal iterative method 
    while true
        for i=2:nx-1
            for j=2:ny-1
                if ~(((i >= nxedge_1) && (i <= nxedge_2)) && ((j==nyedge_1t) || (j==nyedge_1b) || (j==nyedge_2t) || (j==nyedge_2b)))
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
            Ex(i,j) = -1*(vnew(i+1,j)-vnew(i-1,j))/2;
        end
    end
    for i=2:nx-1
        for j=2:ny-1
            if (((i >= nxedge_1) && (i <= nxedge_2)) && ((j==nyedge_1t) || (j==nyedge_2t) ))
                Ey(i,j) = -1*(vnew(i,j+1)-vnew(i,j));
            elseif (((i >= nxedge_1) && (i <= nxedge_2)) && ((j==nyedge_1b) || (j==nyedge_2b) ))
                Ey(i,j) = -1*(vnew(i,j)-vnew(i,j-1));
            else
                Ey(i,j) = -1*(vnew(i,j+1)-vnew(i,j-1))/2;
            end
        end
    end
    [ex,ey]=gradient(vnew');
    ex = -ex;
    ey = -ey;
    %applying Gauss law 
    %
    %find the net flux near the plate
   esum1 = sum(Ey(nxedge_1:nxedge_2,nyedge_2t))
   esum2 = -1*sum(Ey(nxedge_1:nxedge_2,nyedge_2b))
   esum1e = sum(ey(nyedge_2t,nxedge_1:nxedge_2));
   esum2e = -1*sum(ey(nyedge_2b,nxedge_1:nxedge_2))
   
   esum3 =  sum(Ey(nxedge_1:nxedge_2,nyedge_2b))
   esum4 = -1*sum(Ey(nxedge_1:nxedge_2,nyedge_2t))
   
   esum3e = sum(ey(nyedge_2b,nxedge_1:nxedge_2))
   esum4e =  -1*sum(ey(nyedge_2t,nxedge_1:nxedge_2))
      
%     esum3 = -1*Ex(nxedge_1,nyedge_2t);
%     esum4 = -1*Ex(nxedge_1,nyedge_2b);
%     esum5 = Ex(nxedge_2,nyedge_2t);
%     esum6 = Ex(nxedge_2,nyedge_2b);

    %etotal = esum1+esum2+esum3+esum4+esum5+esum6;
    etotal = esum2+esum1;
    etotale = esum2e+esum1e;
    etotal1e = esum3e+esum4e;
    etotal1 = esum3+esum4;
    er = 1; % relative permitivity
    e0 = 8.85e-12; %permitivity of free space
    Dtotal = e0*etotal*er;
    Dtotale = e0*etotale*er;
    Dtotal1e = e0*etotal1e*er;
    Dtotal1 = e0*etotal1*er;
    Qtotal = Dtotal;
    Qtotal1 = Dtotal1;
    cap = Qtotal/(V2-V1);
    cap_ideal = (l*e0*er)/d;
    par_cap = cap - cap_ideal;
    Qtotal - Qtotal1
    % [X,Y]=meshgrid(1:nx,1:ny);
    % contour(X,Y,vnew',31)
    % hold on 
    % [x,y]=meshgrid(1:nx-1,1:ny-1);
    % quiver(x,y,Ex',Ey');
    % hold on
    % startx = 400:1600;
    % streamline(x,y,Ex',Ey',startx);
