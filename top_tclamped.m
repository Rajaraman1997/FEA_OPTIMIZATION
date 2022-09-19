
%%% Run: top_tclamped(280,40,0.5,3,0.5) clamped clamped beam;

function top_tclamped(nelx,nely,volfrac,penal,rmin)
% INITIALIZE
x(1:nely,1:nelx) = volfrac;
loop = 0;
change = 1.;
lambda1=1;
% START ITERATION
while change > 0.005 
    loop = loop + 1;
    xold = x;
    % FE-ANALYSIS
    [lambda1,Un]=FE(nelx,nely,x,penal);
    % OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    [KE] = lk;
    [ME]=lm;
    beta=lambda1;
    for ely = 1:nely
        for elx = 1:nelx
            n1 = (nely+1)*(elx-1)+ely;
            n2 = (nely+1)* elx   +ely;
            Une = Un([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],1);
            if x(ely,elx)<0.1
                c1=6e5;c2=-5e6;
               dc(ely,elx) = -(Une')*((penal*x(ely,elx)^(penal-1)*KE)-(lambda1*((6*c1*x(ely,elx)^5)+(7*c2*x(ely,elx)^6))*ME))*(Une);
            else
               dc(ely,elx) = -(Une')*((penal*x(ely,elx)^(penal-1)*KE)-lambda1*ME)*(Une);
            end
        end
    end
    % FILTERING OF SENSITIVITIES
    [dc]   = check(nelx,nely,rmin,x,dc);
    % DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    [x]    = OC(nelx,nely,x,volfrac,dc);
    % PRINT RESULTS
    change = max(max(abs(x-xold)));
    disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4f',beta) ...
        ' Vol.: ' sprintf('%6.3f',sum(sum(x))/(nelx*nely)) ...
        ' ch.: ' sprintf('%6.3f',change )])
    % PLOT DENSITIES
    colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-4);
end
%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xnew]=OC(nelx,nely,x,volfrac,dc)
l1 = 0; l2 = 100000; move = 0.40;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*max(0,(-dc./lmid)).^0.5))));
    if sum(sum(xnew)) - volfrac*nelx*nely > 0
        l1 = lmid;
    else
        l2 = lmid;
    end
end
%%%%%%%%%% MESH-INDEPENDENCY FILTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dcn]=check(nelx,nely,rmin,x,dc)
dcn=zeros(nely,nelx);
for i = 1:nelx
    for j = 1:nely
        sum=0.0;
        for k = max(i-floor(rmin),1):min(i+floor(rmin),nelx)
            for l = max(j-floor(rmin),1):min(j+floor(rmin),nely)
                fac = rmin-sqrt((i-k)^2+(j-l)^2);
                sum = sum+max(0,fac);
                dcn(j,i) = dcn(j,i) + max(0,fac)*x(l,k)*dc(l,k);
            end
        end
        dcn(j,i) = dcn(j,i)/(x(j,i)*sum);
    end
end
%%%%%%%%%% FE-ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lambda1,Un]=FE(nelx,nely,x,penal)
[KE] = lk;
[ME]=  lm;
K = sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
M=  sparse(2*(nelx+1)*(nely+1), 2*(nelx+1)*(nely+1));
for elx = 1:nelx
    for ely = 1:nely
        n1 = (nely+1)*(elx-1)+ely;
        n2 = (nely+1)* elx   +ely;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof,edof) = K(edof,edof) + x(ely,elx)^penal*KE;
        if x(ely,elx)>0.1
            M(edof,edof) = M(edof,edof) + x(ely,elx)*ME;
        else
            c1=6e5;c2=-5e6;
            M(edof,edof) = M(edof,edof) + ((c1*x(ely,elx)^6)+(c2*x(ely,elx)^7))*ME;
        end
    end
end
M((nely+1)*(nelx+1),(nely+1)*(nelx+1))=M((nely+1)*(nelx+1),(nely+1)*(nelx+1))+125;
%%%% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
%fixeddofs   = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);

%%%% For clamped clamped beam

fixeddofs = union([1:2*(nely+1)],[(2*(nelx+1)*(nely+1)-2*(nely)-1):(2*(nelx+1)*(nely+1))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alldofs     = [1:2*(nely+1)*(nelx+1)];
gdof=2*(nely+1)*(nelx+1);
freedofs    = setdiff(alldofs,fixeddofs);
% SOLVING
modes = zeros(gdof,1);
[V,D] = eigs(K(freedofs,freedofs),M(freedofs,freedofs),1,'smallestabs');
eigenvalues = diag(D);
modes(freedofs,:) = V;
lambda1=eigenvalues;
unorm=norm(modes);
U=modes./unorm;
wn= sqrt(eigenvalues)
norm1=sqrt(((U)')*(M)*(U));
Un=U./(norm1);
%%%%%%%%%% ELEMENT STIFFNESS MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [KE]=lk
E = 25e9;
nu = 0.3;
k=[ 1/2-nu/6   1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 ...
    -1/4+nu/12 -1/8-nu/8  nu/6       1/8-3*nu/8];
KE = E/(1-nu^2)*[ k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];
function [ME]=lm
rho=2500;
m=rho*14*2*0.01;
ME=(m/72).*            [  4 0 2 0 1 0 2 0;
    0 4 0 2 0 1 0 2;
    2 0 4 0 2 0 1 0;
    0 2 0 4 0 2 0 1;
    1 0 2 0 4 0 2 0;
    0 1 0 2 0 4 0 2;
    2 0 1 0 2 0 4 0;
    0 2 0 1 0 2 0 4];
%       ME=(m/32)* [1 0 1 0 1 0 1 0;
%                   0 1 0 1 0 1 0 1;
%                   1 0 1 0 1 0 1 0;
%                   0 1 0 1 0 1 0 1;
%                   1 0 1 0 1 0 1 0;
%                   0 1 0 1 0 1 0 1;
%                   1 0 1 0 1 0 1 0;
%                   0 1 0 1 0 1 0 1];






