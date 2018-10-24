%===========================================================
%             Basic LB Code for Poiseuille Flow
%===========================================================
%        developed by Prof. Giacomo Falcucci, PhD
%          for AC274 Class @ HARVARD UNIVERSITY
%                          2017
%___________________________________________________________
%
%vars
clear all
lx=100;
ly=40;
w=[4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36];
ex=[0,1,0,-1,0,1,-1,-1,1];
ey=[0,0,1,0,-1,1,1,-1,-1];
tau=1;
g=10^-4;
kr = 1;
cs=1/sqrt(3);
tstep=5000;
tplot=100;



% arrays
rho=ones(lx,ly);
u=zeros(lx,ly);
v=zeros(lx,ly);
storeu = zeros(1000,lx,ly);
storev = zeros(1000,lx,ly);
f_in=zeros(9,lx,ly);
opp=[1,4,5,2,3,8,9,6,7];
pshape = zeros(1000,9,lx,ly);

%% initialization

for kk=1:9
   f_in(kk,:,:)=w(kk)*rho; 
end

f_out=f_in;

% main loop

for tt=1:tstep
    rest = cputime;
    %% macroscopic
    
    for ii=1:lx
       for jj=1:ly
           rho(ii,jj)=0;
           u(ii,jj)=0;
           v(ii,jj)=0;
           
           for kk=1:9
                rho(ii,jj)=rho(ii,jj)+f_in(kk,ii,jj);
                u(ii,jj)=u(ii,jj)+f_in(kk,ii,jj)*ex(kk);
                v(ii,jj)=v(ii,jj)+f_in(kk,ii,jj)*ey(kk);
           end
           u(ii,jj)=u(ii,jj)/rho(ii,jj);
           v(ii,jj)=v(ii,jj)/rho(ii,jj);
       end
    end

%% collision
for ii=1:lx
   for jj=1:ly
       u(ii,jj)=u(ii,jj) + g*tau/rho(ii,jj);
       for kk=1:9
          cu=(1/cs^2)*(ex(kk)*u(ii,jj) + ey(kk)*v(ii,jj));
          feq=rho(ii,jj)*w(kk)*(1 + cu + cu^2 - (u(ii,jj).^2 + v(ii,jj).^2)/(2*cs^2));
          f_in(kk,ii,jj)=f_in(kk,ii,jj) + (1/tau)*(feq - f_in(kk,ii,jj));
       end
        
   end
end

%% prepare the next time step

  fout=f_in; 
%%bbboundary

%north & south
for ii=1:lx
    for kk=1:9
      fout(kk,ii,ly)=f_in(opp(kk),ii,ly);
      fout(kk,ii,1)=f_in(opp(kk),ii,1);
    end
end
%% streaming

for kk=1:9
   f_in(kk,:,:)=circshift(fout(kk,:,:),[0,ex(kk),ey(kk)]);
end

if (mod(tt,100) == 1)
    pshape(kr,:,:,:) = f_in;
    storeu(kr,:,:) = u;
    storev(kr,:,:) = v;
    kr = kr + 1;

end
display(cputime-rest)
if (mod(tt,tplot)==1)
   imagesc(u') 
   axis xy
   colorbar
   drawnow
end
end
