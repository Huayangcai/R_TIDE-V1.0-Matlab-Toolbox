function [cof,per,fband]=Rtide_pso(z,q,t,fu,id,fuv,ff,Qc,fband)
n=length(fu);% NUMBER OF tidal constituents 
[M,N]=size(z);%length of data and N stations
%optimization for p and taoQ,pr(n),PSD optimization based on the least square method
NN = 40;                         % initial swarm size  
ger = 100;                       % maximum iteration number       
limit = [0, 100];                % minimum and maximum value of every particle x(:,1:d)  
vlimit = [-0.8, 0.8];           % minimum and maximum value of every particle velocity  
w = 0.8;                        % inertia weight  
wmax=0.9;wmin=0.4;              % the maximum and minimum of inertia weight 
c1 = 0.5;                       % particle learning rate  
c2 = 0.7;                       % swarm learning rate 
ii=1:NN;   
% apart from p,tao related to discharge stage£¬tidal constituent bonds are adopted to derive optimized parameters£¬and d
d=0;j=0;
 for i=1:length(fband)
       m=find(fu>fband(i,1) & fu<fband(i,2));
       if isempty(m)
           j=j+1;
           i1(j)=i;
       else
           d=d+1;
       end
 end
 d = d+2;  % particle (variable x) size 
 if j>0
 fband(i1,:)=[];
 end

for k=1:N   %N station
    %initial swarm position x and velocity v for swarm size NN and 5 particle
    k
     x(ii,1:d)=1.3;x(ii,2)=ii;
     v(ii,1:d)=vlimit(1) + (vlimit(2) - vlimit(1)) * rand(NN,d);

  xm = x;                          % initial history best position of every particle  
  ym = zeros(1, d);                % initial history best position of swarm  
  if id==1
    fxm = zeros(NN, 1);              % initial history best fitness of every particle 
    fym = -inf;                      % initial history best fitness of swarm 
  else
    fxm = ones(NN, 1)*10;              % initial history best fitness of every particle 
    fym = inf;                      % initial history best fitness of swarm 
  end

  iter = 1;    %initial iteration
  record = zeros(ger, 1);            
 while iter <= ger  
      fx = Rtide_sixcof(x,q,z(:,k),t,fu,fband,id,fuv,ff,Qc);  % % local fitness of all particles in swarm,object function R2   
     for i = 1:NN
         if id==1
            if fxm(i) < fx(i)  
                fxm(i) = fx(i);     % updated history best finess of every particle  
                xm(i,:) = x(i,:);   % updated history best position of every particle  
            end   
         else
           if fxm(i) > fx(i)  
               fxm(i) = fx(i);     % updated history best finess of every particle  
               xm(i,:) = x(i,:);   % updated history best position of every particle  
           end   
         end
     end  
 if id==1
    if fym < max(fxm)  
        [fym, nmax] = max(fxm);   % updated history best fitness of swarm 
        ym = xm(nmax, :);      % updated history best position of swarm   
    end
 else
    if fym > min(fxm)  
        [fym, nmin] = min(fxm);   % updated history best fitness of swarm 
        ym = xm(nmin, :);      % updated history best position of swarm   
    end  
 end
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, NN, 1) - x);% update velocity of all particle in swarm   
  % limit the velocity of boundary ,this choice may lead the result change
    %v(v > vlimit(2)) = vlimit(2);  
    %v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% updated position of all particle in swarm  
    % limit the position of boundary  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;  
    iter = iter+1;  
 end 
 cof(k,:)=ym;  %define best position of swarm
 per(k)=fym;   %define best fitness(maximum R2) of swarm
end
