function [xp,ymin]=fmin_Q_tide(z,q,t,fu,fband,f,v,Qc)
    z1=z;
    m=find(~isnan(z1));
    zz=z1(m);qq=q(m);n=length(zz);
    A=[];b=[];Aeq=[];beq=[];nonlcon=[];
    options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');
    x0(1:2) =[1.2, 1];x0(3:11)=1.1;
    lb(1:2)=[0.2,0];ub(1:2)=[2,24];lb(3:11)=0;ub(3:11)=2;
        fun = @(x)Qztidefun(qq,zz,t,fu,fband,x,f,v,Qc);
        [xp(:,1),ymin]= fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

