function Qc=R_Qtidal(xin,Q,T,dH,Qc)
% critical river discharge
zz=xin;
[pks,locs]=findpeaks(zz);
IndMin=find(diff(sign(diff(zz)))>0)+1;
IndMax=find(diff(sign(diff(zz)))<0)+1;
ind=sort([IndMin;IndMax]);
dz=zz(ind(2:end))-zz(ind(1:end-1));
Tz=(T(ind(2:end))-T(ind(1:end-1)))*24;
Q1=Q(ind(1:end-1));T1=T(ind(1:end-1));
j=0;
for i=T(1):T(end)
    t1=i;t2=t1+3;
        j=j+1;
    m=find(T1>=t1&T1<t2);
    if length(m)==0
        dzmax(j)=0;
        qmax(j)=Qc;
        continue
    else
        dzmax(j)=max(dz(m));
        qmax(j)=max(Q1(m));
    end
end
m=find(qmax>mean(Q)*1.5 &dzmax<dH);
if length(m)>0
Qc=min(qmax(m));
end
