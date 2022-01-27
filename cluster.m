function ain=cluster(ain,clusang);
% CLUSTER: Clusters angles in rows around the angles in the first 
% column. CLUSANG is the allowable ambiguity (usually 360 degrees but
% sometimes 180).

ii=(ain-ain(:,ones(1,size(ain,2))))>clusang/2;
ain(ii)=ain(ii)-clusang;
ii=(ain-ain(:,ones(1,size(ain,2))))<-clusang/2;
ain(ii)=ain(ii)+clusang;


%----------------------------------------------------------------------