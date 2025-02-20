function res=MahDis_James_accel(iFR,IdxL,lambda)
if nargin<3, lambda=0.5; end
% if ~iscell(IdxL{1}), IdxL{1}=IdxL; end;
res=[];
%matlabpool

% for l=1:length(lambda)
    % compute regularized Mahal distances between main sets
%     res{l}.dEuc{1}=[];
%     res{l}.dMah{1}=[];
%     res{l}.dMahR{1}=[];
    res.Euc=[];
    res.Mah=[];
%     res{l}.dMahRR{1}=[];
%     res{l}.CVErrR{1}=[];
%     res{l}.CVErrRR{1}=[];
    if length(IdxL)>1
        for i=1:length(IdxL)-1
            for j=i+1:length(IdxL)
%            for j=i+1
                X1=iFR(IdxL{1,i},:);
                X2=iFR(IdxL{1,j},:);
                [res.Euc(i,j),~, ...
                    res.Mah(i,j),~]=CompDisMs(X1,X2,lambda);
%                 [res{l}.CVErrR{1}(i,j),res{l}.CVErrRR{1}(i,j)]=CompCVErr(X1,X2,lambda(l));
            end;
        end;
    end;
    % compute regularized Mahal distances between sub-sets
%     for k=1:length(IdxL)
%         res{l}.dEuc{k+1}=[];
%         res{l}.dMah{k+1}=[];
%         res{l}.dMahR{k+1}=[];
%         res{l}.dMahRR{k+1}=[];
%         res{l}.CVErrR{k+1}=[];
%         res{l}.CVErrRR{k+1}=[];
%         if length(IdxL{k})>1
%             for i=1:length(IdxL{k})-1
%                 for j=i+1:length(IdxL{k})
%                     X1=iFR(:,unique(IdxL{k}{i}))';
%                     X2=iFR(:,unique(IdxL{k}{j}))';
%                     [res{l}.dEuc{k+1}(i,j),res{l}.dMah{k+1}(i,j), ...
%                         res{l}.dMahR{k+1}(i,j),res{l}.dMahRR{k+1}(i,j)]=CompDisMs(X1,X2,lambda(l));
%                     [res{l}.CVErrR{k+1}(i,j),res{l}.CVErrRR{k+1}(i,j)]=CompCVErr(X1,X2,lambda(l));
%                 end;
%             end;
%         end;
%     end;
% end;
%matlabpool close


function [dE,dM,dMR,dMRR]=CompDisMs(X1,X2,lambda)
[Cpool,CpoolR,CpoolRR,m1,m2,sel]=CovMtx(X1,X2,lambda);
dm=m1-m2;
dE=norm(dm);
dM=sqrt(dm(sel)*Cpool^-1*dm(sel)');
dMR=sqrt(dm(sel)*CpoolR^-1*dm(sel)');
dMRR=sqrt(dm*CpoolRR^-1*dm');


function [CVErrR,CVErrRR]=CompCVErr(X1,X2,lambda)
% leave-one-out prediction error
n1=size(X1,1); n2=size(X2,1);
yR=[]; yRR=[];
warning off
for ko=1:n1+n2
%parfor ko=1:n1+n2
    Xo1=X1(setdiff(1:n1,ko),:);
    Xo2=X2(setdiff(n1+1:n1+n2,ko)-n1,:);
    [Cpool,CpoolR,CpoolRR,m1,m2,sel]=CovMtx(Xo1,Xo2,lambda);
    if ko>n1, x=X2(ko-n1,:); else x=X1(ko,:); end;
    dR1=(x(sel)-m1(sel))*CpoolR^-1*(x(sel)-m1(sel))';
    dR2=(x(sel)-m2(sel))*CpoolR^-1*(x(sel)-m2(sel))';
    yR(ko)=sign(dR1-dR2);
    dRR1=(x-m1)*CpoolRR^-1*(x-m1)';
    dRR2=(x-m2)*CpoolRR^-1*(x-m2)';
    yRR(ko)=sign(dRR1-dRR2);
end;
warning on
CVErrR=(length(find(yR(1:n1)>0))/n1+ ...
    length(find(yR(n1+1:n1+n2)<0))/n2)/2;
CVErrRR=(length(find(yRR(1:n1)>0))/n1+ ...
    length(find(yRR(n1+1:n1+n2)<0))/n2)/2;


function [Cpool,CpoolR,CpoolRR,m1,m2,sel]=CovMtx(X1,X2,lambda)
n1=size(X1,1); n2=size(X2,1);
C1=cov(X1); C2=cov(X2);
CC=((n1-1)*C1+(n2-1)*C2)./(n1+n2-2);
% eliminate columns with zero variance
v=diag(CC); r=rank(CC);
sel=1:length(v);
if r<length(v)
    [vv,k]=sort(v,'descend');
    sel=k(1:r);
end;
Cpool=CC(sel,sel);
CpoolR=lambda*diag(diag(Cpool))+(1-lambda)*Cpool;
CpoolRR=lambda*mean(diag(CC))*eye(size(CC))+(1-lambda)*CC;
if n1>1, m1=mean(X1); else m1=X1; end;
if n2>1, m2=mean(X2); else m2=X2; end;


% (c) 2009 Daniel Durstewitz, CIMH & ICN, Univ. of Heidelberg, Germany
