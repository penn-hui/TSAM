%% this function is to extend the function hmmpara with k-mer emission
function [tr,et,alpStat]=hmmparak(Seq,t,au,ad,basicE)
% seq--the target sequence
% t--the k-mer,k=1 or k=2

[N,a]=size(Seq);
l=length(Seq{1,1});
Fa=FullAarry(basicE,t);
[L,a]=size(Fa);
alpStat=zeros(N,l+1-t);
for i=1:N
    seq=Seq{i,1};
    for j=1:l+1-t
        mer=seq(1,j:j+t-1);
        index=find(strcmp(Fa,mer));
        alpStat(i,j)=index;
    end
end
b=ones(1,l+2-t);
c=zeros(l+2-t,1);
tr=[c,diag(b)];
tr(l+3-t,1)=0;
et=zeros(l+3-t,L);
et(1,:)=0;
et(l+3-t,:)=0;
for i=2:l+2-t
    for j=1:L
%     p=length(find(alpStat(:,i-1)==j))/N;
    p=(length(find(alpStat(:,i-1)==j))+au)/(N+ad);
    et(i,j)=p;
    end
end
