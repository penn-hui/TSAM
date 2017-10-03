%% 2017-3.22 transfer spacer to its target sequence
function tseqs=Spa2tar(seqs)
[N,a]=size(seqs);
tseqs=cell(N,1);
for i=1:N
    seq=seqs{i,1};
    l=length(seq);
    tseq=[];
    for j=1:l
        ind=l+1-j;
        %ind=j;
        if strcmp(seq(1,ind),'A')
            tseq=[tseq,'T'];
        elseif strcmp(seq(1,ind),'G')
            tseq=[tseq,'C'];
        elseif strcmp(seq(1,ind),'C')
            tseq=[tseq,'G'];
        elseif strcmp(seq(1,ind),'T')
            tseq=[tseq,'A'];
        elseif strcmp(seq(1,ind),'N')
            tseq=[tseq,'N'];
        end
    end
    tseqs{i,1}=tseq;
end