%% 2017-7-20 this function is to map the cut position to the protein
function [seqinfo,pos_subseq,geneinfo]=MapGeno2Pro(fafile,struct)

%fafile--the fasta file of a given gene with all the exon, intro,
%chromosome info
% name--gene name
% position--the cut position

%[seqs,leninfo]=ReadFastaFiles(pathname,key_file,key_head)
path_all=fafile;
seq_all = fastaread(path_all);
%struct=['exon';'intron';'utr3';'utr5';'GRC','cds','cdna'];
N=length(seq_all);
%seq_final=cell(N,1);
seqinfo=cell(N,6);
exon_pos_info=[];


for j=1:N
    seq=seq_all(j);
    seq_name=seq.Header;
    seqinfo{j,1}=seq_name;
    seqinfo{j,2}=seq.Sequence;
    seqinfo{j,3}=length(seq.Sequence);
    % infer the type of the subsequence
    for i=1:8
        pos=strfind(seq_name,struct(i,1));
        if length(pos)==1
            seqinfo{j,4}=struct{i,1};
            if i<3
                exon_pos_info(j,1)=1;
            else
                exon_pos_info(j,1)=0;
            end
        end
        
    end
end

if length(exon_pos_info)==0 && length(find(strcmp(seqinfo,'GRC')))==0
    pos_subseq=[];
    geneinfo=[];
else
    subseq=seqinfo(find(exon_pos_info==1),:);
    chr_seq=seqinfo{find(strcmp(seqinfo(:,4),'GRC')),2};
    
    [L,a]=size(subseq);
    pos_subseq=zeros(L,4);
    for i=1:L
        pos=strfind(chr_seq,subseq{i,2});
        if length(pos)==1
            if strcmp(subseq{i,4},'exon')
                pos_subseq(i,1)=1;
            elseif strcmp(subseq{i,4},'intron')
                pos_subseq(i,1)=0;
            end
            pos_subseq(i,2)=pos(1,1);
            pos_subseq(i,3)=pos(1,1)+subseq{i,3}-1;
            pos_subseq(i,4)=subseq{i,3};
            index=find(strcmp(seqinfo(:,1),subseq(i,1)));
            seqinfo{index,5}=pos(1,1);
            seqinfo{index,6}=pos(1,1)+subseq{i,3}-1;
        end
    end
    
    gene_name_in=strfind(seqinfo{1,1},' ');
    geneinfo=cell(2,6);
    geneinfo{1,1}='gene_name';
    geneinfo{1,2}='chromosome_name';
    geneinfo{1,3}='gene_start';
    geneinfo{1,4}='gene_end';
    geneinfo{1,5}='strand';
    geneinfo{1,6}='genome version';
    % genome=e1+i1+e2+i2+...
    % cDNA=e1+e2+...
    % cds=cDNA-utr5-utr3;
    % transcription start site=exon1 pos1
    % geneinfo{1,7}='gene_TSS';
    % % translation start site=utr5 pos end+1
    % geneinfo{1,8}='gene_tss';
    geneinfo{2,1}=seqinfo{1,1}(1,1:gene_name_in-1);
    chromosome_name=seqinfo{N,1};
    seq_se_in=strfind(chromosome_name,':');
    seq_chr_in=strfind(chromosome_name,' ');
    seq_chr=chromosome_name(1,1:seq_chr_in(1,1)-1);
    geneinfo{2,2}=seq_chr;
    seq_start=str2num(chromosome_name(1,seq_se_in(length(seq_se_in)-2)+1:seq_se_in(length(seq_se_in)-1)-1));
    seq_end=str2num(chromosome_name(1,seq_se_in(length(seq_se_in)-1)+1:seq_se_in(length(seq_se_in))-1));
    geneinfo{2,3}=seq_start;
    geneinfo{2,4}=seq_end;
    geneinfo{2,5}=str2num(chromosome_name(1,seq_se_in(length(seq_se_in))+1:end));
    geneinfo{2,6}=chromosome_name(1,seq_se_in(length(seq_se_in)-4)+1:seq_se_in(length(seq_se_in)-3)-1);
end



