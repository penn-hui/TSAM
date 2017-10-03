%% this function is to search spacer sequences from given genome sequences
function spacers=seaSpacer(gSeq,ginfo,fn,bn,PAM)
% gSeq--genomic sequence
% ginfo--genomic informations
% PAM
% fn--forward extend number
% bn--backward extend number
[N,a]=size(gSeq);
if length(PAM)==0
    PAM='GG';
end

spacers=[];
for i=1:N
    % obtain its antisense sequences
    seq=gSeq(i,1);
    if length(seq{1,1})>0
        antigSeq=Spa2tar(seq);
        ind=strfind(seq,PAM);
        ind_a=strfind(antigSeq,PAM);
        n1=length(ind{1,1});
        n2=length(ind_a{1,1});
        spacer_sense=[];
        spacer_antisense=[];
        if length(ginfo)==0
            for j=1:n1
                spacer_info=cell(1,6);
                if ind{1,1}(1,j)-21-fn>=1 && ind{1,1}(1,j)+1+bn<=length(seq{1,1})
                    start_p=ind{1,1}(1,j)-21;
                    start_p_ex=ind{1,1}(1,j)-21-fn;
                    end_p=ind{1,1}(1,j)-2;
                    end_p_ex=ind{1,1}(1,j)+1+bn;
                    spacer=seq{1,1}(1,start_p:end_p);
                    spacer_exd=seq{1,1}(1,start_p_ex:end_p_ex);
                    spacer_info{1,1}=spacer;
                    spacer_info{1,2}=spacer_exd;
                    spacer_info{1,3}=gSeq{i,2};
                    spacer_info{1,4}=seq{1,1}(1,ind{1,1}(1,j)-1:ind{1,1}(1,j)+1);
                    spacer_info{1,5}=start_p;
                    spacer_info{1,6}=end_p;
                    spacer_sense=[spacer_sense;spacer_info];
                end
            end
            for j=1:n2
                spacer_info=cell(1,6);
                if ind_a{1,1}(1,j)-21-fn>=1 && ind_a{1,1}(1,j)+1+bn<=length(antigSeq{1,1})
                    start_p=ind_a{1,1}(1,j)-21;
                    start_p_ex=ind_a{1,1}(1,j)-21-fn;
                    end_p=ind_a{1,1}(1,j)-2;
                    end_p_ex=ind_a{1,1}(1,j)+1+bn;
                    spacer=antigSeq{1,1}(1,start_p:end_p);
                    spacer_ex=antigSeq{1,1}(1,start_p_ex:end_p_ex);
                    % spacer (N20NGG)
                    spacer_info{1,1}=spacer;
                    % N4N20NGGN3
                    spacer_info{1,2}=spacer_ex;
                    % gene name
                    spacer_info{1,3}=gSeq{i,2};
                    % PAM
                    spacer_info{1,4}=antigSeq{1,1}(1,ind_a{1,1}(1,j)-1:ind_a{1,1}(1,j)+1);
                    % start pos
                    spacer_info{1,5}=length(antigSeq{1,1})-end_p+1;
                    % end pos
                    spacer_info{1,6}=length(antigSeq{1,1})-start_p+1;
                    spacer_antisense=[spacer_antisense;spacer_info];
                end
            end
            spacers=[spacer_sense;spacer_antisense];
            spacers=sortrows(spacers,5);
        else
            pos_info=ginfo{i,2};
            pos_info=sortrows(pos_info,2);
            len_gene=ginfo{i,1}{find(strcmp(ginfo{i,1}(:,4),'GRC')),3};
            len_transcript=ginfo{i,1}{find(strcmp(ginfo{i,1}(:,4),'cdna')),3};
            len_protein=ginfo{i,1}{find(strcmp(ginfo{i,1}(:,4),'peptide')),3};
            if length(find(strcmp(ginfo{i,1}(:,4),'utr5')))==0
                len_utr5=0;
            else
                len_utr5=ginfo{i,1}{find(strcmp(ginfo{i,1}(:,4),'utr5')),3};
            end
            for j=1:n1
                spacer_info=cell(1,13);
                if ind{1,1}(1,j)-21-fn>=1 && ind{1,1}(1,j)+1+bn<=length(seq{1,1})
                    start_p=ind{1,1}(1,j)-21;
                    start_p_ex=ind{1,1}(1,j)-21-fn;
                    end_p=ind{1,1}(1,j)-2;
                    end_p_ex=ind{1,1}(1,j)+1+bn;
                    spacer=seq{1,1}(1,start_p:end_p);
                    spacer_exd=seq{1,1}(1,start_p_ex:end_p_ex);
                    spacer_info{1,1}=spacer;
                    spacer_info{1,2}=spacer_exd;
                    spacer_info{1,3}=gSeq{i,2};
                    spacer_info{1,4}=seq{1,1}(1,ind{1,1}(1,j)-1:ind{1,1}(1,j)+1);
                    spacer_info{1,5}=start_p;
                    spacer_info{1,6}=end_p;
                    % cut positions in genome
                    cut_g=ind{1,1}(1,j)-5;
                    cut_g_per=((cut_g-1)/length(seq{1,1}))*100;
                    spacer_info{1,7}=cut_g;
                    spacer_info{1,8}=cut_g_per;
                    % cut position in transcripts
                    exon_index_s=find(pos_info(:,2)<=cut_g);
                    exon_index_e=find(pos_info(:,3)>=cut_g);
                    %exon_index=find(pos_info(:,2)==max(exon_index_s) & pos_info(:,3)==min(exon_index_e));
                    exon_index=intersect(exon_index_s,exon_index_e);
                    if pos_info(exon_index,1)==0 % located at intron
                        cut_t=0;
                        cut_t_per=0;
                        % cut position in protein
                        cut_p=0;
                        cut_p_per=0;
                    else % located at exon
                        intros=find(pos_info(1:exon_index,1)==0);
                        if length(intros)==0 % cut position at the first exon
                            cut_t=cut_g;
                        else
                            cut_t=cut_g-sum(pos_info(intros,4));
                        end
                        cut_t_per=((cut_t-1)/len_transcript)*100;
                        cut_p=ceil((cut_t-len_utr5-1)/3);
                        cut_p_per=(cut_p/len_protein)*100;
                        if cut_p_per>100 || cut_p_per<0
                            % cut at the 3'utr or 5' utr
                            cut_p=0;
                            cut_p_per=0;
                        end
                    end
                    spacer_info{1,9}=cut_t;
                    spacer_info{1,10}=cut_t_per;
                    spacer_info{1,11}=cut_p;
                    spacer_info{1,12}=cut_p_per;
                    spacer_info{1,13}=1;
                    spacer_sense=[spacer_sense;spacer_info];
                end
                
            end
            
            for j=1:n2
                spacer_info=cell(1,13);
                if ind_a{1,1}(1,j)-21-fn>=1 && ind_a{1,1}(1,j)+1+bn<=length(antigSeq{1,1})
                    start_p=ind_a{1,1}(1,j)-21;
                    start_p_ex=ind_a{1,1}(1,j)-21-fn;
                    end_p=ind_a{1,1}(1,j)-2;
                    end_p_ex=ind_a{1,1}(1,j)+1+bn;
                    spacer=antigSeq{1,1}(1,start_p:end_p);
                    spacer_ex=antigSeq{1,1}(1,start_p_ex:end_p_ex);
                    % spacer (N20NGG)
                    spacer_info{1,1}=spacer;
                    % N4N20NGGN3
                    spacer_info{1,2}=spacer_ex;
                    % gene name
                    spacer_info{1,3}=gSeq{i,2};
                    % PAM
                    spacer_info{1,4}=antigSeq{1,1}(1,ind_a{1,1}(1,j)-1:ind_a{1,1}(1,j)+1);
                    % start pos
                    spacer_info{1,5}=length(antigSeq{1,1})-end_p+1;
                    % end pos
                    spacer_info{1,6}=length(antigSeq{1,1})-start_p+1;
                    % cut pos
                    cut_g=1+len_gene-ind_a{1,1}(1,j)+4;
                    spacer_info{1,7}=cut_g;
                    % cut percent
                    spacer_info{1,8}=((cut_g-1)/len_gene)*100;
                    % cut position in transcripts
                    exon_index_s=find(pos_info(:,2)<=cut_g);
                    exon_index_e=find(pos_info(:,3)>=cut_g);
                    exon_index=intersect(exon_index_s,exon_index_e);
                    %if length(cut_exon)
                    %exon_index=find(pos_info(:,2)==pos_info(max(exon_index_s),2) & pos_info(:,3)==pos_info(min(exon_index_e),3));
                    if pos_info(exon_index,1)==0 % located at intron
                        cut_t=0;
                        cut_t_per=0;
                        % cut position in protein
                        cut_p=0;
                        cut_p_per=0;
                    else % located at exon
                        intros=find(pos_info(1:exon_index,1)==0);
                        if length(intros)==0 % cut position at the first exon
                            cut_t=cut_g;
                        else
                            cut_t=cut_g-sum(pos_info(intros,4));
                        end
                        cut_t_per=((cut_t-1)/len_transcript)*100;
                        cut_p=ceil((cut_t-len_utr5-1)/3);
                        cut_p_per=(cut_p/len_protein)*100;
                        if cut_p_per>100 || cut_p_per<0
                            % cut at the 3'utr or 5' utr
                            cut_p=0;
                            cut_p_per=0;
                        end
                    end
                    %end
                    spacer_info{1,9}=cut_t;
                    spacer_info{1,10}=cut_t_per;
                    spacer_info{1,11}=cut_p;
                    spacer_info{1,12}=cut_p_per;
                    spacer_info{1,13}=0;
                    spacer_antisense=[spacer_antisense;spacer_info];
                end
                spacers=[spacer_sense;spacer_antisense];
                spacers=sortrows(spacers,7);
                
            end
        end
        
    end
end
