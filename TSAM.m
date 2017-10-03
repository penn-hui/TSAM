%% 2017-8-21 this function is to implement the final prediction model TSAM
function Predict_score=TSAM(filename,pretype,featype,sgtype)
% filename--the .fa file containing the fasta format genome sequences
% pretype--prediction type: 1 for regression(default) cut human,mouse; 2 for
%        regression cut zebrafish; 3 for classification cut human,mouse
% featype--feature type: 2 for sequence feature only; 1 for sequence + cut
%           position feature(default); 3 for sequence feature+cut_per_geno
% sgtype--the type of sgRNAs: 1 for exons only; 0 for all including non-coding part
%% reading gene sequence: receive fasta files
load('Models.mat');
struct=Models{5,1};
basicE=Models{4,1};
if strfind(filename,'.fa') || strfind(filename,'.fasta')

    [seqinfo,pos_subseq,geneinfo]=MapGeno2Pro(filename,struct);
    %     if length(pos_subseq)==0 && length(geneinfo{6,2})==0
    %         disp('this sequence has no detail genome information');
    seqs_info{1,1}=seqinfo;
    seqs_info{1,2}=pos_subseq;
    seqs_info{1,3}=geneinfo;
    [N,a]=size(seqinfo);
    seqs=[];
    seq_final=cell(1,2);
    seq_final{1,1}=seqinfo{N,2};
    seq_final{1,2}=geneinfo{2,1};
    seqs=[seqs;seq_final];
    spacers=seaSpacer(seqs,seqs_info,4,3,'');

    [N,a]=size(spacers);
    test_seqs=spacers(:,1:3);
    cut_fea=zeros(N,3);
    for i=1:N
        cut_fea(i,1)=spacers{i,8};
        cut_fea(i,2)=spacers{i,11};
        cut_fea(i,3)=spacers{i,12};
    end
    if sgtype==1 % cut at coding region only
        index=find(cut_fea(:,3)>0);
        test_seqs=test_seqs(index,:);
        cut_fea=cut_fea(index,:);
    end
    save('test_seqs.mat','test_seqs','cut_fea');
    % call the xgboost
    system(['python xgboost_predict.py ',num2str(pretype),' ',num2str(featype)]);
    % load feautures and models
    test_fea_file=['test_fea',num2str(pretype),'_',num2str(featype),'.csv'];
    xgboost_file_inde_score=['xgboost_pre_score',num2str(pretype),'_',num2str(featype),'.csv'];
    xgboost_score=csvread(xgboost_file_inde_score);
    add_col=Models{pretype,1}{featype,1}{2,1};
    test_fea=csvread(test_fea_file);
    svm_model=Models{pretype,1}{featype,1}{1,1};
    pHMM1=Models{pretype,1}{featype,1}{3,1};
    pHMM2=Models{pretype,1}{featype,1}{4,1};
    
    % generate the pHMM feature
    au=1;
    ad=20;
    k=20;
    [n,a]=size(test_seqs);
    if pretype==3 % classification
        % encode sequence
        [tr_s,et_s,alpStat_s1]=hmmparak(spacers(:,1),1,au,ad,basicE);
        [tr_s,et_s,alpStat_s2]=hmmparak(spacers(:,1),2,au,ad,basicE);
        pHMM_fea=zeros(n,4);
        for i=1:n
            [STATES,LOPp1] = hmmdecode(alpStat_s1(i,:),pHMM1{1,1},pHMM1{1,2});
            [STATES,LOPn1] = hmmdecode(alpStat_s1(i,:),pHMM1{2,1},pHMM1{2,2});
            [STATES,LOPp2] = hmmdecode(alpStat_s2(i,:),pHMM2{1,1},pHMM2{1,2});
            [STATES,LOPn2] = hmmdecode(alpStat_s2(i,:),pHMM2{2,1},pHMM2{2,2});
            pHMM_fea(i,1)=LOPp1;
            pHMM_fea(i,2)=LOPn1;
            pHMM_fea(i,3)=LOPp2;
            pHMM_fea(i,4)=LOPn2;
        end
    elseif pretype<3 % regression
        [tr_tr,et_tr,alpStat_tr1]=hmmparak(spacers(:,1),1,au,ad,basicE);
        [tr_tr,et_tr,alpStat_tr2]=hmmparak(spacers(:,1),2,au,ad,basicE);
        test_fea1=zeros(n,k);
        for j=1:n
            %SgRNA=train_seq{j,1};
            for z=1:k
                [STATES,LOPp] = hmmdecode(alpStat_tr1(j,:),pHMM1{z,1},pHMM1{z,2});
                test_fea1(j,z)=LOPp;
            end
        end
        
        test_fea2=zeros(n,k);
        for j=1:n
            %SgRNA=train_seq{j,1};
            for z=1:k
                [STATES,LOPp] = hmmdecode(alpStat_tr2(j,:),pHMM2{z,1},pHMM2{z,2});
                test_fea2(j,z)=LOPp;
            end
        end
        pHMM_fea=[test_fea1,test_fea2];
    end
    test_data=[test_fea(:,add_col),pHMM_fea];
    test_score=zeros(n,1);
    
    if pretype==3
        [predict_lib,acc,dec_values] = svmpredict(test_score,test_data,svm_model,'-b 1');
        libsvm_score=dec_values(:,1);
        final_score=(libsvm_score+xgboost_score)/2;
        predict_score=zeros(n,1);
        predict_score(find(final_score>0.5),1)=1;
    elseif pretype<3
        [libsvm_score,acc,mse] = svmpredict(test_score,test_data,svm_model);
        final_score=(libsvm_score+xgboost_score)/2;
        predict_score=final_score;
    end
    Predict_score=test_seqs;
    for i=1:n
        Predict_score{i,4}=predict_score(i,1);
    end
else
    disp('please provide a .fa format sequence file');
    return;
    
end
%         disp('this file just containing genome information');

