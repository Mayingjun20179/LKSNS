function Comparison_new_methods_integrated_dataset
clc
clear 
result = five_cross();
save result result;
end


%%%%%Comparison with relatively new methods on integrated dataset
function result = five_cross()
seed = 1;
load extracted_interaction.txt;
load protein_ctd;
load extracted_lncRNA_sequence_CT.txt
load extracted_lncRNA_expression.txt
interaction_matrix = extracted_interaction;
[NR,NP] = size(extracted_interaction);
CV=5;
rand('state',seed);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);     
result = zeros(5,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    %%测试数据集和训练数据集的选取
    if (k~=CV)
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    train_set=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_set(test_row_index(i),test_col_index(i))=0;                 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LKSNS
    %proteins’ interaction profile   
    canshu.neighbor_num = floor(0.3*NP);    canshu.sigma = 0.5;  canshu.knum = 9; 
    canshu.flag = 1;   canshu.jingdu = 0.00001;    
    SPC = KSNS_opt(train_set',canshu);      
    RPC = calculate_labels(SPC ,train_set',0.7)'; %%%标签传递  
    %proteins’ CTD
    canshu.neighbor_num = floor(0.5*NP);    canshu.sigma = 0.3;  canshu.knum = 17;  
    canshu.flag = 1;    canshu.jingdu = 0.00001;
    SPX = KSNS_opt(protein_ctd,canshu);     
    RPX = calculate_labels(SPX,train_set',0.1)'; %%%标签传递  
    %lncRNAs’ interaction profile 
    canshu.neighbor_num = floor(0.1*NR);    canshu.sigma = 0.5;  canshu.knum=900;         
    canshu.jingdu = 0.0001;  canshu.flag= 1;
    SRC = KSNS_opt(train_set,canshu);      
    RRC = calculate_labels(SRC,train_set,0.5);
    %lncRNA's expression profile   
    canshu.neighbor_num = floor(0.1*NR);    canshu.a = 0.9;    canshu.b = 2;
    canshu.flag= 2 ;     canshu.jingdu =  0.0001;
    SRP =  KSNS_opt(extracted_lncRNA_expression,canshu);
    RRP = calculate_labels(SRP,train_set,0.3);    
    %lncRNAs’ Sequence composition
    canshu.neighbor_num = floor(0.7*NR);    canshu.flag=2;  canshu.a = 0.1; canshu.b = 4;
    canshu.jingdu = 0.0001;
    SRX = KSNS_opt(extracted_lncRNA_sequence_CT,canshu); %%%计算权重矩阵
    RRX = calculate_labels(SRX,train_set,0.1);    
    predict_matrix_LP = 0.25* RPC + 0.25 * RPX +  0.23 * RRC + 0.1*RRP + 0.17*RRX;
    result(1,:) = result(1,:)+model_evaluate(interaction_matrix,predict_matrix_LP,train_set);   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LPLNP
    %proteins’ interaction profile  
    similairty_matrix=LNS_opt(train_set',0,6,'regulation2'); 
    predict_p_interaction=calculate_labels(similairty_matrix,train_set',0.5)'; 
    %proteins’ CTD
    similairty_matrix=LNS_opt(protein_ctd,0,23,'regulation2'); 
    predict_p_ctd=calculate_labels(similairty_matrix,train_set',0.3)';    
    %lncRNAs’ interaction profile 
    similairty_matrix=LNS_opt(train_set,0,100,'regulation2'); 
    predict_l_interaction=calculate_labels(similairty_matrix,train_set,0.7);    
    %lncRNA's expression profile 
    similairty_matrix=LNS_opt(extracted_lncRNA_expression,0,100,'regulation2');   
    predict_l_exp=calculate_labels(similairty_matrix,train_set,0.9);    
    %lncRNAs’ Sequence composition
    similairty_matrix=LNS_opt(extracted_lncRNA_sequence_CT,0,800,'regulation2'); 
    predict_l_seq=calculate_labels(similairty_matrix,train_set,0.1);
    predict_matrix_LP = 0.4 * predict_p_interaction + 0.1 * predict_p_ctd + 0.3 * predict_l_interaction + 0.19 * predict_l_seq + 0.01 * predict_l_exp;
    result(2,:) = result(2,:)+model_evaluate(interaction_matrix,predict_matrix_LP,train_set);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LPBNI
    LPBNI_score_matrix = LPBNI_opt(train_set);
    result(3,:) = result(3,:)+model_evaluate(interaction_matrix,LPBNI_score_matrix,train_set);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LPIHN
    gamma = 0.9;  beta = 0.9;  delta = 0.9;
    LPIHN_score_matrix = LPIHN_opt(extracted_lncRNA_expression, train_set', gamma, beta, delta);
    result(4,:) = result(4,:)+model_evaluate(interaction_matrix,LPIHN_score_matrix',train_set);  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIRWH
    alpha=0.4; l=2;  r=2;    
    BIRWH_score_matrix = BiRWH_opt(train_set,alpha,l,r);
    result(5,:) = result(5,:)+model_evaluate(interaction_matrix,BIRWH_score_matrix,train_set);  
end
result = result/CV;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  LKSNK
function W = KSNS_opt(feature_matrix,canshu)
    neighbor_num =canshu. neighbor_num;
    nearst_neighbor_matrix = calculate_neighbors(feature_matrix,neighbor_num);
    W = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu);
end
function nearst_neighbor_matrix=calculate_neighbors(feature_matrix,neighbor_num)
%%%%nearst_neighbor_matrix 
  X = feature_matrix;
  N = size(X,1);
  D = pdist2(X,X,'euclidean');  
  D = D+diag(inf*ones(1,N));
  [~, si]=sort(D,2,'ascend');
  nearst_neighbor_matrix=zeros(N,N);
  index=si(:,1:neighbor_num);
  for i=1:N
      nearst_neighbor_matrix(i,index(i,:))=1;    
  end
end
%%%%The process of KSNS
function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu)
X=feature_matrix';  
[~,N] = size(X); 
C = nearst_neighbor_matrix';
rand('state',1);
W = rand(N,N);
W = W./repmat(sum(W),N,1);
G  =jisuan_Kel(X,canshu);
G(isnan(G))=0;
G = G/max(G(:));
WC1 = W'*G*W-2*W*G+G;
WC = sum(diag(WC1))/2;
lata1=4;
lata2=1;
wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
objv = wucha;
error = canshu.jingdu*(1+lata1+lata2);     
we = 1;   
gen=1;
while  gen<1000 && we>error
    FZ = G+lata1*C.*W;
    FM = G*W+lata1*W+lata2*W;    
    W = FZ./(FM+eps).*W;  
    WC1 = W'*G*W-2*W*G+G;
    WC = sum(diag(WC1))/2;
    wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;    
    we = abs(wucha-objv(end));
    objv = [objv,wucha];
    gen = gen+1;
end
W = W./repmat(sum(W),N,1);
W=W';    
end

function G  =jisuan_Kel(X,canshu)
flag = canshu.flag;
KN = length(flag);
for i=1:KN
    if  flag(i)==1  
        G = Gaussian_Kel(X,canshu);
    elseif  flag(i) == 2
        G = GPolyn_Kel(X,canshu); 
    end
end
end

function W = Gaussian_Kel(X,canshu)
%%This program is referenced from Wang B
X = X'; 
sigma = canshu.sigma;
K=canshu.knum;
D = dist2(X,X);
D = ( D + D' )/2;
D= D- diag(diag(D));    
[T,~]=sort(D,2);
[~,n]=size(D);    
TT=mean(T(:,1:K),2)+eps;
Sig=(repmat(TT,1,n)+repmat(TT',n,1) + 1*D)/3;
Sig=Sig.*(Sig>eps)+eps;
W=normpdf(D,0,sigma*Sig);    
W = (W + W')/2;
end
function D = dist2(x, c)
[nx, dimx] = size(x);
[nc, dimc] = size(c);
if dimx ~= dimc
	error('Data dimension does not match dimension of centres')
end
tempx = full(sum(x.^2, 2));
tempc = full(sum(c.^2, 2)');
D = tempx(:, ones(1,nc)) + tempc(ones(1,nx), :) - 2.*(x*(c'));
end


function W  = GPolyn_Kel(X,canshu)
XX = X'*X;
DX = sqrt(diag(XX));
DXX = DX*DX';
X1 = XX./DXX; 
a = canshu.a;
b = canshu.b;
W = (X1+a).^b;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LPLNP(This program is referenced from Zhang Wen)
function W=LNS_opt(feature_matrix,tag,neighbor_num,regulation)
    distance_matrix=calculate_instances(feature_matrix);
    nearst_neighbor_matrix= LPLNP_neighbors(distance_matrix,neighbor_num);
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end
function W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation)
   row_num=size(feature_matrix,1);
   W=zeros(1,row_num);
   if tag==1
       row_num=1;
   end
   for i=1:row_num
       nearst_neighbors=feature_matrix(logical(nearst_neighbor_matrix(i,:)'),:);   
       neighbors_num=size(nearst_neighbors,1);
       G1=repmat(feature_matrix(i,:),neighbors_num,1)-nearst_neighbors;
       G2=repmat(feature_matrix(i,:),neighbors_num,1)'-nearst_neighbors';
       if regulation=='regulation2'
         G_i=G1*G2+eye(neighbors_num);
       end
       if regulation=='regulation1'
         G_i=G1*G2;
       end
       H=2*G_i;
       f=[];
       A=[];
       if isempty(H)
           A;
       end
       
       b=[];
       Aeq=ones(neighbors_num,1)';
       beq=1;
       lb=zeros(neighbors_num,1);
       ub=[];
       options=optimset('Display','off');
       [w,fval]= quadprog(H,f,A,b,Aeq,beq,lb,ub,[],options);
       w=w';
       W(i,logical(nearst_neighbor_matrix(i,:)))=w;     
   end
end

function distance_matrix=calculate_instances(feature_matrix)
    [row_num,col_num]=size(feature_matrix);
    distance_matrix=zeros(row_num,row_num);
    for i=1:row_num
        for j=i+1:row_num
            distance_matrix(i,j)=sqrt(sum((feature_matrix(i,:)-feature_matrix(j,:)).^2));
            distance_matrix(j,i)=distance_matrix(i,j);
        end
        distance_matrix(i,i)=col_num;
    end
end

function nearst_neighbor_matrix = LPLNP_neighbors(distance_matrix,neighbor_num)
  [sv si]=sort(distance_matrix,2,'ascend');
  [row_num,col_num]=size(distance_matrix);
  nearst_neighbor_matrix=zeros(row_num,col_num);
  index=si(:,1:neighbor_num);
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;
  end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LPBNI
function score_matrix = LPBNI_opt(interaction_matrix)
A = interaction_matrix';
[M,N] = size(A);
A1 = A./repmat(sum(A),M,1);
A1(isnan(A1))=0;
B1 = A'./repmat(sum(A'),N,1);
B1(isnan(B1))=0;
C = B1*A1;
score_matrix = C*interaction_matrix;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LPIHN  (This program is referenced from Zhang Wen)
function score_matrix = LPIHN_opt(feature_matrix, interaction_matrix, gamma, beta, delta)
    protein_similarity_matrix = get_protein_similarity_matrix(interaction_matrix);  
    rna_similarity_matrix = get_ran_similarity_matrix(feature_matrix); 
    row_sum_i = sum(interaction_matrix, 2);
    col_sum_i = sum(interaction_matrix);
    %get protein transformation matrix
    protein_num = size(protein_similarity_matrix);
    for i = 1 : protein_num
        protein_similarity_matrix(i, i) = 0;
    end
    row_sum_matrix_p = sum(protein_similarity_matrix, 2);
    sum_diagonal_matrix_p = pinv(diag(row_sum_matrix_p));
    transformation_matrix_p =  sum_diagonal_matrix_p * protein_similarity_matrix .* (1 - gamma);
    transformation_matrix_p(row_sum_i == 0, :) = transformation_matrix_p(row_sum_i == 0, :) ./ (1 - gamma);
    %get lncRNA transformation matrix
    rna_num = size(rna_similarity_matrix);
    for i = 1 : rna_num
        rna_similarity_matrix(i, i) = 0;
    end
    row_sum_matrix_l = sum(rna_similarity_matrix, 2);
    sum_diagonal_matrix_l = pinv(diag(row_sum_matrix_l));
    transformation_matrix_l = sum_diagonal_matrix_l * rna_similarity_matrix .* (1 - gamma);
    transformation_matrix_l(:, col_sum_i == 0) = transformation_matrix_l(:, col_sum_i == 0) ./ (1 - gamma);
    %get protein-lncRNA transformation matrix and lncRNA-protein transformation matrix
    row_sum_matrix_pl = row_sum_i;
    sum_diagonal_matrix_pl = pinv(diag(row_sum_matrix_pl));
    transformation_matrix_pl = sum_diagonal_matrix_pl * interaction_matrix .* gamma;
    %%%
    row_sum_matrix_lp = col_sum_i;
    sum_diagonal_matrix_lp = pinv(diag(row_sum_matrix_lp));
    transformation_matrix_lp = sum_diagonal_matrix_lp * interaction_matrix' .* gamma;
    transformation_matrix_lp(isnan(transformation_matrix_lp)) = 0;
    transformation_matrix_pl(isnan(transformation_matrix_pl)) = 0;
    %get transformation matrix
    transformation_matrix = [transformation_matrix_p transformation_matrix_pl; transformation_matrix_lp transformation_matrix_l];
    %get initial state
    rna_initial_state = eye(rna_num) * (1 - beta);
    protein_initial_state = interaction_matrix *  sum_diagonal_matrix_lp * beta;
    initial_state_matrix = [protein_initial_state; rna_initial_state];
    %get score matrix
    num = size(transformation_matrix, 1);
    score_matrix = (eye(num) / ((eye(num) - (1 - delta) * transformation_matrix'))) * delta * initial_state_matrix;
    score_matrix = score_matrix([1 : protein_num], :);
end

function rna_similarity_matrix = get_ran_similarity_matrix(feature_matrix)
%     rna_similarity_matrix = corrcoef(feature_matrix);  %%PCC相关性（就是相关系数）
%     rna_similarity_matrix(isnan(rna_similarity_matrix)) = 0;
%     rna_similarity_matrix = abs(rna_similarity_matrix);
    covfeature_matrix = feature_matrix*feature_matrix';
    rna_similarity_matrix = covfeature_matrix./(sqrt(diag(covfeature_matrix)*diag(covfeature_matrix)'));  %%PCC相关性（就是相关系数）
    rna_similarity_matrix(isnan(rna_similarity_matrix)) = 0;
    rna_similarity_matrix = abs(rna_similarity_matrix);
end

function protein_similarity_matrix = get_protein_similarity_matrix(interaction_matrix)
    %get intersection matrix
    intersection_matrix = interaction_matrix * interaction_matrix';
    %get denominator matrix
    protein_degree_matrix = sum(interaction_matrix, 2);
    denominator_matrix = sqrt(protein_degree_matrix * protein_degree_matrix');
    %calculate similarity_matrix of protein
    protein_similarity_matrix = intersection_matrix ./ denominator_matrix;
    protein_similarity_matrix(isnan(protein_similarity_matrix)) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BIRWH((This program is referenced from Zou S, Zhang J, Zhang Z)  )
function Rt = BiRWH_opt(A,alpha,l,r)
%bi-random walk on the heterogeneous network to calculate association scores for each lncRNA-protein pair. 
%A: adjacency matrix of the lncRNA-protein association network
%normFun: Laplacian normalization for lncRNA similarity and protein similarity  
[Wrr,Wpp] = Similarity_opt(A);
%Wrr: adjacency matrix of the lncRNA similarity network
%Wpp: adjacency matrix of the protein similarity network
normWrr = normFun(Wrr);
normWpp = normFun(Wpp);
%R0: initial probability
R0 = A/sum(A(:));
Rt = R0;
%bi-random walk on the heterogeneous network
for t=1:max(l,r)
    ftl = 0;
    ftp = 0;
    %random walk on the lncRNA similarity network
    if(t<=l)
        nRtleft = alpha *normWrr * Rt   + (1-alpha)*R0;
        ftl = 1;
    end
    %random walk on the protein similarity network
    if(t<=r)
        nRtright = alpha *   Rt * normWpp  + (1-alpha)*R0;
        ftp = 1;
    end
    %Rt: predictive association scores between each lncRNA-protein pair
    Rt =  (ftl*nRtleft + ftp*nRtright)/(ftl + ftp); 
end  
end

function [ result ] = normFun( M )
%normFun: Laplacian normalization
num = size(M,1);
nM = zeros(num,num);
result = zeros(num,num);
for i = 1:num
    nM(i,i) = sum(M(i,:));
end
for i = 1:num
    rsum = nM(i,i);
    for j = 1:num
        csum = nM(j,j);
        if((rsum==0)||(csum==0))
            result(i,j) = 0;
        else
            result(i,j) = M(i,j)/sqrt(rsum*csum);
        end
    end
end
    
end


function [Wll,Wpp] = Similarity_opt(train_set)
[nl,np] = size(train_set);  
gamall=1;  gamapp=1;
gamal = gamall/mean(diag(train_set*train_set'));
%calculate gamad for Gaussian kernel calculation
gamap = gamapp/mean(diag(train_set'*train_set));    
%calculate Gaussian train_set profile kernel similarity for lncRNA: wll
sA = (sum(train_set.^2, 2));
sB = sA;
Wll = exp(bsxfun(@minus,bsxfun(@minus,2*train_set*train_set', sA), sB')*gamal);  
%calculate Gaussian train_set profile kernel similarity for protein: pkd
train_set1 = train_set';
sA = (sum(train_set1.^2, 2));
sB = sA;
pkp = exp(bsxfun(@minus,bsxfun(@minus,2*train_set1*train_set1', sA), sB')*gamap);
%logistic function transformation for protein similarity
for i=1:np
    for j=1:np
        Wpp(i,j)=1/(1+exp(-15*pkp(i,j)+log(9999)));
    end
end

end


function F=calculate_labels(W,Y,alpha)
    F=(1-alpha)*pinv(eye(size(W,1))-alpha*W)*Y;
end

function result=model_evaluate(interaction_matrix,predict_matrix,train_ddi_matrix)
%%%%  (This program is referenced from Zhang Wen)
    real_score=interaction_matrix(:);
    predict_score=predict_matrix(:);
    index=train_ddi_matrix(:);
    test_index=find(index==0);
    real_score=real_score(test_index);
    predict_score=predict_score(test_index);
    aupr=AUPR(real_score,predict_score);
    auc=AUC(real_score,predict_score);
    [sen,spec,precision,accuracy,f1]=evaluation_metric(real_score,predict_score);
    result=[aupr,auc,sen,spec,precision,accuracy,f1];
end

function area=AUC(real,predict)
    max_value=max(predict);
    min_value=min(predict);
    threshold=min_value+(max_value-min_value)*(1:999)/1000;
    threshold=threshold';
    threshold_num=length(threshold);
    tn=zeros(threshold_num,1);
    tp=zeros(threshold_num,1);
    fn=zeros(threshold_num,1);
    fp=zeros(threshold_num,1);
    for i=1:threshold_num
        tp_index=logical(predict>=threshold(i) & real==1);
        tp(i,1)=sum(tp_index);

        tn_index=logical(predict<threshold(i) & real==0);
        tn(i,1)=sum(tn_index);

        fp_index=logical(predict>=threshold(i) & real==0);
        fp(i,1)=sum(fp_index);

        fn_index=logical(predict<threshold(i) & real==1);
        fn(i,1)=sum(fn_index);
    end

    sen=tp./(tp+fn);
    spe=tn./(tn+fp);
    y=sen;
    x=1-spe;
    [x,index]=sort(x);
    y=y(index,:);
    [y,index]=sort(y);
    x=x(index,:);

    area=0;
    x(threshold_num+1,1)=1;
    y(threshold_num+1,1)=1;
    area=0.5*x(1)*y(1);
    for i=1:threshold_num
        area=area+(y(i)+y(i+1))*(x(i+1)-x(i))/2;
    end
end
function area=AUPR(real,predict)
    max_value=max(predict);
    min_value=min(predict);

    threshold=min_value+(max_value-min_value)*(1:999)/1000;

    threshold=threshold';
    threshold_num=length(threshold);
    tn=zeros(threshold_num,1);
    tp=zeros(threshold_num,1);
    fn=zeros(threshold_num,1);
    fp=zeros(threshold_num,1);

    for i=1:threshold_num
        tp_index=logical(predict>=threshold(i) & real==1);
        tp(i,1)=sum(tp_index);
        tn_index=logical(predict<threshold(i) & real==0);
        tn(i,1)=sum(tn_index);
        fp_index=logical(predict>=threshold(i) & real==0);
        fp(i,1)=sum(fp_index);
        fn_index=logical(predict<threshold(i) & real==1);
        fn(i,1)=sum(fn_index);
    end
    sen=tp./(tp+fn);
    precision=tp./(tp+fp);
    recall=sen;
    x=recall;
    y=precision;
    [x,index]=sort(x);
    y=y(index,:);
    area=0;
    x(1,1)=0;
    y(1,1)=1;
    x(threshold_num+1,1)=1;
    y(threshold_num+1,1)=0;
    area=0.5*x(1)*(1+y(1));
    for i=1:threshold_num
        area=area+(y(i)+y(i+1))*(x(i+1)-x(i))/2;
    end
end

function [sen,spec,precision,accuracy,f1]=evaluation_metric(interaction_score,predict_score)
    max_value=max(predict_score);
    min_value=min(predict_score);
    threshold=min_value+(max_value-min_value)*(1:999)/1000;
    for i=1:999
       predict_label=(predict_score>threshold(i));
       [temp_sen(i),temp_spec(i),temp_precision(i),temp_accuracy(i),temp_f1(i)]=classification_metric(interaction_score,predict_label);
    end
    [max_score,index]=max(temp_f1);
    sen=temp_sen(index);
    spec=temp_spec(index);
    precision=temp_precision(index);
    accuracy=temp_accuracy(index);
    f1=temp_f1(index);
end

function [sen,spec,precision,accuracy,f1]=classification_metric(real_label,predict_label)
    tp_index=find(real_label==1 & predict_label==1);
    tp=size(tp_index,1);
    tn_index=find(real_label==0 & predict_label==0);
    tn=size(tn_index,1);
    fp_index=find(real_label==0 & predict_label==1);
    fp=size(fp_index,1);
    fn_index=find(real_label==1 & predict_label==0);
    fn=size(fn_index,1);
    accuracy=(tn+tp)/(tn+tp+fn+fp);
    sen=tp/(tp+fn);
    recall=sen;
    spec=tn/(tn+fp);
    precision=tp/(tp+fp);
    f1=2*recall*precision/(recall+precision);
end
