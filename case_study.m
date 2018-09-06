function case_study
clc
clear 
load YCdata.mat;
TR_LJM =  YCdata.TR_LJM;
NP = size(TR_LJM,2);
M=20;
%%%%%%%%%%LKSNS
[PLP,RLP] = LKSNS_LP(TR_LJM); 
W = 0.75;
P = W*PLP+(1-W)*RLP;
[RANK,real_LJ]=  jisuan_opt(P,M,YCdata);
LKSNS_result = cell(20,2);
for i=1:20
    LKSNS_result{i,1} = real_LJ{i}(8:end);
    index = find(RANK==i);    
    if length(index)>0
        LKSNS_result{i,2}='T';
    else
        LKSNS_result{i,2}= '';
    end
end
% xlswrite('LKSNS_result.xls',LKSNS_result)

%%%%%%%%%%LPLNP
[PLP,RLP] = LPLNP_LP(TR_LJM); 
W = 0.65;
P = W*PLP+(1-W)*RLP;
[RANK,real_LJ]=  jisuan_opt(P,M,YCdata);
LPLNP_result = cell(20,2);
for i=1:20
    LPLNP_result{i,1} = real_LJ{i}(8:end);
    index = find(RANK==i);    
    if length(index)>0
        LPLNP_result{i,2}='T';
    else
        LPLNP_result{i,2}='';
    end
end
% xlswrite('LPLNP_result.xls',LPLNP_result)




%%%%%%%LPIHN
LP = LPIHN_LP(TR_LJM,[0.7000,0.9000,0.9000]);
[RANK,real_LJ] =  jisuan_opt(LP,M,YCdata);
LPIHN_result = cell(20,2);
for i=1:20
    LPIHN_result{i,1} = real_LJ{i}(8:end);
    index = find(RANK==i);    
    if length(index)>0
        LPIHN_result{i,2}='Y';
    else
        LPIHN_result{i,2}='';
    end
end
% xlswrite('LPIHN_result.xls',LPIHN_result)


%%%%%%%LPBNI
LP = LPBNI_LP(TR_LJM);
[RANK,real_LJ] =  jisuan_opt(LP,M,YCdata);
LPBNI_result = cell(20,2);
for i=1:20
    LPBNI_result{i,1} = real_LJ{i}(8:end);
    index = find(RANK==i);    
    if length(index)>0
        LPBNI_result{i,2}='Y';
    else
        LPBNI_result{i,2}='';
    end
end
% xlswrite('LPBNI_result.xls',LPBNI_result)

%%%%%%BIRWH
LP = BIRWH_LP(TR_LJM,[0.5,2,2]);
[RANK,real_LJ] =  jisuan_opt(LP,M,YCdata);
BIRWH_result = cell(20,2);
for i=1:20
    BIRWH_result{i,1} = real_LJ{i}(8:end);
    index = find(RANK==i);    
    if length(index)>0
        BIRWH_result{i,2}='Y';
    else
        BIRWH_result{i,2}='';
    end
end
% xlswrite('BIRWH_result.xls',BIRWH_result)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LKSNS model
function [PLP,RLP] = LKSNS_LP(TR_LJM)
interaction_matrix = TR_LJM;
[NR,NP] = size(interaction_matrix);
%%% proteins°Ø interaction profile
canshu.neighbor_num = 2;canshu.knum  = 4;  canshu.sigma  = 0.5;
canshu.flag = 1;   canshu.jingdu = 0.00001;    
SP = Label_Propagation(interaction_matrix',canshu);      
PLP = calculate_labels(SP,interaction_matrix',0.6);
PLP = PLP';
%%% lncRNAs°Ø interaction profile
canshu.lata1 = 4 ;  canshu.lata2 = 1 ; canshu.neighbor_num = floor(0.9*NR);    
canshu.knum  = floor(0.8*NR);  canshu.sigma  = 0.8;
canshu.jingdu = 0.0001;  canshu.flag= 1;
SR = Label_Propagation(interaction_matrix,canshu);      
RLP = calculate_labels(SR,interaction_matrix,0.1);  
end  
%%%Compute the similarity matrix
function W=Label_Propagation(feature_matrix,canshu)
    neighbor_num =canshu. neighbor_num;
    nearst_neighbor_matrix = calculate_neighbors(feature_matrix,neighbor_num);
    W = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu);
end

%%%%Calculate the neighbors
function nearst_neighbor_matrix=calculate_neighbors(feature_matrix,neighbor_num)
  X = feature_matrix;
  N = size(X,1);
  D = pdist2(X,X,'euclidean');  
  D = D+diag(inf*ones(1,N));
  [~, si]=sort(D,2,'ascend');
  nearst_neighbor_matrix=zeros(N,N);
  index=si(:,1:neighbor_num);
  for i=1:N
      nearst_neighbor_matrix(i,index(i,:))=1;     %%%º∆À„¡⁄æ”æÿ’Û
  end
end
%%%%%The process of KSNS
function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu)
lata1 = 4;
lata2 = 1;
X=feature_matrix';  
[M,N] = size(X);
C = nearst_neighbor_matrix';
rand('state',1);
W = rand(N,N);
W = W./repmat(sum(W),N,1);
G = gaussian_Kel(X,canshu);
G = G/max(G(:));
WC1 = W'*G*W-2*W*G+G;
WC = sum(diag(WC1))/2;
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

function W = gaussian_Kel(X,canshu)
%%This program is referenced from Wang B
X = X'; 
sigma = canshu.sigma;
K = canshu.knum;
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
function F=calculate_labels(W,Y,alpha)
    F=(1-alpha)*pinv(eye(size(W,1))-alpha*W)*Y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LPBNI
function LP = LPBNI_LP(interaction_matrix)
A = interaction_matrix';
[M,N] = size(A);
A1 = A./repmat(sum(A),M,1);
A1(isnan(A1))=0;
B1 = A'./repmat(sum(A'),N,1);
B1(isnan(B1))=0;
C = B1*A1;
LP = C*interaction_matrix;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LPIHN
%This program is referenced from Zhang wen
function LP = LPIHN_LP(TR_LJM,abc)
a = abc(1);
b = abc(2);
c = abc(3);
predict_matrix_hrwr = heterogeneous_matrix(TR_LJM', a, b, c);
LP = predict_matrix_hrwr';
end  
 
function score_matrix = heterogeneous_matrix(interaction_matrix, gamma, beta, delta)
    protein_similarity_matrix = get_protein_similarity_matrix(interaction_matrix);
    rna_similarity_matrix = get_ran_similarity_matrix(interaction_matrix);
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
function rna_similarity_matrix = get_ran_similarity_matrix(interaction_matrix)
    rna_similarity_matrix = corrcoef(interaction_matrix);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%LPLNP
%This program is referenced from Zhang wen
function [PLP,RLP] = LPLNP_LP(TR_LJM)
interaction_matrix = TR_LJM;
[NR,NP] = size(interaction_matrix);
%%%%proteins°Ø interaction profile
similairty_matrix=LPLNP_Propagation(interaction_matrix',0,11,'regulation2');    
predict_l_interaction=calculate_labels(similairty_matrix,interaction_matrix',0.1);
PLP =predict_l_interaction';   
%%%lncRNAs°Ø interaction profile
similairty_matrix=LPLNP_Propagation(interaction_matrix,0,floor(NR*0.7),'regulation2');    
RLP =calculate_labels(similairty_matrix,interaction_matrix,0.1);
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
function nearst_neighbor_matrix=LPLNP_neighbors(distance_matrix,neighbor_num)
  [sv si]=sort(distance_matrix,2,'ascend');
  [row_num,col_num]=size(distance_matrix);
  nearst_neighbor_matrix=zeros(row_num,col_num);
  index=si(:,1:neighbor_num);
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;
  end
end
function W=LPLNP_Propagation(feature_matrix,tag,neighbor_num,regulation)
    distance_matrix=calculate_instances(feature_matrix);
    nearst_neighbor_matrix=LPLNP_neighbors(distance_matrix,neighbor_num);
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BIRWH
function LP = BIRWH_LP(TR_LJM,abc)
alpha = abc(1);
l = abc(2);
r = abc(3);
[Wrr,Wpp] = Similarity_opt(TR_LJM);
LP = BiRWHMDA_opt(Wrr,Wpp,TR_LJM,alpha,l,r);
end

function Rt = BiRWHMDA_opt(Wrr,Wpp,A,alpha,l,r)
%BiRWHMDA(0.4,2,2) (This program is referenced from Zou S, Zhang J, Zhang Z)  
%bi-random walk on the heterogeneous network to calculate association scores for each microbe-disease pair. 
%Wrr: adjacency matrix of the lncRNA similarity network
%Wpp: adjacency matrix of the protein similarity network
%A: adjacency matrix of the lncRNA-protein association network
%normFun: Laplacian normalization for lncRNA similarity andprotein similarity
normWrr = normFun(Wrr);
normWpp = normFun(Wpp);
%R0: initial probability
R0 = A/sum(A(:));
Rt = R0;
%bi-random walk on the heterogeneous network
for t=1:max(l,r)    
    ftl = 0;
    ftr = 0;
    %random walk on the lncRNA similarity network
    if(t<=l)
        nRtleft = alpha * Rt * normWrr + (1-alpha)*R0;
        ftl = 1;
    end
    %random walk on theprotein similarity network
    if(t<=r)
        nRtright = alpha *  normWpp * Rt + (1-alpha)*R0;
        ftr = 1;
    end
    Rt =  (ftl*nRtleft + ftr*nRtright)/(ftl + ftr);    
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
function [Wrr,Wpp] = Similarity_opt(train_set)
[nd,nl] = size(train_set); 
gamall=1;
gamadd=1;
%calculate gamal for Gaussian kernel calculation
for i=1:nl
    sl(i)=norm(train_set(:,i))^2;
end
gamal=nl/sum(sl')*gamall;
%calculate gamad for Gaussian kernel calculation
for i=1:nd
    sd(i)=norm(train_set(i,:))^2;
end
gamad=nd/sum(sd')*gamadd;
%calculate Gaussian train_set profile kernel similarity for lncRNA: Wrr
Wrr = zeros(nl,nl);
for i=1:nl
    for j=1:nl
        Wrr(i,j)=exp(-gamal*(norm(train_set(:,i)-train_set(:,j)))^2);
    end
end
%calculate Gaussian train_set profile kernel similarity for disease: pkd
for i=1:nd
    for j=1:nd
        pkd(i,j)=exp(-gamad*(norm(train_set(i,:)-train_set(j,:)))^2);
    end
end 
%   Wpp = pkd;
%logistic function transformation for disease similarity
for i=1:nd
    for j=1:nd
        Wpp(i,j)=1/(1+exp(-15*pkd(i,j)+log(9999)));
    end
end
end


function [RANK,real_LJ]=  jisuan_opt(P,M,YCdata)
%%%%input
%P   The interaction score matrix
%M   The number of interaction with the highest score
%%%%output
%RANK     The confirmed indexs in the first M interactions
%%LJ      The name of the first M interactions
TQ_LJ_stasebate = YCdata.TQ_LJ_stasebate; 
TQ_LJ_inter3 = YCdata.TQ_LJ_inter3;
ZL_GY_inter2 = YCdata.ZL_GY_inter2;
name = YCdata.name;
TR_LJM =  YCdata.TR_LJM;
UR = name.UR;
UP = name.UP;
%%%%%Calculate the first M interactions
PLJ = jisuan_LJ(P,M,TR_LJM);
M = size(PLJ,1);
YC = cell(M,2);
for i=1:M
    YC{i,1} = UR{PLJ(i,1)};
    YC{i,2} = UP{PLJ(i,2)};
end
%%%%%Judge each connection using NPInter v3.0 and stasebate
T=0;
RANK = [];
real_LJ=[];
for i=1:M    
    index = find(strcmp(ZL_GY_inter2(:,1), YC{i,1}));
    real_LJ = [real_LJ;{[YC{i,1},'--',YC{i,2}]}];
    flag=0;
    for j=1:length(index)
        a = [ZL_GY_inter2{index(j),2},'--',YC{i,2}];
        index1 = find(strcmp(TQ_LJ_inter3, a));
        if length(index1)>0
            T=T+1;
            RANK = [RANK;i];
            flag=1;            
            break;
        end
    end
    if  flag == 0 
        a = [YC{i,1},'--',YC{i,2}];
        index1 = find(strcmp(TQ_LJ_stasebate, a));
        if length(index1)>0
            T=T+1;
            RANK = [RANK;i];
        end
    end
end


end

    
function PLJ = jisuan_LJ(P,M,TR_LJM)
%%%%output
%PLJ The first M interactions
[NR,~] = size(P);
index=TR_LJM(:);   
test_index = find(index==0);    
predict_score=P(:);
predict_score = predict_score(test_index);
[PXP,index] = sort(predict_score,'descend');
xy = test_index(index(1:M)); 
rx = zeros(M,1);  
py = zeros(M,1);  
for i=1:M
    rx(i) = mod(xy(i),NR);
    if rx(i)==0
        rx(i)=NR;
    end
    py(i) = (xy(i)-rx(i))/NR+1;
end
PLJ = [rx,py];
end   
