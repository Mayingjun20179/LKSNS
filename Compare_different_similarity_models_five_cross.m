function Compare_different_similarity_models_five_cross
clc
clear 
ar = 0.1:0.1:0.9;
[cosine,FSWeight,Gaussian,Jaccard,LNS,KSNS]= five_cross(ar);
save cosine cosine;
save FSWeight FSWeight;
save Gaussian Gaussian;
save Jaccard Jaccard;
save LNS LNS;
save KSNS KSNS;
end


function [cosine,FSWeight,Gaussian,Jaccard,LNS,KSNS]= five_cross(ar)
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
L = length(ar);
cosine.PC = zeros(L,7); cosine.PX = zeros(L,7); cosine.RC = zeros(L,7); cosine.RP = zeros(L,7); cosine.RX = zeros(L,7);
FSWeight.PC = zeros(L,7); FSWeight.PX = zeros(L,7); FSWeight.RC = zeros(L,7); FSWeight.RP = zeros(L,7); FSWeight.RX = zeros(L,7);
Gaussian.PC = zeros(L,7); Gaussian.PX = zeros(L,7); Gaussian.RC = zeros(L,7); Gaussian.RP = zeros(L,7); Gaussian.RX = zeros(L,7);
Jaccard.PC = zeros(L,7); Jaccard.PX = zeros(L,7); Jaccard.RC = zeros(L,7); Jaccard.RP = zeros(L,7); Jaccard.RX = zeros(L,7);
LNS.PC = zeros(L,7); LNS.PX = zeros(L,7); LNS.RC = zeros(L,7); LNS.RP = zeros(L,7); LNS.RX = zeros(L,7);
KSNS.PC = zeros(L,7); KSNS.PX = zeros(L,7); KSNS.RC = zeros(L,7); KSNS.RP = zeros(L,7); KSNS.RX = zeros(L,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    %%extracts test datasets and training datasets
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
    %%%%%%%%%%%%%%%%%%%%%%%%    proteins¡¯ interaction profile
    %%cosine
    cosine_SPC = similarity_opt(train_set',1); 
    result = Label_result(cosine_SPC,train_set',ar,interaction_matrix');
    cosine.PC = cosine.PC + result;
    %%FSWeight
    FSWeight_SPC = similarity_opt(train_set',2);      
    result = Label_result(FSWeight_SPC,train_set',ar,interaction_matrix');
    FSWeight.PC = FSWeight.PC + result;
    %%Gaussian
    Gaussian_SPC = similarity_opt(train_set',3);   
    result = Label_result(Gaussian_SPC,train_set',ar,interaction_matrix');
    Gaussian.PC = Gaussian.PC + result;
    %%Jaccard
    Jaccard_SPC = similarity_opt(train_set',4);      
    result = Label_result(Jaccard_SPC,train_set',ar,interaction_matrix');
    Jaccard.PC = Jaccard.PC + result;
    %%LNS  
    LNS_SPC = similarity_opt(train_set',5,0,6,'regulation2');      
    result = Label_result(LNS_SPC,train_set',ar,interaction_matrix');
    LNS.PC = LNS.PC + result;      
    %%KSNS
    canshu.neighbor_num = floor(0.3*NP);  canshu.sigma = 0.5;  canshu.knum = 9;    
    canshu.jingdu = 0.00001;    canshu.flag = 1; 
    KSNS_SPC = similarity_opt(train_set',6,[],[],[],canshu);      
    result = Label_result(KSNS_SPC,train_set',ar,interaction_matrix');
    KSNS.PC = KSNS.PC + result;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%proteins¡¯ CTD
    %%cosineine
    cosine_SPX = similarity_opt(protein_ctd,1);      
    result = Label_result(cosine_SPX,train_set',ar,interaction_matrix');
    cosine.PX = cosine.PX + result;
    %%FSWeight  (NULL)

    %%Gaussian
    Gaussian_SPX = similarity_opt(protein_ctd,3);      
    result = Label_result(Gaussian_SPX,train_set',ar,interaction_matrix');
    Gaussian.PX = Gaussian.PX + result;
    %%Jaccard  (NULL)

    %%LNS    
    LNS_SPX = similarity_opt(protein_ctd,5,0,23,'regulation2');      
    result = Label_result(LNS_SPX,train_set',ar,interaction_matrix');
    LNS.PX = LNS.PX + result;  
    
    %%KSNS
    canshu.neighbor_num = floor(0.5*NP);    canshu.sigma = 0.3;  canshu.knum = 17;  
    canshu.flag = 1;    canshu.jingdu = 0.00001; 
    KSNS_SPX = similarity_opt(protein_ctd,6,[],[],[],canshu);   
    result = Label_result(KSNS_SPX,train_set',ar,interaction_matrix');
    KSNS.PX = KSNS.PX + result;
  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNAs¡¯ interaction profile
    %%cosine
    cosine_SRC = similarity_opt(train_set,1); 
    result = Label_result(cosine_SRC,train_set,ar,interaction_matrix);
    cosine.RC = cosine.RC + result;   
    %%FSWeight
    FSWeight_SRC = similarity_opt(train_set,2);      
    result = Label_result(FSWeight_SRC,train_set,ar,interaction_matrix);
    FSWeight.RC = FSWeight.RC + result; 
    %Gaussian
    Gaussian_SRC = similarity_opt(train_set,3);      
    result = Label_result(Gaussian_SRC,train_set,ar,interaction_matrix);
    Gaussian.RC = Gaussian.RC + result; 
    %%Jaccard
    Jaccard_SRC = similarity_opt(train_set,4);      
    result = Label_result(Jaccard_SRC,train_set,ar,interaction_matrix);
    Jaccard.RC = Jaccard.RC + result; 
    %%LNS_opt    
    LNS_SRC = similarity_opt(train_set,5,0,100,'regulation2');   
    result = Label_result(LNS_SRC,train_set,ar,interaction_matrix);
    LNS.RC = LNS.RC + result; 
    %%KSNS
    canshu.neighbor_num = floor(0.1*NR);    canshu.sigma = 0.5;  canshu.knum=900;         
    canshu.jingdu = 0.0001;  canshu.flag= 1;
    KSNS_SRC = similarity_opt(train_set,6,[],[],[],canshu);     
    result = Label_result(KSNS_SRC,train_set,ar,interaction_matrix);
    KSNS.RC = KSNS.RC + result;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNA's expression profile
    %%cosine
    cosine_SRP = similarity_opt(extracted_lncRNA_expression,1); 
    result = Label_result(cosine_SRP,train_set,ar,interaction_matrix);
    cosine.RP = cosine.RP + result;  
    %%FSWeight  £¨NULL£©

    %Gaussian
    Gaussian_SRP = similarity_opt(extracted_lncRNA_expression,3);      
    result = Label_result(Gaussian_SRP,train_set,ar,interaction_matrix);
    Gaussian.RP = Gaussian.RP + result;  
    %%Jaccard   £¨NULL£©
    
    
    %%LNS
    LNS_SRP = similarity_opt(extracted_lncRNA_expression,5,0,100,'regulation2');     
    result = Label_result(LNS_SRP,train_set,ar,interaction_matrix);
    LNS.RP = LNS.RP + result;   
    %%KSNS
    canshu.neighbor_num = floor(0.1*NR);    canshu.a = 0.9;    canshu.b = 2;
    canshu.flag= 2 ;  canshu.jingdu =  0.0001;
    KSNS_SRP = similarity_opt(extracted_lncRNA_expression,6,[],[],[],canshu);      
    result = Label_result(KSNS_SRP,train_set,ar,interaction_matrix);
    KSNS.RP = KSNS.RP + result;   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNAs¡¯ Sequence composition
    %%cosine
    cosine_SRX = similarity_opt(extracted_lncRNA_sequence_CT,1);      
    result = Label_result(cosine_SRX,train_set,ar,interaction_matrix);
    cosine.RX = cosine.RX + result; 
    %%FSWeight  £¨NULL£©
    
    
    %%Gaussian
    Gaussian_SRX = similarity_opt(extracted_lncRNA_sequence_CT,3);      
    result = Label_result(Gaussian_SRX,train_set,ar,interaction_matrix);
    Gaussian.RX = Gaussian.RX + result; 
    %%Jaccard   £¨NULL£©
    
    
    %%LNS_opt    
    LNS_SRX = similarity_opt(extracted_lncRNA_sequence_CT,5,0,800,'regulation2');     
    result = Label_result(LNS_SRX,train_set,ar,interaction_matrix);
    LNS.RX = LNS.RX + result; 
    %%KSNS
    canshu.neighbor_num = floor(0.7*NR);    canshu.a = 0.1; canshu.b = 4;
    canshu.jingdu = 0.0001;  canshu.flag=2;
    KSNS_SRX = similarity_opt(extracted_lncRNA_sequence_CT,6,[],[],[],canshu);      
    result = Label_result(KSNS_SRX,train_set,ar,interaction_matrix);
    KSNS.RX = KSNS.RX + result;     
    
end
cosine.PC = cosine.PC/CV; cosine.PX = cosine.PX/CV ; cosine.RC = cosine.RC/CV ; cosine.RP = cosine.RP/CV; cosine.RX = cosine.RX/CV;
FSWeight.PC = FSWeight.PC/CV; FSWeight.PX = FSWeight.PX/CV; FSWeight.RC = FSWeight.RC/CV; FSWeight.RP = FSWeight.RP/CV; FSWeight.RX = FSWeight.RX/CV;
Gaussian.PC = Gaussian.PC/CV; Gaussian.PX = Gaussian.PX/CV; Gaussian.RC = Gaussian.RC/CV; Gaussian.RP = Gaussian.RP/CV; Gaussian.RX = Gaussian.RX/CV;
Jaccard.PC = Jaccard.PC/CV; Jaccard.PX = Jaccard.PX/CV; Jaccard.RC = Jaccard.RC/CV; Jaccard.RP = Jaccard.RP/CV; Jaccard.RX = Jaccard.RX/CV;
LNS.PC = LNS.PC/CV; LNS.PX = LNS.PX/CV; LNS.RC = LNS.RC/CV; LNS.RP = LNS.RP/CV; LNS.RX = LNS.RX/CV;
KSNS.PC = KSNS.PC/CV; KSNS.PX = KSNS.PX/CV; KSNS.RC = KSNS.RC/CV; KSNS.RP = KSNS.RP/CV; KSNS.RX = KSNS.RX/CV;

end

function result = Label_result(S,train_set,ar,interaction_matrix)
%%%%Calculat the prediction scores of different propagation parameters 
L = length(ar);
result = zeros(L,7);
for i=1:L
    score = calculate_labels(S,train_set,ar(i));  
    result(i,:) = model_evaluate(interaction_matrix,score,train_set);
end
end

function F=calculate_labels(W,Y,alpha)
%Label propagation  (The program is referenced from the Zhang W)
    F=(1-alpha)*pinv(eye(size(W,1))-alpha*W)*Y;
end
    
function S = similarity_opt(X,flag,tag,neighbor_num,regulation,canshu)
if flag==1
    S = cos_opt(X);
elseif flag==2
    S = FSWeight_opt(X);
elseif flag==3
    S = Gaussian_opt(X);
elseif flag==4
    S = Jaccard_opt(X); 
elseif flag==5
    S = LNS_opt(X,tag,neighbor_num,regulation); 
elseif flag == 6
    S = KSNS_opt(X,canshu);  
end
end

%%%%%%Cosine similarity
function S = cos_opt(X)
%%X feature matrix £¨Refer to zhou zhi hua£©
W = X*X';
DX = sqrt(diag(W))*sqrt(diag(W))';
W = W./DX;


W(isnan(W))=0;
D = sum(W,2);
index = find(D==0);
D = ones(size(D))./sqrt(D);  
D(index)=0;
D(isnan(D))=0;
D = diag(D);
S = D*W*D;

end


%%%%%FSWeight similarity
function S = FSWeight_opt(X)
N = size(X,1);
W = zeros(N);
for i=1:N
    for j=1:N
        JI= sum(X(i,:).*X(j,:)); 
        P = X(i,:)-X(j,:);
        N1 = length(find(P==1));
        N2 = length(find(P==-1));
        W1 = 2*JI/(N1+2*JI+1);
        W2 = 2*JI/(N2+2*JI+1);
        W(i,j)  = W1*W2;
    end
end
W(isnan(W))=0;
D = sum(W,2);
index = find(D==0);
D = ones(size(D))./sqrt(D);   
D(index)=0;
D(isnan(D))=0;
D = diag(D);
S = D*W*D;
end


function S = Gaussian_opt(X)
beta = mean(diag(X*X'));
sA = (sum(X.^2, 2));
sB = (sum(X.^2, 2));
W = exp(bsxfun(@minus,bsxfun(@minus,2*X*X', sA), sB')/beta);

W(isnan(W))=0;
[row,col]=size(W);
for i=1:row
    W(i,i)=0;
end
for round=1:20
    D=diag(sum(W,2));
    D1=pinv(sqrt(D));
    W=D1*W*D1;
end
S=W;
end


function S = Jaccard_opt(X)
N = size(X,1);
W = zeros(N);
for i=1:N
    for j=1:N
        FZ = X(i,:).*X(j,:);  
        FZ = sum(FZ);
        FM = X(i,:)+X(j,:);    
        FM = length(find(FM>0)); 
        W(i,j)  = FZ/FM;
    end
end
W(isnan(W))=0;
D = sum(W,2);
index = find(D==0);
D = ones(size(D))./sqrt(D);   
D(index)=0;
D(isnan(D))=0;
D = diag(D);
S = D*W*D;
end

function W = LNS_opt(feature_matrix,tag,neighbor_num,regulation)
    distance_matrix=calculate_instances(feature_matrix);
    nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num);
    W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W=optimization_similairty_matrix(feature_matrix,nearst_neighbor_matrix,tag,regulation)
%%%%%%The program is referenced from the Zhang W
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

function nearst_neighbor_matrix=calculate_neighbors(distance_matrix,neighbor_num)
  [sv si]=sort(distance_matrix,2,'ascend');
  [row_num,col_num]=size(distance_matrix);
  nearst_neighbor_matrix=zeros(row_num,col_num);
  index=si(:,1:neighbor_num);
  for i=1:row_num
       nearst_neighbor_matrix(i,index(i,:))=1;
  end
end

%%%%%%%%%%% KSNS similarity
function W = KSNS_opt(feature_matrix,canshu)
    neighbor_num =canshu. neighbor_num;
    nearst_neighbor_matrix = KSNS_neighbors(feature_matrix,neighbor_num);
    W = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu);
end

function nearst_neighbor_matrix=KSNS_neighbors(feature_matrix,neighbor_num)
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


function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu)
lata1=4;
lata2=1;
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
wucha = WC + norm(W.*(1-C),'fro')^2*lata1/2 +  norm(W,'fro')^2*lata2/2;
objv = wucha;
error = canshu.jingdu*(1+lata1+lata2);    
we = 1;   
gen=1;
while  gen<1000 && we>error
    %gen
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
% toc
W=W';    
end


function W  =jisuan_Kel(X,canshu)
flag = canshu.flag;
KN = length(flag);
for i=1:KN
if  flag(i)==1       
    W = Gaussian_Kel(X,canshu);
elseif  flag(i) == 2 
    W = GPolyn_Kel(X,canshu); 
end
end


end

function W = Gaussian_Kel(X,canshu)
%%%%% The program is referenced from the  Wang B
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




function result=model_evaluate(interaction_matrix,predict_matrix,train_ddi_matrix)
%%%%%%%The program is referenced from the Zhang W
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % plot(x,y)
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
   
    
    
    
