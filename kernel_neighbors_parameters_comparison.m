%The prediction results of LKSNS model under different kernel functions,
% number of neighbors and propagation parameters
function kernel_neighbors_parameters_comparison
clc
clear 
num = 0.1:0.2:0.9;    %Neighbors proportion
ar =  0.1:0.1:0.9;    %Transmission parameters
Nn = length(num);
Na = length(ar);
numar = [];
for i=1:Nn
    for j=1:Na
        numar = [numar;[num(i),ar(j)]];     
    end
end
N = size(numar,1);

%%%%%%proteins°Ø interaction profile
Gaussian =[];  sca_exp =[]; 
poly =[];  cos_poly =[]; 
FSWeight =[]; 
for i=1:N
    jieguo=computer_PC(numar(i,1),numar(i,2));
    Gaussian = [Gaussian;jieguo(1,:)];
    sca_exp = [sca_exp;jieguo(2,:)];
    poly =[poly;jieguo(3,:)];
    cos_poly =[cos_poly;jieguo(4,:)];
    FSWeight =[FSWeight;jieguo(5,:)];        
end
resultPC.Gaussian = Gaussian;
resultPC.sca_exp = sca_exp;
resultPC.poly = poly;
resultPC.cos_poly = cos_poly;
resultPC.FSWeight = FSWeight;
save resultPC resultPC;

%%%%%%proteins°Ø CTD
Gaussian =[];  sca_exp =[]; 
poly =[];  cos_poly =[]; 
cosine =[]; 
parfor i=1:N
    jieguo=computer_PX(numar(i,1),numar(i,2));
    Gaussian = [Gaussian;jieguo(1,:)];
    sca_exp = [sca_exp;jieguo(2,:)];
    poly =[poly;jieguo(3,:)];
    cos_poly =[cos_poly;jieguo(4,:)];
    cosine =[cosine;jieguo(5,:)];        
end
resultPX.Gaussian = Gaussian;
resultPX.sca_exp = sca_exp;
resultPX.poly = poly;
resultPX.cos_poly = cos_poly;
resultPX.cosine = cosine;
save resultPX resultPX;

%%%%%%lncRNAs°Ø interaction profile
Gaussian =[];  sca_exp =[]; 
poly =[];  cos_poly =[]; 
FSWeight =[]; 
parfor i=1:N
    result = computer_RC(numar(i,1),numar(i,2));
    Gaussian = [Gaussian;result.Gaussian];
    sca_exp = [sca_exp;result.sca_exp];
    poly =[poly;result.poly];
    cos_poly =[cos_poly;result.cos_poly];
    FSWeight =[FSWeight;result.FSWeight];        
end
resultRC.Gaussian = Gaussian;
resultRC.sca_exp = sca_exp;
resultRC.poly = poly;
resultRC.cos_poly = cos_poly;
resultRC.FSWeight = FSWeight;
save resultRC resultRC;

%%%%%%lncRNA's expression profile
Gaussian =[];  sca_exp =[]; 
poly =[];  cos_poly =[]; 
cosine =[]; 
parfor i=1:N
    result = computer_RP(numar(i,1),numar(i,2));
    Gaussian = [Gaussian;result.Gaussian];
    sca_exp = [sca_exp;result.sca_exp];
    poly =[poly;result.poly];
    cos_poly =[cos_poly;result.cos_poly];
    cosine =[cosine;result.cosine];        
end
resultRP.Gaussian = Gaussian;
resultRP.sca_exp = sca_exp;
resultRP.poly = poly;
resultRP.cos_poly = cos_poly;
resultRP.cosine = cosine;
save resultRP resultRP;

%%%%%%lncRNAs°Ø Sequence composition
Gaussian =[];  sca_exp =[]; 
poly =[];  cos_poly =[]; 
cosine =[]; 
parfor i=1:N
    result = computer_RX(numar(i,1),numar(i,2));
    Gaussian = [Gaussian;result.Gaussian];
    sca_exp = [sca_exp;result.sca_exp];
    poly =[poly;result.poly];
    cos_poly =[cos_poly;result.cos_poly];
    cosine =[cosine;result.cosine];        
end
resultRX.Gaussian = Gaussian;
resultRX.sca_exp = sca_exp;
resultRX.poly = poly;
resultRX.cos_poly = cos_poly;
resultRX.cosine = cosine;
save resultRX resultRX;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%proteins°Ø interaction profile
function result = computer_PC(num,ar)
seed = 1;
load extracted_interaction.txt;
interaction_matrix = extracted_interaction;
CV=5;
rand('state',seed);
[row,col]=size(interaction_matrix);
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);                                                   
[~,NP] = size(extracted_interaction);
result=zeros(5,7);
for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    if (k~=CV)
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    train_set=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_set(test_row_index(i),test_col_index(i))=0;                 
    end
    %%%%Gaussian
    canshu.flag = 1; 
    canshu.gama = 1.2; canshu.jingdu = 0.00001;  
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(train_set',canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(1,:) = result(1,:) +model_evaluate(interaction_matrix,Rt,train_set);     
    %sca_exp
    canshu.flag = 2; 
    canshu.sigma = 0.5;  canshu.knum = 9;   canshu.jingdu = 0.00001;  
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(train_set',canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(2,:) = result(2,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %poly
    canshu.flag = 3; 
    canshu.a =0.1;  canshu.b= 2;   canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(train_set',canshu,lata1,lata2);     
    Rt = calculate_labels(SP,train_set',ar)';
    result(3,:) = result(3,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %%cos-poly
    canshu.flag = 4; 
    canshu.a =0.2;  canshu.b= 4;   canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(train_set',canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(4,:) = result(4,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %%%%FSweight
    canshu.flag = 5; 
    canshu.jingdu = 0.00001;  lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(train_set',canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(5,:) = result(5,:)+model_evaluate(interaction_matrix,Rt,train_set);    
end
result=result/CV;
disp('result')
result
result = [num*ones(5,1),ar*ones(5,1),result];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%proteins°Ø CTD
function result = computer_PX(num,ar)
seed = 1;
load extracted_interaction.txt;
load protein_ctd.mat;
interaction_matrix = extracted_interaction;
CV=5;
rand('state',seed);
[row,col]=size(interaction_matrix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);                                                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,NP] = size(extracted_interaction);
result=zeros(5,7);
for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
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
    %%%%gaussian
    canshu.flag = 1; 
    canshu.gama = 24; canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(protein_ctd,canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(1,:) = result(1,:) +model_evaluate(interaction_matrix,Rt,train_set);
    %sca_exp
    canshu.flag = 2; 
    canshu.sigma = 0.3 ;  canshu.knum = 17;   canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(protein_ctd,canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(2,:) = result(2,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %poly
    canshu.flag = 3; 
    canshu.a = 0.1;  canshu.b= 4;   canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(protein_ctd,canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(3,:) = result(3,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %%cos_poly
    canshu.flag = 4; 
    canshu.a =0.1;  canshu.b= 4;   canshu.jingdu = 0.00001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NP);  
    SP = Label_Propagation(protein_ctd,canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(4,:) = result(4,:)+model_evaluate(interaction_matrix,Rt,train_set);
    %%%%cosine
    canshu.flag = 6; 
    lata1 = 4;  lata2= 1;  canshu.jingdu = 0.00001;
    SP = Label_Propagation(protein_ctd,canshu,lata1,lata2);      
    Rt = calculate_labels(SP,train_set',ar)';
    result(5,:) = result(5,:)+model_evaluate(interaction_matrix,Rt,train_set);    
end
result=result/CV;
disp('result')
result
result = [num*ones(5,1),ar*ones(5,1),result];    
end  
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNAs°Ø Sequence composition
function result = computer_RX(num,ar)
seed = 1;
load extracted_interaction.txt;
load extracted_lncRNA_sequence_CT.txt;
interaction_matrix = extracted_interaction;
CV=5;
rand('state',seed);
[row,col]=size(interaction_matrix);
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);                                                  
[NR,~] = size(extracted_interaction);
Na = length(ar);
Gaussian=zeros(Na,7);  sca_exp=zeros(Na,7); 
poly=zeros(Na,7);  cos_poly=zeros(Na,7); 
cosine=zeros(Na,7);

for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    if (k~=CV)
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    train_set=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_set(test_row_index(i),test_col_index(i))=0;                 
    end
    %%%%Gaussian
    canshu.flag=1;
    canshu.gama = 1;    canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRX = Label_Propagation(extracted_lncRNA_sequence_CT,canshu,lata1,lata2);          
    Gaussian = Gaussian + biaoqian_opt(SRX,train_set,ar,interaction_matrix);
    %sca_exp
    canshu.flag = 2; 
    canshu.sigma = 0.9 ;  canshu.knum = 900;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRX = Label_Propagation(extracted_lncRNA_sequence_CT,canshu,lata1,lata2);        
    sca_exp = sca_exp+ biaoqian_opt(SRX,train_set,ar,interaction_matrix);
    %poly
    canshu.flag = 3; 
    canshu.a = 0.1  ;  canshu.b=4 ;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRX = Label_Propagation(extracted_lncRNA_sequence_CT,canshu,lata1,lata2);     
    poly = poly + biaoqian_opt(SRX,train_set,ar,interaction_matrix);    
    %%cos_poly
    canshu.flag = 4; 
    canshu.a = 0.1;  canshu.b= 4;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRX = Label_Propagation(extracted_lncRNA_sequence_CT,canshu,lata1,lata2);      
    cos_poly = cos_poly + biaoqian_opt(SRX,train_set,ar,interaction_matrix);
    %%%%cosine
    canshu.flag = 6; 
    lata1 = 4;  lata2= 1  ;  canshu.jingdu = 0.0001;
    SRX = Label_Propagation(extracted_lncRNA_sequence_CT,canshu,lata1,lata2);     
    cosine = cosine+ biaoqian_opt(SRX,train_set,ar,interaction_matrix);   
end
disp('Gaussian')
Gaussian=[num*ones(Na,1),ar(:),Gaussian/5];
Gaussian
disp('sca_exp')
sca_exp=[num*ones(Na,1),ar(:),sca_exp/5];
sca_exp
disp('poly')
poly=[num*ones(Na,1),ar(:),poly/5];
poly
disp('cos_poly')
cos_poly=[num*ones(Na,1),ar(:),cos_poly/5];
cos_poly
disp('cosine')
cosine = [num*ones(Na,1),ar(:),cosine/5];
cosine

result.Gaussian = Gaussian;
result.sca_exp = sca_exp;
result.poly = poly;
result.cos_poly = cos_poly;
result.cosine = cosine;
end

function jieguo = biaoqian_opt(SR,train_set,ar,interaction_matrix)
Na = length(ar);
jieguo = zeros(Na,7);
for i=1:Na
    Rt = calculate_labels(SR,train_set,ar(i));
    jieguo(i,:) = model_evaluate(interaction_matrix,Rt,train_set);
end
end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNAs°Ø interaction profile    
function result = computer_RC(num,ar)
seed = 1;
load extracted_interaction.txt;
interaction_matrix = extracted_interaction;
CV=5;
rand('state',seed);
[row,col]=size(interaction_matrix);
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);                                                    
[NR,~] = size(extracted_interaction);
Na = length(ar);
Gaussian=zeros(Na,7);  sca_exp=zeros(Na,7); 
poly=zeros(Na,7);  cos_poly=zeros(Na,7); 
FSWeight=zeros(Na,7);

for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    if (k~=CV)
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    train_set=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_set(test_row_index(i),test_col_index(i))=0;                 
    end
    %%%%Gaussian
    canshu.flag=1;
    canshu.gama = 0.3;    canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SR = Label_Propagation(train_set,canshu,lata1,lata2);               
    Gaussian = Gaussian + biaoqian_opt(SR,train_set,ar,interaction_matrix);     
    %sca_exp
    canshu.flag = 2; 
    canshu.sigma = 0.5;  canshu.knum= 900;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SR = Label_Propagation(train_set,canshu,lata1,lata2);          
    sca_exp = sca_exp+ biaoqian_opt(SR,train_set,ar,interaction_matrix);    
    %poly
    canshu.flag = 3; 
    canshu.a = 0.9;  canshu.b= 2;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SR = Label_Propagation(train_set,canshu,lata1,lata2);      
    poly = poly + biaoqian_opt(SR,train_set,ar,interaction_matrix);    
    %%poly
    canshu.flag = 4; 
    canshu.a = 0.9;  canshu.b= 2;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SR = Label_Propagation(train_set,canshu,lata1,lata2);     
    cos_poly = cos_poly + biaoqian_opt(SR,train_set,ar,interaction_matrix);        
    %%%%FSweight
    canshu.flag = 5; 
    canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SR = Label_Propagation(train_set,canshu,lata1,lata2);      
    FSWeight = FSWeight+ biaoqian_opt(SR,train_set,ar,interaction_matrix);   
end
disp('Gaussian')
Gaussian=[num*ones(Na,1),ar(:),Gaussian/5];
Gaussian
disp('sca_exp')
sca_exp=[num*ones(Na,1),ar(:),sca_exp/5];
sca_exp
disp('poly')
poly=[num*ones(Na,1),ar(:),poly/5];
poly
disp('cos_poly')
cos_poly=[num*ones(Na,1),ar(:),cos_poly/5];
cos_poly
disp('FSWeight')
FSWeight = [num*ones(Na,1),ar(:),FSWeight/5];
FSWeight

result.Gaussian = Gaussian;
result.sca_exp = sca_exp;
result.poly = poly;
result.cos_poly = cos_poly;
result.FSWeight = FSWeight;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lncRNA's expression profile
function result = computer_RP(num,ar)
seed = 1;
load extracted_interaction.txt;
load extracted_lncRNA_expression.txt;
interaction_matrix = extracted_interaction;
CV=5;
rand('state',seed);
[row,col]=size(interaction_matrix);
[row_index,col_index]=find(interaction_matrix==1);
link_num=sum(sum(interaction_matrix)); 
rand('state',seed);
random_index=randperm(link_num);
size_of_CV=round(link_num/CV);                                                 
[NR,~] = size(extracted_interaction);
Na = length(ar);
Gaussian=zeros(Na,7);  sca_exp=zeros(Na,7); 
poly=zeros(Na,7);  cos_poly=zeros(Na,7); 
cosine=zeros(Na,7);
for k=1:CV
    fprintf('begin to implement the cross validation:round =%d/%d\n', k, CV);
    if (k~=CV)
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):(size_of_CV*k)));
    else
        test_row_index=row_index(random_index((size_of_CV*(k-1)+1):end));
        test_col_index=col_index(random_index((size_of_CV*(k-1)+1):end));
    end
    train_set=interaction_matrix;
    test_link_num=size(test_row_index,1);
    for i=1:test_link_num
        train_set(test_row_index(i),test_col_index(i))=0;                 
    end
    %%%%Gaussian
    canshu.flag=1;
    canshu.gama = 7;    canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRP = Label_Propagation(extracted_lncRNA_expression,canshu,lata1,lata2);    
    Gaussian = Gaussian + biaoqian_opt(SRP,train_set,ar,interaction_matrix);     
    %sca_exp
    canshu.flag = 2; 
    canshu.sigma = 0.5 ;  canshu.knum = 900;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRP = Label_Propagation(extracted_lncRNA_expression,canshu,lata1,lata2);     
    sca_exp = sca_exp+ biaoqian_opt(SRP,train_set,ar,interaction_matrix);    
    %poly
    canshu.flag = 3; 
    canshu.a = 0.1;  canshu.b= 2;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRP = Label_Propagation(extracted_lncRNA_expression,canshu,lata1,lata2); 
    poly = poly + biaoqian_opt(SRP,train_set,ar,interaction_matrix);    
    %%cos_poly
    canshu.flag = 4; 
    canshu.a = 0.9;  canshu.b= 2.0;   canshu.jingdu = 0.0001; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  
    SRP = Label_Propagation(extracted_lncRNA_expression,canshu,lata1,lata2);   
    cos_poly = cos_poly + biaoqian_opt(SRP,train_set,ar,interaction_matrix);        
    %cosine
    canshu.flag = 6; 
    lata1 = 4;  lata2= 1  ; canshu.neighbor_num = floor(num*NR);  canshu.jingdu = 0.0001; 
    SRP = Label_Propagation(extracted_lncRNA_expression,canshu,lata1,lata2);  
    cosine = cosine+ biaoqian_opt(SRP,train_set,ar,interaction_matrix);   
end
disp('Gaussian')
Gaussian=[num*ones(Na,1),ar(:),Gaussian/5];
Gaussian
disp('sca_exp')
sca_exp=[num*ones(Na,1),ar(:),sca_exp/5];
sca_exp
disp('poly')
poly=[num*ones(Na,1),ar(:),poly/5];
poly
disp('cos_poly')
cos_poly=[num*ones(Na,1),ar(:),cos_poly/5];
cos_poly
disp('cosine')
cosine = [num*ones(Na,1),ar(:),cosine/5];
cosine

result.Gaussian = Gaussian;
result.sca_exp = sca_exp;
result.poly = poly;
result.cos_poly = cos_poly;
result.cosine = cosine;

end


%%%Calculate KSNS similarity
function W=Label_Propagation(feature_matrix,canshu,lata1,lata2)
    neighbor_num =canshu. neighbor_num;
    nearst_neighbor_matrix = calculate_neighbors(feature_matrix,neighbor_num);
    W = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu,lata1,lata2);
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
      nearst_neighbor_matrix(i,index(i,:))=1;     %%%º∆À„¡⁄æ”æÿ’Û
  end
end

%%The KSNS iteration process
function [W,objv] = jisuanW_opt(feature_matrix,nearst_neighbor_matrix,canshu,lata1,lata2)
X=feature_matrix';  
[M,N] = size(X);
C = nearst_neighbor_matrix';
rand('state',1);
 W = rand(N,N);
 W = W./repmat(sum(W),N,1);
G = jisuan_Kel(X,canshu);
G = G/max(G(:));
G(isnan(G))=0;
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
W(isnan(W))=0;
end    
%%%%%Kernel Function
function K = jisuan_Kel(X,canshu)
flag = canshu.flag;
X = X';
if flag==1  
    K = gaosi_kel(X,canshu);
elseif flag==2
    K = sca_exp_kel(X,canshu); 
elseif flag==3 
    K = poly_kel(X,canshu); 
elseif flag==4
    K = cos_poly(X,canshu); 
elseif flag==5 
    K = FSweight_kel(X); 
elseif flag==6
    K = cos_kel(X); 
end
end
function K = gaosi_kel(X,canshu)
gama = canshu.gama;
D = pdist2(X,X);
L = pdist2(X,zeros(size(X)),'euclidean');
L = L(:,1);
sig = mean(L);
K = exp(-D/sig*gama);
end

function K = sca_exp_kel(X,canshu)
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
K=normpdf(D,0,sigma*Sig);    
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


function K = poly_kel(X,canshu)
a = canshu.a;
b = canshu.b;
K = (X*X'+a).^b;
end

function  K = cos_poly(X,canshu)
X = X';
XX = X'*X;
DX = sqrt(diag(XX));
DXX = DX*DX';
X1 = XX./DXX; 
a = canshu.a;
b = canshu.b;
K = (X1+a).^b;

K(isnan(K)) = 0; 
end


function K = FSweight_kel(X)
N = size(X,1);
K = zeros(N);
for i=1:N
    for j=1:N
        JI= sum(X(i,:).*X(j,:)); 
        P = X(i,:)-X(j,:);
        N1 = length(find(P==1));
        N2 = length(find(P==-1));
        W1 = 2*JI/(N1+2*JI+1);
        W2 = 2*JI/(N2+2*JI+1);
        K(i,j)  = W1*W2;
    end
end
K(isnan(K))=0;
end

function K  = cos_kel(X)
XX = X*X';
DX = sqrt(diag(XX));
DXX = DX*DX';
K = XX./DXX; 
K(isnan(K))=0;
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
    
    
    
   
    
    

    
    


