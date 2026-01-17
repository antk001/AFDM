function Bitsre=Receiver(N_DFnT,M,Block_Num,N,C,Equal,Symbols1,H0,SNR)
%1=ZF, 2=MMSE
% %%%%%%%%%%%%%%%%%%%%%%%%%%5
% M=16;%QAM 
% Block_Num=2;
% N=128;
% C=4;

%%%%%%%%%%%%%%%%%%%%%%%%%%5
P=N+C;
R=[zeros(N,P-N),eye(N)];
Symbols2=zeros(N,1,Block_Num);
%% 去除循环前缀
for a=1:Block_Num
    Symbols2(:,:,a)=R*Symbols1(:,:,a);
end
%% 构造和应用DFT矩阵 
%% Discrete DA Transform at Transmitter
% c1=1/(2*N_DFnT);
% c2=1/(2*N_DFnT);
c1=5/(2*N_DFnT);
c2=1.414;
W_IDFT=zeros(N_DFnT);
for a=1:N_DFnT
   for b=1:N_DFnT
       W_IDFT(a,b)=exp(1i*2*pi*(a-1)*(b-1)/N_DFnT);
   end 
end
W_DFT=W_IDFT*1/sqrt(N_DFnT);
WH_DFT=conj(W_IDFT);

T_1=zeros(N_DFnT);
for a=1:N_DFnT
   for b=1:N_DFnT
       if a==b
       T_1(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c1);
       end
   end 
end
T_1=T_1*exp(-1i*pi/4);
T_1H=conj(T_1);

T_2=zeros(N_DFnT);
for a=1:N_DFnT
    for b=1:N_DFnT
        if a==b
           T_2(a,b)=exp(1i*2*pi*(a-1)*(b-1)*c2);
        end
    end 
end
T_2 =T_2*1/sqrt(N_DFnT);
T_2H=conj(T_2);

Q_IDFnT0=T_2H*W_IDFT*T_1H;
P_DFnT0=conj(Q_IDFnT0);
  


Symbols3=zeros(size(Symbols2));
% 经1/√N归一化的 IDFT 矩阵
for count=1:Block_Num
    sample=P_DFnT0*Symbols2(:,:,count); 
    Symbols3(:,:,count)=sample;
end

%% 构造均衡矩阵
S=eye(N);
T=[S(2*N-P+1:N,:);S];
H=R*H0*T;
D=P_DFnT0*H*Q_IDFnT0; %%对角矩阵 

%%版本9
D1=(D').*D;
if Equal==1
    G=inv(D);
end
if Equal==2
   G=(D').*((D1+(1/SNR).*eye(size(D)))^(-1));
   %G=(D').*(pinv(D1+(1/SNR).*eye(size(D))));
  % G=pinv((D)'*D + eye(size(D)).*(1/SNR))*(D)';
end
G=pinv((D)'*D + eye(size(D)).*(1/SNR))*(D)';
%% 均衡
Symbols4=zeros(size(Symbols3));
for count=1:Block_Num
    Symbols4(:,:,count)=G*Symbols3(:,:,count);
end

 
%% 解调
%%BPSK
for count=1:Block_Num    
    Symbols6(:,:,count)=pskdemod(Symbols4(:,:,count),M,pi/2);  % Qpsk
% Symbols6(:,:,count)=qamdemod(Symbols4(:,:,count)/sqrt(1/10),M);
end
% if M==4
%     Symbols6=qamdemod(Symbols5/sqrt(1/2),M);
% end
% if M==16
%     Symbols6=qamdemod(Symbols5/sqrt(1/10),M);
% end
% if M==64
%     Symbols6=qamdemod(Symbols5/sqrt(1/42),M);
% end
Bitsre=zeros(1,N*Block_Num*log2(M));
start=1;
for count=1:Block_Num
    for k=1:N
        dec=dec2bin(Symbols6(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end
%end








 






