function [Bits,Symbols0]=Transmitter(N_DFnT,M,Block_Num,N,C)
P=N+C;                              
%% Initalize Symbols, Unit Power    
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)        
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
%% BPSK
 Bits3=pskmod(Bits2,M,pi/2);%bpsk
% if M==4
%     Bits3=qammod(Bits2,M)*sqrt(0.5);%4QAM
% end
% if M==16
%     Bits3=qammod(Bits2,M)*sqrt(1/10);
% end
% if M==64
%     Bits3=qammod(Bits2,M)*sqrt(1/42);
% end
%% Discrete DA Transform at Transmitter
%case1 ：AFDM
c1=1/(2*N_DFnT);
c2=1/(2*N_DFnT);
c1=5/(2*N_DFnT);
c2=1.414;
c1=9/(2*N_DFnT);
c2=0.4;
%case2 ：OFDM
c1=0;
c2=0;
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
 
%% IDFnT 
Symbols=zeros(N,1,Block_Num);
Block=1;
if mod(N,2)==0
    for k=1:N:length(Bits3)
        Symbols(:,:,Block)=Q_IDFnT0*Bits3(k:k+N-1).';
        Block=Block+1;
    end
end
if mod(N,2)==1
    for k=1:N:length(Bits3)
        Symbols(:,:,Block)=Q_IDFnT0*Bits3(k:k+N-1).';
        Block=Block+1;
    end
end
%% Cyclic Prefix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
Symbols0=zeros(P,1,Block_Num);
for a=1:Block_Num
    Symbols0(:,:,a)=T*Symbols(:,:,a);
end
end









