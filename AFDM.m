clear; clc; close all;
N=4; %子载波数
N_DFnT=4;
L=2; %信道长度
Block_Num=2; %块数
%M=4; %QAM 调制
C=2; %循坏前缀长度
P=N+C;
loop_Num=10000;%1000000
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
%%
total=zeros(1,21,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,21,6);  
for dB=0:2:40
    disp(dB);
    SNR=10^(dB/10);%EB/N0
    count=1;
   % for Equal=1:2
     for Equal=1
      %  for U=1:3
         for U=1
            M=2^U;
            for loop=1:loop_Num
                [Bits,Symbols0]=Transmitter(N_DFnT,M,Block_Num,N,C);
                [H0,Symbols1]=Channel(Symbols0,L,N,Block_Num,SNR);
                Bitsre=Receiver(N_DFnT,M,Block_Num,N,C,Equal,Symbols1,H0,SNR);             
                ratio(1,dB/2+1,count)=sum(Bits~=Bitsre)/(Block_Num*N*log2(M));
                total(1,dB/2+1,count)=total(1,dB/2+1,count)+ratio(1,dB/2+1,count);
            end
            count=count+1;
        end
    end
end
total=total/loop_Num;
figure()
box on; hold on;
plot(0:2:40,total(:,:,1),'bo-');
set(gca,'Yscale','log');
ylim([1e-6 1]);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('AFDM/MMSE BPSK')
grid on;