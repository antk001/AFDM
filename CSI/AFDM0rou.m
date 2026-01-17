clear; clc; close all;
N=16;
N_DFnT=16;
L=2; %Channel Length
Block_Num=2; %Block Number
%M=4; %Modulation QAM
C=2; %Len Cyclic Prefix 
P=N+C;
loop_Num=10;%1000000
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
%%
total=zeros(1,21,6);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,21,6);  %
for rou=logspace(-2, 0, 11); %不完美CSI
% for w= logspace(-2, 0, 11)
myfilename = strcat('OFDM-CSI=', num2str(rou),'-MMSE-BPSK-', num2str(N),'DFT(0.1).txt');
for dB=6
    disp(dB);
    SNR=10^(dB/10);%EB/N0
    count=1;
   % for Equal=1:2
     for Equal=2
      %  for U=1:3
         for U=1
            M=2^U;
            for loop=1:loop_Num
                [Bits,Symbols0]=Transmitter(N_DFnT,M,Block_Num,N,C);
                [h,Symbols1]=Channel(Symbols0,L,N,Block_Num,SNR);
                Bitsre=Receiver(rou,N_DFnT,M,Block_Num,N,C,Equal,Symbols1,h,SNR);             
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
plot(rou,total(:,:,1),'bo-');
set(gca,'Yscale','log');
ylim([1e-6 1]);
xlabel('rou (dB)');
ylabel('BER');
legend(' CSI c12=0 OFDM/MMSE BPSK')
grid on;


%save BER
%myfilename = strcat('AFDM-CSI=', num2str(rou),'-MMSE-BPSK-', num2str(N),'DFT(0.1).txt');
fid1=fopen(myfilename,'a');
  formatSpec = '∆=%4.4f BER =';
  fprintf(fid1,formatSpec,dB);
 fprintf(fid1,' %f',total(:,:,1));
 fprintf(fid1,'\r\n');
 fclose(fid1);
end