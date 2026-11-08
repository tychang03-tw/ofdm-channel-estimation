% OFDM 通道估測（早期版本）：線性內插、DFT、傳統 LS 比較
% 最大通道延遲 100 點、衰減係數 α = 0.05，傳統 LS 的修正參數 λ = 10^-1
% 輸出：BER vs. SNR
clear
clc
num = 20 ;       %幾組 OFDM symbol 
length1 = 2048;  %長度
                 %20*2048矩陣
cp = 128;      
signal=zeros(num, length1);
c = zeros(1, num * (length1 + cp));
fft_signal=zeros(num, length1);
ypilot=zeros(num,141);
pilot=zeros(num,141);
H=zeros(num,length1);
H2=zeros(num,length1);
H3=zeros(num,141);
H1_r=zeros(num,length1);
X=zeros(num,length1);
num1=zeros(num,841);
num2=zeros(num,840);
num4=zeros(num,841);
num5=zeros(num,840);
db=1:20;

for i = 1:num
    x = sign(randn(1,841)) + 1i * sign(randn(1, 841));
    x(1209:2048)=sign(randn(1,840)) + 1i * sign(randn(1, 840));
    signal(i,:)=x;
     
    %pilot
    for a=1:71
        pilot(i,a)=signal(i,1+12*(a-1));
    end
    for b=1:70
        pilot(i,1+70+b)=signal(i,1209+12*(b-1));
    end

    %ifft
    y=signal(i,1:length1);
    ifft_signal = ifft(y);

    %cp
    guard=ifft_signal((length1 -cp+1) :length1);
    ifft_cpsignal = [guard,ifft_signal];%後16位元加上x64=80位元
    c(1, 1 + 2176*(i-1) : 2176*i) = ifft_cpsignal;
end

%20個通道
for j=1:20
    n=randn(1,100) + 1i * randn(1, 100);
    g=exp(-0.05*[0:99]).*n;
    C=conv(c,g);
    E=(C*C')/43619;

    %20個錯誤率
    for snr=1:length(db)
        noise=randn(1, 43619)+1i.*randn(1, 43619);
        N =(0.5 .* 10^((-db(snr)) / 10))*E;
        Y = C + noise .*sqrt(N);

        %移除保護區間
        for i = 1:num
            resignal=Y(1,129+2176*(i-1):2176*i);
            %fft
            fft_signal(i,:)=fft(resignal);
            %pilot
            for a=1:71
                ypilot(i,a)=fft_signal(i,1+12*(a-1));
            end
            for b=1:70
                ypilot(i,1+70+b)=fft_signal(i,1209+12*(b-1));
            end
        end
        %case1 線性內插 求完整H
        %生成Hpilot
        for i=1:num
            for a=1:71
                H(i,1+12*(a-1))=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H(i,1209+12*(b-1))=ypilot(i,71+b)/pilot(i,71+b);
            end
        end

        for i=1:num
            for q=1:70
                for p=1:11
                    H(i,1+12*(q-1)+p)=((12-p)/12).*H(i,1+12*(q-1))+(p/12).*H(i,1+12*q);
                end
            end
            for r=1:69
                for p=1:11
                    H(i,1209+12*(r-1)+p)=((12-p)/12).*H(i,1209+12*(r-1))+(p/12).*H(i,1209+12*r);
                end
            end
            for p=1:11
                H(i,2037+p)=((12-p)/12).*H(i,2037)+(p/12).*H(i,1);
            end
        end
        X([1:num],[1:length1])=fft_signal([1:num],[1:length1])./H([1:num],[1:length1]);         
      
        %case2 DFT   把每一橫行放進DFT
        %生成Hpilot
        for i=1:num
            for a=1:71
                H2(i,1+12*(a-1))=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H2(i,1209+12*(b-1))=ypilot(i,71+b)/pilot(i,71+b);
            end
        end

        for i = 1:num
            D(i,[1:length1]) = ifft(H2(i,[1:length1]));
        end

        D([1:num],[142:length1])=0;
        T([1:num],[1:length1])=12.*D;

        for i = 1:num
            Do(i,[1:length1]) = fft(T(i,[1:length1]));
        end
        DFT([1:num],[1:length1])=fft_signal([1:num],[1:length1])./Do([1:num],[1:length1]);

        % %case3 LS法
        Le=100;
        f1 = exp(-1i * 2 * pi * (12.*((0:70).')) * (0:Le-1) / 2048);% k =12*( 0:70&71:140 ) ;N = 2048;
        f2 = exp(-1i * 2 * pi * (1209+12.*((0:69).')) * (0:Le-1) / 2048);
        f = [f1;f2];

        %生成Hpilot
        for i=1:num
            for a=1:71
                H3(i,a)=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H3(i,71+b)=ypilot(i,71+b)/pilot(i,71+b);
            end
        end

        h_ls=zeros(Le,num);

        for i=1:num
            hn([1:141],1)=H3(i,[1:141]);
            h_ls([1:Le],i)=inv(f'*f+10^-1.*eye(Le))*f'*hn;
            h_LS=h_ls.';
        end

        h_LS([1:num],[142:length1])=0;

        for i = 1:num
            h_LS1(i,[1:length1]) = fft(h_LS(i,[1:length1]));
        end
        h_LS2([1:num],[1:length1])=fft_signal([1:num],[1:length1])./h_LS1([1:num],[1:length1]);

        %case1 錯誤率  每一dB的平均錯誤率

        X([1:20],1+12*([1:71]-1))=0;
        X([1:20],1209+12*([1:70]-1))=0;
        signal([1:20],1+12*([1:71]-1))=0;
        signal([1:20],1209+12*([1:70]-1))=0;

        X_real = sign(real(X));
        X_imag = sign(imag(X));

        for i = 1:num
            num1(i, :) = (X_real(i, 1:841) ~= real(signal(i, 1:841))) + (X_imag(i, 1:841) ~= imag(signal(i, 1:841)));
            num2(i, :) = (X_real(i, 1209:2048) ~= real(signal(i, 1209:2048))) + (X_imag(i, 1209:2048) ~= imag(signal(i, 1209:2048)));
        end
        num3=sum(sum(num1))+sum(sum(num2));
        error_rate1 = num3 /(1540*2*20) ;
        error1(:,snr)=error_rate1;

        %case2 錯誤率

        DFT([1:20],1+12*([1:71]-1))=0;
        DFT([1:20],1209+12*([1:70]-1))=0;
        signal([1:20],1+12*([1:71]-1))=0;
        signal([1:20],1209+12*([1:70]-1))=0;

        DFT_real = sign(real(DFT));
        DFT_imag = sign(imag(DFT));
        for i = 1:num
            num4(i, :) = (DFT_real(i, 1:841) ~= real(signal(i, 1:841))) + (DFT_imag(i, 1:841) ~= imag(signal(i, 1:841)));
            num5(i, :) = (DFT_real(i, 1209:2048) ~= real(signal(i, 1209:2048))) + (DFT_imag(i, 1209:2048) ~= imag(signal(i, 1209:2048)));
        end
        num6=sum(sum(num4))+sum(sum(num5));
        error_rate2 = num6 /(1540*2*20) ;
        error2(:,snr)=error_rate2;

        %case3 錯誤率

        h_LS2([1:20],1+12*([1:71]-1))=0;
        h_LS2([1:20],1209+12*([1:70]-1))=0;
        signal([1:20],1+12*([1:71]-1))=0;
        signal([1:20],1209+12*([1:70]-1))=0;

        h_LS2_real = sign(real(h_LS2));
        h_LS2_imag = sign(imag(h_LS2));
        for i = 1:num
            num7(i, :) = (h_LS2_real(i, 1:841) ~= real(signal(i, 1:841))) + (h_LS2_imag(i, 1:841) ~= imag(signal(i, 1:841)));
            num8(i, :) = (h_LS2_real(i, 1209:2048) ~= real(signal(i, 1209:2048))) + (h_LS2_imag(i, 1209:2048) ~= imag(signal(i, 1209:2048)));
        end
        num9=sum(sum(num7))+sum(sum(num8));
        error_rate8 = num9 /(1540*2*20) ;
        error8(:,snr)=error_rate8;
    end

    %case1
    error3(j,:)=error1;
    %case2
    error4(j,:)=error2;
    %case3
    error9(j,:)=error8;
end

%把每個通道的1~20dB的錯誤率取平均
%內插法
Linear_e=sum(error3)./20;
%DFT法
DFT_e=sum(error4)./20;
%LS法
LS_e=sum(error9)./20;
semilogy(db,Linear_e, 'b-',db,DFT_e, 'r-',db,LS_e,'g-');
legend('Linear_e','DFT_e','LS_e');
xlabel('(dB)');
ylabel('BER');
grid;

% a1=g;
% k=1:2048;
% a1(2048)=0;
% a1f=fft(a1);
% plot(k,abs(a1f),k,abs(h_LS1(1,:)));





