% 主程式 OFDM 通道估測：四種方法 BER 比較（最大通道延遲 100 點）
% 20 個 QPSK OFDM 符元（N = 2048、CP = 128、每 12 個子載波一個領航），20 組隨機多路徑通道 x SNR 1～20 dB；比較線性內插、DFT、傳統 LS（λ = 10^-5）、時域 LS（λ = 10^-5 與無修正參數）。衰減係數 α 在 g = exp(-α*[0:99]).*n 修改
% 輸出：BER vs. SNR
clear
clc

OFDM_symbol_num = 20 ;
OFDM_symbol_length = 2048;
cp = 128;

xsignal = zeros(OFDM_symbol_num, OFDM_symbol_length);
c = zeros(1, OFDM_symbol_num * (OFDM_symbol_length + cp));

pilot = zeros(OFDM_symbol_num,141);
ypilot = zeros(OFDM_symbol_num,141);
ysignal = zeros(OFDM_symbol_num, OFDM_symbol_length);
resignal = zeros(OFDM_symbol_num, OFDM_symbol_length);

H = zeros(OFDM_symbol_num,OFDM_symbol_length);
X = zeros(OFDM_symbol_num,OFDM_symbol_length);
num1 = zeros(OFDM_symbol_num,841);
num2 = zeros(OFDM_symbol_num,840);

H2 = zeros(OFDM_symbol_num,OFDM_symbol_length);
D = zeros(OFDM_symbol_num,OFDM_symbol_length);
Do = zeros(OFDM_symbol_num,OFDM_symbol_length);
num4 = zeros(OFDM_symbol_num,841);
num5 = zeros(OFDM_symbol_num,840);

H3 = zeros(OFDM_symbol_num,141);
num7 = zeros(OFDM_symbol_num,841);
num8 = zeros(OFDM_symbol_num,840);

ifft_xpilot = zeros(OFDM_symbol_num,OFDM_symbol_length);
ifft_hpilot = zeros(OFDM_symbol_num,OFDM_symbol_length);
fft_HP_LS = zeros(OFDM_symbol_num,OFDM_symbol_length);
X_LS2 = zeros(OFDM_symbol_num,OFDM_symbol_length);
xpilot = zeros(OFDM_symbol_num,OFDM_symbol_length);
hpilot = zeros(OFDM_symbol_num,OFDM_symbol_length);

num10 = zeros(OFDM_symbol_num,841);
num11 = zeros(OFDM_symbol_num,840);

% 生成QPSK調製的OFDM符號和插入pilot序列
for i = 1:OFDM_symbol_num
    x = sign(randn(1,841)) + 1i * sign(randn(1, 841));
    x(1209:2048) = sign(randn(1,840)) + 1i * sign(randn(1, 840));
    xsignal(i,:) = x;
    %----------------------------------------------------------------
    %case4
    for a = 1:71
        xpilot(i,1+12*(a-1)) = xsignal(i,1+12*(a-1));
    end
    for b = 1:70
        xpilot(i,1209+12*(b-1)) = xsignal(i,1209+12*(b-1));
    end
    ifft_xpilot(i,1:OFDM_symbol_length) = ifft(xpilot(i,1:OFDM_symbol_length));
    for k = 0:99
        shift = [0 k];%shift 點數為通道點數
        x_bar([1:2048],k+1,i) = circshift(ifft_xpilot(i,[1:2048]),shift);
    end
    %----------------------------------------------------------------
    %pilot
    for a = 1:71
        pilot(i,a) = xsignal(i,1+12*(a-1));
    end
    for b = 1:70
        pilot(i,1+70+b) = xsignal(i,1209+12*(b-1));
    end
    %ifft
    ifft_xsignal = ifft(xsignal(i,1:OFDM_symbol_length));
    %cp
    guard = ifft_xsignal((OFDM_symbol_length -cp+1):OFDM_symbol_length);
    ifft_cpxsignal = [guard,ifft_xsignal];%後16位元加上x64=80位元
    c(1, 1 + 2176*(i-1) : 2176*i) = ifft_cpxsignal;
end
%20個通道
for j = 1:20
    n = randn(1,100) + 1i * randn(1, 100);
    g = exp(-0.001*[0:99]).*n;
    C = conv(c,g);
    E = (C*C')/43619;
    db = 1:20;
    %20個錯誤率
    for snr = 1:length(db)
        noise = randn(1, 43619)+1i.*randn(1, 43619);
        N = (0.5 .* 10^((-db(snr)) / 10))*E;
        Y = C + noise .*sqrt(N);
        %移除保護區間
        for i = 1:OFDM_symbol_num
            resignal=Y(1,129+2176*(i-1):2176*i);
            %fft
            ysignal(i,:)=fft(resignal);
            %pilot
            for a=1:71
                ypilot(i,a)=ysignal(i,1+12*(a-1));
            end
            for b=1:70
                ypilot(i,1+70+b)=ysignal(i,1209+12*(b-1));
            end
        end
        %----------------------------------------------------------------
        %case1 線性內插 求完整H
        %生成Hpilot
        for i=1:OFDM_symbol_num
            for a=1:71
                H(i,1+12*(a-1))=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H(i,1209+12*(b-1))=ypilot(i,71+b)/pilot(i,71+b);
            end
        end
        for i=1:OFDM_symbol_num
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
        X([1:OFDM_symbol_num],[1:OFDM_symbol_length])=ysignal([1:OFDM_symbol_num],[1:OFDM_symbol_length])./H([1:OFDM_symbol_num],[1:OFDM_symbol_length]);
        %case1 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X([1:20],1+12*([1:71]-1))=0;
        X([1:20],1209+12*([1:70]-1))=0;
        X_real = sign(real(X));
        X_imag = sign(imag(X));
        for i = 1:OFDM_symbol_num
            num1(i, :) = (X_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (X_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num2(i, :) = (X_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (X_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num3=sum(sum(num1))+sum(sum(num2));
        error_rate_LINEAR = num3 /(1540*2*20) ;
        errorLINEAR(:,snr)=error_rate_LINEAR;
        %----------------------------------------------------------------
        %case2 DFT   把每一橫行放進DFT
        %生成Hpilot
        for i=1:OFDM_symbol_num
            for a=1:71
                H2(i,1+12*(a-1))=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H2(i,1209+12*(b-1))=ypilot(i,71+b)/pilot(i,71+b);
            end
        end
        for i = 1:OFDM_symbol_num
            D(i,[1:OFDM_symbol_length]) = ifft(H2(i,[1:OFDM_symbol_length]));
        end
        D([1:OFDM_symbol_num],[142:OFDM_symbol_length])=0;
        T([1:OFDM_symbol_num],[1:OFDM_symbol_length])=12.*D;
        for i = 1:OFDM_symbol_num
            Do(i,[1:OFDM_symbol_length]) = fft(T(i,[1:OFDM_symbol_length]));
        end
        X_DFT([1:OFDM_symbol_num],[1:OFDM_symbol_length])=ysignal([1:OFDM_symbol_num],[1:OFDM_symbol_length])./Do([1:OFDM_symbol_num],[1:OFDM_symbol_length]);
        %case2 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_DFT([1:20],1+12*([1:71]-1))=0;
        X_DFT([1:20],1209+12*([1:70]-1))=0;
        DFT_real = sign(real(X_DFT));
        DFT_imag = sign(imag(X_DFT));
        for i = 1:OFDM_symbol_num
            num4(i, :) = (DFT_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (DFT_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num5(i, :) = (DFT_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (DFT_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num6=sum(sum(num4))+sum(sum(num5));
        error_rate_DFT = num6 /(1540*2*20) ;
        errorDFT(:,snr)=error_rate_DFT;
        %----------------------------------------------------------------
        % %case3 LS法
        Le=100;
        f1 = exp(-1i * 2 * pi * (12.*((0:70).')) * (0:Le-1) / 2048);% k =12*( 0:70&71:140 ) ;N = 2048;
        f2 = exp(-1i * 2 * pi * (1209+12.*((0:69).')) * (0:Le-1) / 2048);
        f = [f1;f2];  %f矩陣的設置是N*Le的矩陣 N為領航信號的點數 Le的最佳設置為最大延遲點數(通道點數)
        % L(通道點數) < Le < N(領航信號H的個數)
        %生成Hpilot
        for i=1:OFDM_symbol_num
            for a=1:71
                H3(i,a)=ypilot(i,a)/pilot(i,a);
            end
            for b=1:70
                H3(i,71+b)=ypilot(i,71+b)/pilot(i,71+b);
            end
        end
        h_ls=zeros(Le,OFDM_symbol_num);
        for i=1:OFDM_symbol_num
            hn([1:141],1)=H3(i,[1:141]);
            h_ls([1:Le],i)=inv(f'*f+10^-5.*eye(Le))*f'*hn;
            h_LS=h_ls.';
        end
        h_LS([1:OFDM_symbol_num],[101:OFDM_symbol_length])=0;%估出來的點數為100點 我們補0到2048點 做FFT
        for i = 1:OFDM_symbol_num
            fft_h_LS(i,[1:OFDM_symbol_length]) = fft(h_LS(i,[1:OFDM_symbol_length]));
        end
        X_LS1([1:OFDM_symbol_num],[1:OFDM_symbol_length])=ysignal([1:OFDM_symbol_num],[1:OFDM_symbol_length])./fft_h_LS([1:OFDM_symbol_num],[1:OFDM_symbol_length]);
        %case3 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS1([1:20],1+12*([1:71]-1))=0;
        X_LS1([1:20],1209+12*([1:70]-1))=0;

        h_LS1_real = sign(real(X_LS1));
        h_LS1_imag = sign(imag(X_LS1));
        for i = 1:OFDM_symbol_num
            num7(i, :) = (h_LS1_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (h_LS1_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num8(i, :) = (h_LS1_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (h_LS1_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num9=sum(sum(num7))+sum(sum(num8));
        error_rate_LS1 = num9 /(1540*2*20) ;
        errorLS1(:,snr)=error_rate_LS1;
        %----------------------------------------------------------------
        %case4
        for i=1:OFDM_symbol_num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-5.*eye(100))*x_bar(:,:,i)'*yy;
            hp_ls=hp.';
            %沒有修正參數
            hp0([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i))*x_bar(:,:,i)'*yy;
            hp_ls0=hp0.';

        end
        hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        hp_ls0([1:20],[101:2048])=0;
        for i = 1:OFDM_symbol_num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
            fft_HP_LS0(i,[1:2048]) = fft(hp_ls0(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
        X_LS20([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS0([1:20],[1:2048]);
        %case4 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS2([1:20],1+12*([1:71]-1))=0;
        X_LS2([1:20],1209+12*([1:70]-1))=0;
        X_LS20([1:20],1+12*([1:71]-1))=0;
        X_LS20([1:20],1209+12*([1:70]-1))=0;
        HP_LS2_real = sign(real(X_LS2));
        HP_LS2_imag = sign(imag(X_LS2));
        HP_LS2_real0 = sign(real(X_LS20));
        HP_LS2_imag0 = sign(imag(X_LS20));
        for i = 1:OFDM_symbol_num
            num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
            num15(i, :) = (HP_LS2_real0(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag0(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num16(i, :) = (HP_LS2_real0(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag0(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num12=sum(sum(num10))+sum(sum(num11));
        error_rate_LS2 = num12 /(1540*2*20) ;
        errorLS2(:,snr)=error_rate_LS2;
        num18=sum(sum(num15))+sum(sum(num16));
        error_rate_LS20 = num18 /(1540*2*20) ;
        errorLS20(:,snr)=error_rate_LS20;
        %----------------------------------------------------------------
    end
    %case1
    error1(j,:)=errorLINEAR;
    %case2
    error2(j,:)=errorDFT;
    %case3
    error3(j,:)=errorLS1;
    %case4
    error4(j,:)=errorLS2;
    %case5
    error5(j,:)=errorLS20;

    clc
    done = repmat('■', 1, j);
    nondone = repmat('□', 1, 20 - j);
    fprintf('%s%s\r', done, nondone);
end

%內插法
Linear_e=sum(error1)./20;
%DFT法
DFT_e=sum(error2)./20;
%LS法
LS_e=sum(error3)./20;
%LS法2
LS2_e=sum(error4)./20;
LS20_e=sum(error5)./20;

figure;
semilogy(db,LS2_e, 'b-',db,LS20_e, 'c-');
legend('加上修正參數','沒有加修正參數');
title('時域LS法錯誤率曲線圖');
xlabel('SNR');
ylabel('BER');
grid;

figure;
semilogy(db,Linear_e, 'g-',db,DFT_e, 'k-',db,LS_e,'m-',db,LS2_e, 'b-');
title('錯誤率曲線圖');
legend('線性內插','DFT','LS','時域 LS','Location','southwest');
xlabel('SNR');
ylabel('BER');
grid;

%----------------------------------------------------------------
% 頻譜響應估測圖（預設不執行，需要時將下方程式碼取消註解）
% 藍線：真實通道頻率響應（通道 g 補零至 2048 點後做 FFT）
% 橙線：時域 LS 法估測的頻率響應（第 1 個 OFDM 符元）
%
% 註：迴圈結束後保留的是最後一個通道、SNR = 20 dB 時的結果
%----------------------------------------------------------------
a1=g;
k=1:2048;
a1(2048)=0;
a1f=fft(a1);
figure;
plot(k,abs(a1f),k,abs(H(1,:)));
title('頻譜響應估測圖');
ylabel('頻譜大小響應');
xlabel('第n個子載波');
grid;






