% 線性內插法與 DFT 法：不同最大通道延遲點數比較
% 衰減係數 α = 0.001，最大通道延遲 100 / 80 / 30 點，DFT 法時域保留前 141 點
% 輸出：線性內插法、DFT 法各一張 BER vs. SNR（三種延遲點數）
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
%20個通道100點
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
       
    end
    %case1
    error1(j,:)=errorLINEAR;
    %case2
    error2(j,:)=errorDFT;
    clc
    done = repmat('■', 1, j);
    nondone = repmat('□', 1, 20 - j);
    fprintf('%s%s\r', done, nondone);
end

%內插法
Linear_e1=sum(error1)./20;
%DFT法
DFT_e1=sum(error2)./20;





%20個通道80點
for j = 1:20
    n = randn(1,80) + 1i * randn(1, 80);
    g = exp(-0.001*[0:79]).*n;
    C = conv(c,g);
    E = (C*C')/43599;
    db = 1:20;
    %20個錯誤率
    for snr = 1:length(db)
        noise = randn(1, 43599)+1i.*randn(1, 43599);
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
       
    end
    %case1
    error1(j,:)=errorLINEAR;
    %case2
    error2(j,:)=errorDFT;
    clc
    done = repmat('■', 1, j);
    nondone = repmat('□', 1, 20 - j);
    fprintf('%s%s\r', done, nondone);
end

%內插法
Linear_e2=sum(error1)./20;
%DFT法
DFT_e2=sum(error2)./20;




%20個通道30點
for j = 1:20
    n = randn(1,30) + 1i * randn(1, 30);
    g = exp(-0.001*[0:29]).*n;
    C = conv(c,g);
    E = (C*C')/43549;
    db = 1:20;
    %20個錯誤率
    for snr = 1:length(db)
        noise = randn(1, 43549)+1i.*randn(1, 43549);
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
       
    end
    %case1
    error1(j,:)=errorLINEAR;
    %case2
    error2(j,:)=errorDFT;
    clc
    done = repmat('■', 1, j);
    nondone = repmat('□', 1, 20 - j);
    fprintf('%s%s\r', done, nondone);
end

%內插法
Linear_e3=sum(error1)./20;
%DFT法
DFT_e3=sum(error2)./20;



figure;
semilogy(db,Linear_e1,'g-',db,Linear_e2,'b-',db,Linear_e3,'r-');
title('線性內插法錯誤率曲線圖');
legend('最大通道延遲點數100點','最大通道延遲點數80點','最大通道延遲點數30點');
xlabel('SNR');
ylabel('BER');
grid;

figure;
semilogy(db,DFT_e1,'g-',db,DFT_e2,'b-',db,DFT_e3,'r-');
title('DFT法錯誤率曲線圖');
legend('最大通道延遲點數100點','最大通道延遲點數80點','最大通道延遲點數30點');
xlabel('SNR');
ylabel('BER');
grid;
