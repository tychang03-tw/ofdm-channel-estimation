% 傳統 LS 與時域 LS：不同最大通道延遲點數比較
% 衰減係數 α = 0.001，最大通道延遲 100 / 80 / 30 點，兩種方法的修正參數皆為 λ = 10^-5
% 輸出：傳統 LS、時域 LS 各一張 BER vs. SNR（三種延遲點數）
clear
clc
num = 20 ;       %幾組 OFDM symbol
length1 = 2048;  %長度  20*2048矩陣
cp = 128;
xsignal=zeros(num, length1);
c = zeros(1, num * (length1 + cp));

pilot=zeros(num,141);
ypilot=zeros(num,141);
ysignal=zeros(num, length1);
resignal=zeros(num, length1);

H=zeros(num,length1);
X=zeros(num,length1);
num1=zeros(num,841);
num2=zeros(num,840);

H2=zeros(num,length1);
D=zeros(num,length1);
Do=zeros(num,length1);
num4=zeros(num,841);
num5=zeros(num,840);

H3=zeros(num,141);
num7=zeros(num,841);
num8=zeros(num,840);

ifft_xpilot=zeros(num,length1);
ifft_hpilot=zeros(num,length1);
fft_HP_LS=zeros(num,length1);
X_LS2=zeros(num,length1);
xpilot=zeros(num,length1);
hpilot=zeros(num,length1);

num10=zeros(num,841);
num11=zeros(num,840);


for i = 1:num
    x = sign(randn(1,841)) + 1i * sign(randn(1, 841));
    x(1209:2048)=sign(randn(1,840)) + 1i * sign(randn(1, 840));
    xsignal(i,:)=x;

    %----------------------------------------------------------------
    %case4
    for a=1:71
        xpilot(i,1+12*(a-1))=xsignal(i,1+12*(a-1));
    end
    for b=1:70
        xpilot(i,1209+12*(b-1))=xsignal(i,1209+12*(b-1));
    end
    ifft_xpilot(i,1:length1) = ifft(xpilot(i,1:length1));


    for k=0:29
        shift=[0 k];%shift 點數為通道點數
        x_bar3([1:2048],k+1,i) =circshift(ifft_xpilot(i,[1:2048]),shift);
    end
    %----------------------------------------------------------------
    %----------------------------------------------------------------
    %case4
    for a=1:71
        xpilot(i,1+12*(a-1))=xsignal(i,1+12*(a-1));
    end
    for b=1:70
        xpilot(i,1209+12*(b-1))=xsignal(i,1209+12*(b-1));
    end
    ifft_xpilot(i,1:length1) = ifft(xpilot(i,1:length1));


    for k=0:79
        shift=[0 k];%shift 點數為通道點數
        x_bar2([1:2048],k+1,i) =circshift(ifft_xpilot(i,[1:2048]),shift);
    end
    %----------------------------------------------------------------
    %----------------------------------------------------------------
    %case4
    for a=1:71
        xpilot(i,1+12*(a-1))=xsignal(i,1+12*(a-1));
    end
    for b=1:70
        xpilot(i,1209+12*(b-1))=xsignal(i,1209+12*(b-1));
    end
    ifft_xpilot(i,1:length1) = ifft(xpilot(i,1:length1));


    for k=0:99
        shift=[0 k];%shift 點數為通道點數
        x_bar1([1:2048],k+1,i) =circshift(ifft_xpilot(i,[1:2048]),shift);
    end
    %----------------------------------------------------------------

    %pilot
    for a=1:71
        pilot(i,a)=xsignal(i,1+12*(a-1));
    end
    for b=1:70
        pilot(i,1+70+b)=xsignal(i,1209+12*(b-1));
    end

    %ifft
    ifft_xsignal = ifft(xsignal(i,1:length1));

    %cp
    guard=ifft_xsignal((length1 -cp+1) :length1);
    ifft_cpxsignal = [guard,ifft_xsignal];%後16位元加上x64=80位元
    c(1, 1 + 2176*(i-1) : 2176*i) = ifft_cpxsignal;
end


%100點通道
for j=1:20
    n=randn(1,100) + 1i * randn(1, 100);
    g=exp(-0.001*[0:99]).*n;
    C=conv(c,g);
    E=(C*C')/43619;

    db=1:20;
    %20個錯誤率
    for snr=1:length(db)
        noise=randn(1, 43619)+1i.*randn(1, 43619);

        N =(0.5 .* 10^((-db(snr)) / 10))*E;
        Y = C + noise .*sqrt(N);

        %移除保護區間
        for i = 1:num
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
        % %case3 LS法
        Le=100;
        f1 = exp(-1i * 2 * pi * (12.*((0:70).')) * (0:Le-1) / 2048);% k =12*( 0:70&71:140 ) ;N = 2048;
        f2 = exp(-1i * 2 * pi * (1209+12.*((0:69).')) * (0:Le-1) / 2048);
        f = [f1;f2];  %f矩陣的設置是N*Le的矩陣 N為領航信號的點數 Le的最佳設置為最大延遲點數(通道點數)
        % L(通道點數) < Le < N(領航信號H的個數)

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
            h_ls([1:Le],i)=inv(f'*f+10^-5.*eye(Le))*f'*hn;
            h_LS=h_ls.';
        end
        h_LS([1:num],[101:length1])=0;%估出來的點數為100點 我們補0到2048點 做FFT

        for i = 1:num
            fft_h_LS(i,[1:length1]) = fft(h_LS(i,[1:length1]));
        end
        X_LS1([1:num],[1:length1])=ysignal([1:num],[1:length1])./fft_h_LS([1:num],[1:length1]);
        %case3 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS1([1:20],1+12*([1:71]-1))=0;
        X_LS1([1:20],1209+12*([1:70]-1))=0;

        h_LS1_real = sign(real(X_LS1));
        h_LS1_imag = sign(imag(X_LS1));
        for i = 1:num
            num7(i, :) = (h_LS1_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (h_LS1_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num8(i, :) = (h_LS1_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (h_LS1_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num9=sum(sum(num7))+sum(sum(num8));
        error_rate_LS1 = num9 /(1540*2*20) ;
        errorLS1(:,snr)=error_rate_LS1;
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        %case4
        for i=1:num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:100],i)=inv(x_bar1(:,:,i)'*x_bar1(:,:,i)+10^-5.*eye(100))*x_bar1(:,:,i)'*yy;
            hp_ls=hp.';

        end
        hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        
        for i = 1:num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);


        %case4 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS2([1:20],1+12*([1:71]-1))=0;
        X_LS2([1:20],1209+12*([1:70]-1))=0;
        HP_LS2_real = sign(real(X_LS2));
        HP_LS2_imag = sign(imag(X_LS2));
        for i = 1:num
            num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num12=sum(sum(num10))+sum(sum(num11));
        error_rate_LS2 = num12 /(1540*2*20) ;
        errorLS2(:,snr)=error_rate_LS2;
        %----------------------------------------------------------------

    end
    %case3
    error3(j,:)=errorLS1;
    %case4
    error4(j,:)=errorLS2;
end

%LS法
LS_e1=sum(error3)./20;
%LS法2
LS2_e1=sum(error4)./20;

%80點通道
for j=1:20
    n=randn(1,80) + 1i * randn(1, 80);
    g=exp(-0.001*[0:79]).*n;
    C=conv(c,g);
    E=(C*C')/43599;

    db=1:20;
    %20個錯誤率
    for snr=1:length(db)
        noise=randn(1, 43599)+1i.*randn(1, 43599);

        N =(0.5 .* 10^((-db(snr)) / 10))*E;
        Y = C + noise .*sqrt(N);

        %移除保護區間
        for i = 1:num
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
        % %case3 LS法
        Le=80;
        f1 = exp(-1i * 2 * pi * (12.*((0:70).')) * (0:Le-1) / 2048);% k =12*( 0:70&71:140 ) ;N = 2048;
        f2 = exp(-1i * 2 * pi * (1209+12.*((0:69).')) * (0:Le-1) / 2048);
        f = [f1;f2];  %f矩陣的設置是N*Le的矩陣 N為領航信號的點數 Le的最佳設置為最大延遲點數(通道點數)
        % L(通道點數) < Le < N(領航信號H的個數)

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
            h_ls([1:Le],i)=inv(f'*f+10^-5.*eye(Le))*f'*hn;
            h_LS=h_ls.';
        end
        h_LS([1:num],[81:length1])=0;%估出來的點數為80點 我們補0到2048點 做FFT
        for i = 1:num
            fft_h_LS(i,[1:length1]) = fft(h_LS(i,[1:length1]));
        end
        X_LS1([1:num],[1:length1])=ysignal([1:num],[1:length1])./fft_h_LS([1:num],[1:length1]);
        %case3 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS1([1:20],1+12*([1:71]-1))=0;
        X_LS1([1:20],1209+12*([1:70]-1))=0;

        h_LS1_real = sign(real(X_LS1));
        h_LS1_imag = sign(imag(X_LS1));
        for i = 1:num
            num7(i, :) = (h_LS1_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (h_LS1_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num8(i, :) = (h_LS1_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (h_LS1_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num9=sum(sum(num7))+sum(sum(num8));
        error_rate_LS1 = num9 /(1540*2*20) ;
        errorLS1(:,snr)=error_rate_LS1;
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        %case4
        for i=1:num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:80],i)=inv(x_bar2(:,:,i)'*x_bar2(:,:,i)+10^-5.*eye(80))*x_bar2(:,:,i)'*yy;
            hp_ls=hp.';

        end
        hp_ls([1:20],[81:2048])=0;% 估出來的點數為80點 我們補0到2048點 做FFT
        
        for i = 1:num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);


        %case4 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS2([1:20],1+12*([1:71]-1))=0;
        X_LS2([1:20],1209+12*([1:70]-1))=0;
        HP_LS2_real = sign(real(X_LS2));
        HP_LS2_imag = sign(imag(X_LS2));
        for i = 1:num
            num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num12=sum(sum(num10))+sum(sum(num11));
        error_rate_LS2 = num12 /(1540*2*20) ;
        errorLS2(:,snr)=error_rate_LS2;
        %----------------------------------------------------------------

    end
    %case3
    error3(j,:)=errorLS1;
    %case4
    error4(j,:)=errorLS2;
end
%LS法
LS_e2=sum(error3)./20;
%LS法2
LS2_e2=sum(error4)./20;

%30點通道
for j=1:20
    n=randn(1,30) + 1i * randn(1, 30);
    g=exp(-0.001*[0:29]).*n;
    C=conv(c,g);
    E=(C*C')/43549;

    db=1:20;
    %20個錯誤率
    for snr=1:length(db)
        noise=randn(1, 43549)+1i.*randn(1, 43549);

        N =(0.5 .* 10^((-db(snr)) / 10))*E;
        Y = C + noise .*sqrt(N);

        %移除保護區間
        for i = 1:num
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
        % %case3 LS法
        Le=30;
        f1 = exp(-1i * 2 * pi * (12.*((0:70).')) * (0:Le-1) / 2048);% k =12*( 0:70&71:140 ) ;N = 2048;
        f2 = exp(-1i * 2 * pi * (1209+12.*((0:69).')) * (0:Le-1) / 2048);
        f = [f1;f2];  %f矩陣的設置是N*Le的矩陣 N為領航信號的點數 Le的最佳設置為最大延遲點數(通道點數)
        % L(通道點數) < Le < N(領航信號H的個數)

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
            h_ls([1:Le],i)=inv(f'*f+10^-5.*eye(Le))*f'*hn;
            h_LS=h_ls.';
        end
        h_LS([1:num],[31:length1])=0;%估出來的點數為30點 我們補0到2048點 做FFT
        for i = 1:num
            fft_h_LS(i,[1:length1]) = fft(h_LS(i,[1:length1]));
        end
        X_LS1([1:num],[1:length1])=ysignal([1:num],[1:length1])./fft_h_LS([1:num],[1:length1]);
        %case3 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS1([1:20],1+12*([1:71]-1))=0;
        X_LS1([1:20],1209+12*([1:70]-1))=0;

        h_LS1_real = sign(real(X_LS1));
        h_LS1_imag = sign(imag(X_LS1));
        for i = 1:num
            num7(i, :) = (h_LS1_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (h_LS1_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num8(i, :) = (h_LS1_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (h_LS1_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num9=sum(sum(num7))+sum(sum(num8));
        error_rate_LS1 = num9 /(1540*2*20) ;
        errorLS1(:,snr)=error_rate_LS1;
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        %case4
        for i=1:num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:30],i)=inv(x_bar3(:,:,i)'*x_bar3(:,:,i)+10^-5.*eye(30))*x_bar3(:,:,i)'*yy;
            hp_ls=hp.';

        end
        hp_ls([1:20],[31:2048])=0;% 估出來的點數為30點 我們補0到2048點 做FFT
        for i = 1:num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);


        %case4 錯誤率
        xsignal([1:20],1+12*([1:71]-1))=0;
        xsignal([1:20],1209+12*([1:70]-1))=0;
        X_LS2([1:20],1+12*([1:71]-1))=0;
        X_LS2([1:20],1209+12*([1:70]-1))=0;
        HP_LS2_real = sign(real(X_LS2));
        HP_LS2_imag = sign(imag(X_LS2));
        for i = 1:num
            num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
            num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        end
        num12=sum(sum(num10))+sum(sum(num11));
        error_rate_LS2 = num12 /(1540*2*20) ;
        errorLS2(:,snr)=error_rate_LS2;
        %----------------------------------------------------------------

    end
    %case3
    error3(j,:)=errorLS1;
    %case4
    error4(j,:)=errorLS2;
end

%LS法
LS_e3=sum(error3)./20;
%LS法2
LS2_e3=sum(error4)./20;

figure;
semilogy(db,LS_e1,'g-',db,LS_e2,'b-',db,LS_e3,'r-');
title('LS法錯誤率曲線圖');
legend('最大通道延遲點數100點','最大通道延遲點數80點','最大通道延遲點數30點');
xlabel('SNR');
ylabel('BER');
grid;

figure;
semilogy(db,LS2_e1,'g-',db,LS2_e2,'b-',db,LS2_e3,'r-');
title('時域LS法錯誤率曲線圖');
legend('最大通道延遲點數100點','最大通道延遲點數80點','最大通道延遲點數30點');
xlabel('SNR');
ylabel('BER');
grid;






