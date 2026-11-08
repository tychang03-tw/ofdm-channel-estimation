% 時域 LS 法：不同修正參數比較
% 最大通道延遲 100 點、α = 0.001；程式含 λ = 10^-1～10^-5 與無修正參數共 6 組，預設只執行 10^-5 與無修正參數，其他參數需手動取消迴圈內對應區塊與檔尾畫圖段落的註解
% 輸出：BER vs. SNR
clear

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

% 生成QPSK調製的OFDM符號和插入pilot序列
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


    for k=0:99
        shift=[0 k];%shift 點數為通道點數
        x_bar([1:2048],k+1,i) =circshift(ifft_xpilot(i,[1:2048]),shift);
    end
    %----------------------------------------------------------------

    %ifft
    ifft_xsignal = ifft(xsignal(i,1:length1));

    %cp
    guard=ifft_xsignal((length1 -cp+1) :length1);
    ifft_cpxsignal = [guard,ifft_xsignal];%後16位元加上x64=80位元
    c(1, 1 + 2176*(i-1) : 2176*i) = ifft_cpxsignal;
end


%20個通道
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
        end
        %----------------------------------------------------------------
        % %10^-1
        % for i=1:num
        %     ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
        %     yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
        %     %有修正參數
        %     hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-1.*eye(100))*x_bar(:,:,i)'*yy;
        %     hp_ls=hp.';
        % end
        % hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        % for i = 1:num
        %     fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        % end
        % X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
        % xsignal([1:20],1+12*([1:71]-1))=0;
        % xsignal([1:20],1209+12*([1:70]-1))=0;
        % X_LS2([1:20],1+12*([1:71]-1))=0;
        % X_LS2([1:20],1209+12*([1:70]-1))=0;
        % HP_LS2_real = sign(real(X_LS2));
        % HP_LS2_imag = sign(imag(X_LS2));
        % for i = 1:num
        %     num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
        %     num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        % end
        % num12=sum(sum(num10))+sum(sum(num11));
        % error_rate_LS2 = num12 /(1540*2*20) ;
        % errorLS2(:,snr)=error_rate_LS2;
        % %----------------------------------------------------------------
        % %----------------------------------------------------------------
        % %10^-2
        % for i=1:num
        %     ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
        %     yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
        %     %有修正參數
        %     hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-2.*eye(100))*x_bar(:,:,i)'*yy;
        %     hp_ls=hp.';
        % end
        % hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        % for i = 1:num
        %     fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        % end
        % X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
        % xsignal([1:20],1+12*([1:71]-1))=0;
        % xsignal([1:20],1209+12*([1:70]-1))=0;
        % X_LS2([1:20],1+12*([1:71]-1))=0;
        % X_LS2([1:20],1209+12*([1:70]-1))=0;
        % HP_LS2_real = sign(real(X_LS2));
        % HP_LS2_imag = sign(imag(X_LS2));
        % for i = 1:num
        %     num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
        %     num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        % end
        % num12=sum(sum(num10))+sum(sum(num11));
        % error_rate_LS2 = num12 /(1540*2*20) ;
        % errorLS3(:,snr)=error_rate_LS2;
        % %----------------------------------------------------------------
        % %----------------------------------------------------------------
        % %10^-3
        % for i=1:num
        %     ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
        %     yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
        %     %有修正參數
        %     hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-3.*eye(100))*x_bar(:,:,i)'*yy;
        %     hp_ls=hp.';
        % end
        % hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        % for i = 1:num
        %     fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        % end
        % X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
        % xsignal([1:20],1+12*([1:71]-1))=0;
        % xsignal([1:20],1209+12*([1:70]-1))=0;
        % X_LS2([1:20],1+12*([1:71]-1))=0;
        % X_LS2([1:20],1209+12*([1:70]-1))=0;
        % HP_LS2_real = sign(real(X_LS2));
        % HP_LS2_imag = sign(imag(X_LS2));
        % for i = 1:num
        %     num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
        %     num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        % end
        % num12=sum(sum(num10))+sum(sum(num11));
        % error_rate_LS2 = num12 /(1540*2*20) ;
        % errorLS4(:,snr)=error_rate_LS2;
        % %----------------------------------------------------------------
        % %----------------------------------------------------------------
        % %10^-4
        % for i=1:num
        %     ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
        %     yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
        %     %有修正參數
        %     hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-4.*eye(100))*x_bar(:,:,i)'*yy;
        %     hp_ls=hp.';
        % end
        % hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        % for i = 1:num
        %     fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        % end
        % X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
        % xsignal([1:20],1+12*([1:71]-1))=0;
        % xsignal([1:20],1209+12*([1:70]-1))=0;
        % X_LS2([1:20],1+12*([1:71]-1))=0;
        % X_LS2([1:20],1209+12*([1:70]-1))=0;
        % HP_LS2_real = sign(real(X_LS2));
        % HP_LS2_imag = sign(imag(X_LS2));
        % for i = 1:num
        %     num10(i, :) = (HP_LS2_real(i, 1:841) ~= real(xsignal(i, 1:841))) + (HP_LS2_imag(i, 1:841) ~= imag(xsignal(i, 1:841)));
        %     num11(i, :) = (HP_LS2_real(i, 1209:2048) ~= real(xsignal(i, 1209:2048))) + (HP_LS2_imag(i, 1209:2048) ~= imag(xsignal(i, 1209:2048)));
        % end
        % num12=sum(sum(num10))+sum(sum(num11));
        % error_rate_LS2 = num12 /(1540*2*20) ;
        % errorLS5(:,snr)=error_rate_LS2;
        % %----------------------------------------------------------------
        %----------------------------------------------------------------
        %10^-5
        for i=1:num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i)+10^-5.*eye(100))*x_bar(:,:,i)'*yy;
            hp_ls=hp.';
        end
        hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        for i = 1:num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
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
        errorLS6(:,snr)=error_rate_LS2;
        %----------------------------------------------------------------
        %----------------------------------------------------------------
        %沒有修正參數
        for i=1:num
            ifft_ysignal(i,[1:2048])=ifft(ysignal(i,[1:2048]));
            yy([1:2048],1)=ifft_ysignal(i,[1:2048]); %出通道後的所有信號(時域)
            %有修正參數
            hp([1:100],i)=inv(x_bar(:,:,i)'*x_bar(:,:,i))*x_bar(:,:,i)'*yy;
            hp_ls=hp.';
        end
        hp_ls([1:20],[101:2048])=0;% 估出來的點數為100點 我們補0到2048點 做FFT
        for i = 1:num
            fft_HP_LS(i,[1:2048]) = fft(hp_ls(i,[1:2048]));
        end
        X_LS2([1:20],[1:2048])=ysignal([1:20],[1:2048])./fft_HP_LS([1:20],[1:2048]);
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
        errorLS7(:,snr)=error_rate_LS2;
        %----------------------------------------------------------------
    end
    % error4(j,:)=errorLS2;
    % error5(j,:)=errorLS3;
    % error6(j,:)=errorLS4;
    % error7(j,:)=errorLS5;
    error8(j,:)=errorLS6;
    error9(j,:)=errorLS7;
end
%把每個通道的錯誤率取平均
%LS法2

% LS2_e=sum(error4)./20;
% figure;
% semilogy(db,LS2_e,'g-');
% title('錯誤率曲線圖');
% legend('時域LS法-修正參數10^-^1');
% xlabel('(dB)');
% ylabel('(BER)');
% grid;
% 
% LS2_e=sum(error5)./20;
% figure;
% semilogy(db,LS2_e,'b-');
% title('錯誤率曲線圖');
% legend('時域LS法-修正參數10^-^2');
% xlabel('(dB)');
% ylabel('(BER)');
% grid;
% 
% LS2_e=sum(error6)./20;
% figure;
% semilogy(db,LS2_e,'r-');
% title('錯誤率曲線圖');
% legend('時域LS法-修正參數10^-^3');
% xlabel('(dB)');
% ylabel('(BER)');
% grid;
% 
% LS2_e=sum(error7)./20;
% figure;
% semilogy(db,LS2_e,'k-');
% title('錯誤率曲線圖');
% legend('時域LS法-修正參數10^-^4');
% xlabel('(dB)');
% ylabel('(BER)');
% grid;

LS2_e1=sum(error8)./20;
LS2_e2=sum(error9)./20;
figure;
semilogy(db,LS2_e1,'m-',db,LS2_e2,'c-');
title('修正參數有無之錯誤率曲線圖');
legend('修正參數10^-^5','無修正參數','location','southwest');
xlabel('(dB)');
ylabel('(BER)');
grid;

% LS2_e=sum(error9)./20;
% figure;
% semilogy(db,LS2_e,'c-');
% title('錯誤率曲線圖');
% legend('時域LS法-修正參數10^-^6');
% xlabel('(dB)');
% ylabel('(BER)');
% grid;










