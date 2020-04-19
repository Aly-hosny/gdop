
% TX = [0,0;0,3;3,3;3,0;0,0]
% rx = [0,0;3,3;2,1;1,2]
function [GDOP,tx] = gdop(TX,rx,step)
    maxDim = max(TX);
    minDim = min(TX);
    txX = minDim(1):step:maxDim(1);
    txY = minDim(2):step:maxDim(2);
    [txx,txy] = meshgrid(txX,txY);
    tx = [reshape(txx,[],1),reshape(txy,[],1)];
    G_rxa=20;  %overall receiver gain
    P_tx= 80;  %transmit power
    G_tx= 20 ;  %trasnmitter gain
    freq=40;   %Transmit freq in GHz
    % ranges_rx = coder.nullcopy(zeros(size(tx,1),size(rx,1)));
    % GDOP =(zeros(size(tx,1),1));
    N_rx=-70;

    ranges_rx = pdist2(tx,rx,'euclidean');
    pathLoss = 20*log(((4*3.14*ranges_rx)/(0.3/freq))+ 1e-16);
    P_rx=P_tx+G_rxa+G_tx-pathLoss;
    SNR=P_rx-N_rx;
    sigma = (2*exp((-2/15)*SNR));
    
    for x = 1:size(tx,1) 
        H = ((tx(x,:) - rx(2:end,:))./ranges_rx(x,2:end)') -((tx(x,:) - rx(1,:))/ranges_rx(x,1));
        Q1=sigma(x,1)*ones(size(rx,1)-1,size(rx,1)-1);
        Q2=diag(sigma(x,2:end));
        Q =(Q1+Q2);
        B=pinv(H'*pinv(Q)*H );
        S=svd(B);
        GDOP(x)=sqrt(sum(S));

    end
end