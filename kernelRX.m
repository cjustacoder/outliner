%% kernel RX
% global or local mean and variance
% color space in rgb lab 
% for local different outer window size, inner window size
clear all;
close all;
%% read in the image
im = imread('test1.jpg');
im = im2double(im);
[H,L,~] = size(im);

%% hardcoding ground true
gnd_true =false(H, L);
gnd_true(382:386,819:823) = 1;
gnd_true(917:921,862:866) = 1;
gnd_true(919:923,1555:1559) = 1;
figure('Name','groundtrue')
imshow(gnd_true);
%% global rgb, cannot do global for KRX since K is too large. 
% global mean and global variance
% set threshold
% thrh1 = 30; % 30 for test1
% rx_rgb = global_RX(im);
% score_1 = rx_rgb/max(rx_rgb,[],'all');
% [X1,Y1,T1,AUC1] = perfcurve(gnd_true(:),score_1(:),1);
% label1 = rx_rgb > thrh1;
% figure('Name', 'GRX-RGB')
% imshow(label1)
%% apply markov (to be done)
% n_mask = ones(3);
% n_mask(2,2) = 0; % create a mask to perform Ising potential
% label_1 = cell(1,7);
% label_1{1} = label1;
% gama = 0.5*10^(-4); % 1*10-4
% T = 0.75;
% for i =1:6
%     anom_normal_cnt = -8+2.*conv2(double(label_1{i}),n_mask,'same'); % Qf-Qb 
%     % count the number of anomalous 
%     if i ==1 %% for parameter gama adjustment
%         P_t = gama.*exp((1/T).*anom_normal_cnt);
%     end
%     label_1{i+1} = P<gama.*exp((1/T).*anom_normal_cnt); % improb is the things unchange
% end
% figure('Name','GRX with MRF')
% imshow(label_1{7})
%% global lab
% im_lab = rgb2lab(im);
% rx_lab = global_RX(im_lab);
% thrh2 = 30; % same as rgb? 20 no FN but much FP, 30 less FP, has FN
% score_2 = rx_lab/max(rx_lab(:));
% [X2,Y2,T2,AUC2] = perfcurve(gnd_true(:),score_2(:),1);
% label2 = rx_lab > thrh2;
% figure('Name', 'GRX-Lab')
% imshow(label2)
%% local rgb 7*7 15*15 kernel 
% fail nonsense output and extremely long runtime
% outerw = 7;
% innerw = 3;
% c = 3; % need to be adjusted
% % number of background pixel
% p_num = (outerw*2+1)^2-(innerw*2+1)^2;
% [H,L,~] = size(im);
% v = zeros(H,L);
% for i = 1+outerw: H-outerw
%     for j = 1+outerw: L-outerw
%         r = im(i,j,:);
%         r = reshape(r,3,1);
%         
%         im_copy = im(i-outerw:i+outerw,j-outerw:j+outerw,:); % outer window
% %         %temp = []; % remove inner window
% %         temp = im_copy(1:outerw-innerw,:,:);
% %         temp = cat(2,temp,im_copy(2+innerw+outerw:end,:,:));
% %         temp2 = im_copy(outerw-innerw+1:innerw+outerw+1,1:outerw-innerw,:);
% %         temp2 = cat(2,temp2,im_copy(outerw-innerw+1:innerw+outerw+1,outerw+innerw+2:end,:));
%         im_copy(i-innerw:i+innerw,j-innerw:j+innerw,:) = NaN; % clean inner window cannot
%         % delete so set to NaN
%         % construct K matrix
%         temp = reshape(im_copy,[],size(im_copy,3),1);
%         temp = temp(all(~isnan(temp),2),:); % for nan - rows
% %         B = reshape(A,[],size(A,2),1);
%         K = zeros(p_num,p_num);
% %         temp = reshape(temp,[],size(temp,3),1);
%         for p = 1: p_num
%             for q = 1: p_num
%                 K(p,q) = exp(-(norm(temp(p,:)-temp(q,:)))^2/c);
%             end
%         end
%         % go on on gama
%         gama = zeros(p_num,1);
%         for p = 1:p_num
%             gama(p) = exp(-(norm(temp(p,:)-r))^2/c);
%         end
%         w = ones(p_num,1)*1/p_num; % step 2
%         mu = K*w;
%         e = ones(p_num,1);
%         K_hat = K-mu*e'-e*mu'+e*w'*mu*e';
%         [V,D,U] = svd(K_hat); % pca
%         for p = 1:p_num
%             if (D(p,p) < 10e-8*D(1,1))
% %                 record = p;
%                 break;
%             end
%         end
%         D_bar = D(1:p-1,1:p-1);
%         V_bar = V(:,1:p-1);
%         mu_bar = mu - e*w'*mu;
%         gama_bar = gama - e*w'*gama;
%         v(i,j) = norm(D_bar\V_bar'*(gama_bar-mu_bar));
%     end
% end
%% local normal RX
%% local rgb 7*7 15*15 no kernel 
% extremely long runtime
outerw = 7;
innerw = 3;
% number of background pixel
p_num = (outerw*2+1)^2-(innerw*2+1)^2;
[H,L,~] = size(im);
rx = zeros(H,L);
parfor i = 1+outerw: H-outerw
    for j = 1+outerw: L-outerw
        r = im(i,j,:);
        r = reshape(r,3,1);      
        im_copy = im(i-outerw:i+outerw,j-outerw:j+outerw,:); % outer window
        im_copy(i-innerw:i+innerw,j-innerw:j+innerw,:) = NaN; % clean inner window cannot
        % delete so set to NaN
        temp = reshape(im_copy,[],size(im_copy,3),1);
        temp = temp(all(~isnan(temp),2),:); % for nan - rows
        % calculate mean of outer window
        out_mean = mean(temp);
        % calculate covariance of outer window
        out_cov = cov(temp);
        rx(i,j) = ((r - out_mean')'/out_cov)*(r - out_mean');
    end
end

%% local rgb 7*7 15*15 kernel 
% outerw = 7;
% innerw = 3;
% c = 3; % need to be adjusted
% % number of background pixel
% p_num = (outerw*2+1)^2-(innerw*2+1)^2;
% [H,L,~] = size(im);
% v = zeros(H,L);
% parfor i = 1+outerw: H-outerw
%     for j = 1+outerw: L-outerw
%         r = im(i,j,:);
%         r = reshape(r,3,1);
%         
%         im_copy = im(i-outerw:i+outerw,j-outerw:j+outerw,:); % outer window
%         %temp = []; % remove inner window
%         temp = im_copy(1:outerw-innerw,:,:);
%         temp = cat(2,temp,im_copy(2+innerw+outerw:end,:,:));
%         temp2 = im_copy(outerw-innerw+1:innerw+outerw+1,1:outerw-innerw,:);
%         temp2 = cat(2,temp2,im_copy(outerw-innerw+1:innerw+outerw+1,outerw+innerw+2:end,:));
% %         im_copy(i-innerw:i+innerw,j-innerw:j+innerw,:) = NaN; % clean inner window cannot
%         % delete so set to NaN
%         % construct K matrix
%         temp = reshape(temp,[],size(temp,3),1);
%         temp2 = reshape(temp2,[],size(temp2,3),1);
%         temp = cat(1, temp, temp2);
%         K = zeros(p_num,p_num);
% %         temp = reshape(temp,[],size(temp,3),1);
%         [m,~] = size(temp);
%         for p = 1: m
%             for q = 1: m
%                 K(p,q) = exp(-(norm(temp(p,:)-temp(q,:)))^2/c);
%             end
%         end
%         gama = zeros(m,1);
%         for p = 1:m
%             gama(p) = exp(-(norm(temp(p,:)-r))^2/c);
%         end
%         w = ones(p_num,1)*1/p_num;
%         mu = K*w;
%         e = ones(p_num,1);
%         K_hat = K-mu*e'-e*mu'+e*w'*mu*e';
%         [U,S,V] = svd(K_hat); % pca
%         [m,~] = size(S);
%         for p = 1:m
%             if (S(p,p) < 10e-8*S(1,1))
% %                 record = p;
%                 break;
%             end
%         end
%         S_bar = S(1:p-1,1:p-1);
%         V_bar = U(:,1:p-1);
%         mu_bar = mu - e*w'*mu;
%         gama_bar = gama - e*w'*gama;
%         v(i,j) = norm(S_bar\V_bar'*(gama_bar-mu_bar));
%     end
% end
        
        

%% local rgb 7*7 19*19
%% local rgb 7*7 23*23
%% local lab 7*7 15*15
%% local lab 7*7 19*19
%% local lab 7*7 23*23

%% function area
function rx = global_RX(im)
[H,L,~] = size(im);
temp = reshape(im,[],size(im,3),1);
C = cov(temp);
im_mean = mean(temp);
rx = zeros(H,L);
parfor i = 1: H
    for j = 1:L
        temp = reshape(im(i,j,:),3,1);
        rx(i,j) = ((temp - im_mean)'/C)*(temp - im_mean);
    end
end
end

% function rx = local_RX(im,innerw, outerw)
% [H,L,~] = size(im);
% rx = zeros(H,L);
% parfor i = 1: H
%     for j = 1: L
%         
% end
