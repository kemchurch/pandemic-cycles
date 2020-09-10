function DF_alpha = function_DF_alpha_int(SI,alpha,para,N)
% Set para(7)=para(8)=0 so we can recycle the old cheb-threshold function
para(7)=0; para(8)=0;
% Preprocessing
[~,I]=convert_1to2_int(SI,N);
I = reshape(I,[N+1,4]);
% Output
DF_alpha = [intval(zeros(8*(N+1),1));...
    cheb_threshold(I(:,2),para,alpha,1,N);...
    cheb_threshold(I(:,4),para,alpha,3,N)];
end