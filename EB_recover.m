function [mu,p_alpha,p_beta] = EB_recover(target,base_mat,delta_mu,orth_flag)
% EB_recover - recover data using Empirical Bayesian (type 2 maximum likelihood)
%
%
% Syntax:  [mu,alpha,beta] = EB_recover(target,base_mat)
%
% Inputs:
%    targert - compressed data. M*1
%    base_mat - base matrix. The bases used to fit the target, M*N
%    delta_mu - stop criteria
%    orth_flag - true if base_mat is orthogonal matrix; false if not. 
%                Using orthogonal matrix will reduce computation complexity
%
% Outputs:
%    mu - expectation of the coefficients
%    p_alpha - precision of coefficients
%    p_beta - precision of noise in data
%    
% Example: 
%    [mu,p_alpha,p_beta] = EB_recover(target,base_mat)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: EMG_CS_decode

% Author: Chaohua Wu
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Jan 13 2014; Last revision: 

%------------- BEGIN CODE --------------
% pathor = pwd;
% cd /home/chaohua/Documents/CCA_for_ERP/cvx
% cvx_setup
% cd(pathor)


if nargin < 2
	error('at least 2 parameters!!');
end

if nargin == 2
	orth_flag = false;
	delta_mu = 10^-5;
end

if nargin ==3
	orth_flag = false;
end

if nargin > 3
	error('too many parameters');
end

M = length(target);

[M1, N] = size(base_mat);

if M~=M1
	error('dimension do not match');
end

%%%% initialize
p_alpha = 1;
p_beta = 1;

d_mu = 1;
iter = 0;
mu_pre = 3*ones(N,1);
if orth_flag == true
   
   %Lambda = (p_alpha+p_beta)*eye(N);
   while (d_mu > delta_mu) && (iter < 500)
   
   iter = iter+1;

   mu = p_beta*((p_alpha+p_beta)^-1*eye(N))*base_mat'*target;

   p_gamma = N*p_beta/(p_beta+p_alpha);

   p_alpha = p_gamma/(mu'*mu);

   p_beta = (M - p_gamma)*((target-base_mat*mu)'*(target-base_mat*mu))^-1;

   d_mu = norm(mu-mu_pre,2)/norm(mu_pre,2);

   mu_pre = mu;

   end




else

   while (d_mu > delta_mu) && (iter < 50)
   
   iter = iter+1;

   pb = p_beta*base_mat'*base_mat;

   Lambda = p_alpha*eye(N) + pb;

   mu = p_beta*Lambda^-1*base_mat'*target;

   eig_pb = eig(pb);

   p_gamma = sum(eig_pb./(eig_pb+p_alpha));

   p_alpha = p_gamma/(mu'*mu);

   p_beta = (M - p_gamma)*((target-base_mat*mu)'*(target-base_mat*mu))^-1;

   d_mu = norm(mu-mu_pre,2)/norm(mu_pre,2);

   mu_pre = mu;

   end


end






