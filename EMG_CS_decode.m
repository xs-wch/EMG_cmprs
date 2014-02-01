function recover_data = EMG_CS_decode(CS_data,sensematrix,base_matrix)
% EMG_CS_decode - recover data from compressed result
%
%
% Syntax:  recover_data = EMG_CS_decode(CS_data,sensematrix,base_matrix)
%
% Inputs:
%    CS_data - compressed data
%    sensematrix - sensing matrix
%    base_matrix - base matrix. The signal can be expressed sparsely with this base. 
%                  Wavelet base matrix will be used here 
%    
%
% Outputs:
%    recover_data - data recoverd
%    
% Example: 
%    recover_data = EMG_CS_decode(CS_data,sensematrix,base_matrix)
%
% Other m-files required: cvx toolbox
% Subfunctions: none
% MAT-files required: none
%
% See also: EMG_CS_encode

% Author: Chaohua Wu
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Jan 12 2014; using lasso embeded in MATLAB
% Jan 3 2014; Last revision: 

%------------- BEGIN CODE --------------
% pathor = pwd;
% cd /home/chaohua/Documents/CCA_for_ERP/cvx
% cvx_setup
% cd(pathor)


if nargin ~= 3
	error('3 parameters are required');
end

[frame_outlen, frame_len] = size(sensematrix);

wcoef_len = size(base_matrix,2);

frame_num = length(CS_data)/frame_outlen;

recover_data = zeros(frame_num*frame_len,1);

matrix_de = sensematrix*base_matrix;

%options.RegX = 10;
%options.RegType = 1;
for i = 1:frame_num
    
	temp_CSdata = CS_data((i-1)*frame_outlen+1:i*frame_outlen);
    
    tempW = lasso(matrix_de,temp_CSdata,'Lambda',43);
    %%%% 256->128 lambda = 6-10
    %%%% 256-64   lambda = 43-52
    
%   tempW = LS(matrix_de',temp_CSdata',options); 
% 	cvx_begin 
%  		variable temp(wcoef_len)
%  		minimize(norm(temp,1))
%  		subject to
%  		   matrix_de*temp == temp_CSdata
% 	cvx_end

	recover_data((i-1)*frame_len+1:i*frame_len) = base_matrix*tempW;
end






