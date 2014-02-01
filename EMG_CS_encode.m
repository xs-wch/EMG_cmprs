function [CS_data,sensematrix] = EMG_CS_encode(data,frame_len,frame_outlen)
% EMG_CS_encode - compressing EMG using CS
%
%
% Syntax:  [CS_data,sensematrix] = EMG_CS_encode(data,frame_len,CS_ratio)
%
% Inputs:
%    data - EMG data, one channel, length should be K*2^n
%    frame_len - length of each frame to compress, length must be 2^n
%    CS_ratio - compress ratio
%    
%
% Outputs:
%    CS_data - compressed data
%    sensematrix - sensing matrix (here use gaussian distributed random matrix)
% Example: 
%    [CS_data,sensematrix] = EMG_CS_encode(data,frame_len,CS_ratio)
%
% Other m-files required: 
% Subfunctions: none
% MAT-files required: none
%
% See also: EMG_CS_decode

% Author: Chaohua Wu
% Department of Biomedical Engineering, Tsinghua University
% email: xs.wuchaohua@outlook.com
%
% Jan 3 2014; Last revision: 

%------------- BEGIN CODE --------------

if nargin == 1
	frame_len = 128;
	frame_outlen = frame_len/2;
end 

if nargin == 2
	frame_outlen = floor(frame_len/2);
end

if (nargin < 1)||(nargin > 3)
	error('1-3 parameters are required');
end

%frame_outlen = floor(frame_len/CS_ratio);
frame_num = length(data)/frame_len;
rng(75)
sensematrix = randn(frame_outlen,frame_len);

CS_data = zeros(frame_num*frame_outlen,1);

for i = 1:frame_num
	CS_data((i-1)*frame_outlen +1 : i*frame_outlen) = sensematrix*data((i-1)*frame_len +1 : i*frame_len);
end


