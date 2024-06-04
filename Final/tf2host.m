function tf2host(myTF, filename)
% function tf2host(myTF, filename)
% this function takes any matlab system block and converts 
% to a regular transfer function [num, den] stuff
% and sends it to a format our Labview Code will Accept
% note: it sends 
%  filter length
%  Ts
%  num
%  den
%
% Written by Kevin Chu
% Modified by Yigang Wang
% August 11, 2009


% just makes sure that there is a file extension
if(length(filename)<5 || filename(end-3)~='.')
    filename = [filename '.dat'];
end

[num,den,Ts] = tfdata(myTF,'v');

relOrder = length(den)-length(num);

num = [zeros(1,relOrder), num];

fid = fopen(filename,'w+');
u = [length(den),Ts,num,den];
fwrite(fid,u,'float64', 'ieee-le');
fclose(fid);

% file2tar(filename,'');

