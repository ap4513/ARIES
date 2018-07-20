function [l, I] = LoadOOData(fname)

% started = false;
% N = 1000;
% l = zeros(N, 1);
% I = zeros(N, 1);
% n = 0;
% fid = fopen(fname);
% while ~feof(fid)
% 	data = fgetl(fid);
% 	if started
% 		data = sscanf(data, '%f\t%f');
% 		n = n+1;
% 		if n>N
% 			I = [I; zeros(N, 1)];
% 			l = [l; zeros(N, 1)];
% 			N = 2*N;
% 		end
% 		l(n) = data(1);
% 		I(n) = data(2);
% 	elseif strcmp(data, '>>>>>Begin Spectral Data<<<<<')
% 		started = true;
% 	end
% end
% fclose(fid);
% 
% l = l(1:n);
% I = I(1:n);
temp = textscan(regexp(fileread(fname), ...
	'(?<=>>>>>Begin.*Data<<<<<).*(?=>>>>>End.*Data<<<<<)', 'match', 'once'), ...
	'%f\t%f');
l = temp{1};
I = cell2mat(temp(2:end));
end

