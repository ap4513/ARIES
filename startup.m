fprintf('QPI Modelliing library\n');
arrayfun(@(n) fprintf('='), 1:80);
fprintf('\n');
pth = fullfile(pwd, 'lib');
disp(['Path added: << ' pth ' >>']);
addpath(add_to_path(pth));
