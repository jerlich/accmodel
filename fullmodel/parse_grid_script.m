function [X dX] = parse_grid_script(output_dir, dir_name, output_fname, Nruns)
% function [X dX] = parse_MC_script(output_dir, dir_name, output_fname, Nruns)
% output_dir = '~/output/MC/20000trials/';
% output_fname = 'job_MC_grid_task_%05i.mat';
% Nruns = 4000;

X = zeros(Nruns, 11);
dX = zeros(Nruns, 10);

for i = 1:Nruns,
	try
		load(sprintf([output_dir dir_name '/' output_fname], i));

		X(i, 1) = LL;
		X(i, 2:10) = param1;
		X(i, 11) = i;

		dX(i, 1) = i;
		dX(i, 2:10) = LL_grad;
	catch
		fprintf('run %i could not be loaded or read.\n', i);
	end;
end;