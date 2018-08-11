% This script will read in an Oligomer*.txt file,
% analyze all of the configurations in that file and
% create separate numerical and graphical population
% maps for a user-specified oligomer size.

oligorows = 21;                                             % specify matrix dimensions
oligocols = 13;

j = input('Size of oligomer: ');                            % ask user for size of oligomer to investigate

fid = fopen( ['Oligomer' num2str(j) '.txt'] );              % open infile

if fid > 0                                                  % make sure the file is not empty
    block = 1;                                              % start counter to keep track of configurations
    while ~feof(fid)                                        % read infile
        mystruct(block).title = fscanf(fid, '%s', 1);
        mystruct(block).oligo  = ...                        % fscanf fills array in column order, so transpose results
          fscanf(fid, '%f', [oligocols, oligorows])';
        block = block + 1;                                  % add one to counter; loop until end of infile
    end
else                                                        % terminate if oligomer file doesn't exist
    fprintf('\n\nFile Oligomer%d.txt not found.\n\n',j)
    return
end

fclose(fid);                                                % close the infile

M = zeros(oligorows,oligocols);                             % initialization of final matrix

% start matrix manipulations:
for x=1:block-2
   mystruct(x).oligo = abs(mystruct(x).oligo);              % take abs to accomodate for hydrolyzed sites
   M = M + mystruct(x).oligo;                               % add up all matrices
end

diary ('OligoAvgTotal.txt')                                 % write corresponding outfile
j
M = M/(block-2)                                             % take matrix average
%bar3 (M)                                                    % produce population bar graph and save
%saveas (gcf,['OligoAvg' num2str(j) '.fig']);
diary off                                                   % close outfile
HeatMap(M)
