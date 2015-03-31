function data = read_ansys_output(filelist, varargin)
% produces a .mat file given a list of AQWA .LIS file(s)
% The list can be a 2D char array or a cell array, or it can be a single
% filename

% the output is the hydrodynamic parameters for all 6 DOF stored in a
% structure.

% Version 1.0 - 3/18/2015 - by BL

if iscell(filelist)
    numFiles = length(filelist);
    data = read_file(filelist{1});
    if numFiles > 1
        for ii = 2:numFiles
            dataNew = read_file(filelist{ii});
            fld = fieldnames(dataNew);
            for jj = 1:length(fld)
                if all([~strcmp(fld{jj}, 'K'), ~strcmp(fld{jj}, 'inertia')])
                    data.(fld{jj}) = [data.(fld{jj}) ; ...
                        dataNew.(fld{jj})]; %#ok<*AGROW> SMALL ARRAY, SHORT LOOP
                end
            end
            clear dataNew
        end
    end
elseif ischar(filelist) && numel(filelist) >= length(filelist)
    numFiles = size(filelist,1);
    data = read_file(filelist(1,:));
    if numFiles > 1
        for ii = 2:numFiles
            dataNew = read_file(filelist(ii, :));
            fld = fieldnames(dataNew);
            for jj = 1:length(fld)
                if all([~strcmp(fld{jj}, 'K'), ~strcmp(fld{jj}, 'inertia')])
                    data.(fld{jj}) = [data.(fld{jj}) ; ...
                        dataNew.(fld{jj})]; %#ok<*AGROW> SMALL ARRAY, SHORT LOOP
                end
            end
            clear dataNew
        end
    end
else
    error('Bad Filelist')
end

% look for duplicate frequency data points

uniqueData.inertia = data.inertia;
uniqueData.K = data.K;
[uniqueData.freq, idx, ~] = unique(data.freq, 'rows');

fld = fieldnames(data);
for ii = 1:length(fld)
    if all([~strcmp(fld{ii}, 'K'), ~strcmp(fld{ii}, 'freq'), ~strcmp(fld{ii}, 'inertia')])
        uniqueData.(fld{ii}) = data.(fld{ii})(idx, :);
    end
end

data = uniqueData;

% save to a mat file if given a filename
if nargin > 1
    save(varargin{1}, '-struct', 'data');
end



function [data] = read_file(filename)
fileStr = fileread(filename);

% Get Added Mass
addMassStr = regexp(fileStr, '(ADDED MASS-VAR).*?\r\n1\r\n', 'match');
numbers = textscan(addMassStr{1}, '%n', 'Headerlines', 5);
addedMassData = reshape(numbers{1}(1:end-1), 14 , length(numbers{1}(1:end-1))/14)';

freq = addedMassData(:,2);
A =  addedMassData(:,3:14);
clear numbers addedMassData addMassStr

% Damping Coefficients
dampStr = regexp(fileStr, '(DAMPING-VAR).*?\r\n1\r\n', 'match');
numbers = textscan(dampStr{1}, '%n', 'Headerlines', 5);
dampData = reshape(numbers{1}(1:end-1), 14 , length(numbers{1}(1:end-1))/14)';

C = dampData(:,3:14);
clear numbers dampData dampStr


% Diffraction Forces
diffStr = regexp(fileStr, '(DIFF[A-Z \-\\]*PERIOD).*?\r\n1\r\n', 'match');
numbers = textscan(diffStr{1}, '%n', 'Headerlines', 7);
numbers{1}(3) = [];
numbers{1}(end) = [];
diffData = reshape(numbers{1}, 14, length(numbers{1}) / 14)';

diffMag   = diffData(:,3:2:13);
diffPhase = diffData(:,4:2:14);
clear numbers diffStr diffData
% fid = fopen('test.txt', 'w');
% fprintf(fid, diffStr{1});



% Froude Krylov Forces
fkStr = regexp(fileStr, '(FROUDE KRYLOV FOR[A-Z \-\\]*PERIOD).*?\r\n1\r\n', 'match');
numbers = textscan(fkStr{1}, '%n', 'Headerlines', 7);
numbers{1}(3) = [];
numbers{1}(end) = [];
fkData = reshape(numbers{1}, 14, length(numbers{1}) / 14)';


fkMag   = fkData(:,3:2:13);
fkPhase = fkData(:,4:2:14);
clear numbers fkStr fkData


% Diffraction + Froude Krylov Forces
dfkStr = regexp(fileStr, '(FROUDE KRYLOV \+ DIFF[A-Z \-\\]*PERIOD).*?\r\n1\r\n', 'match');
numbers = textscan(dfkStr{1}, '%n', 'Headerlines', 7);
numbers{1}(3) = [];
numbers{1}(end) = [];
dfkData = reshape(numbers{1}, 14, length(numbers{1}) / 14)';


dfkMag   = dfkData(:,3:2:13);
dfkPhase = dfkData(:,4:2:14);
clear numbers 



% get the stiffness matrix
stiffStr = regexp(fileStr, '[ \t]*STIFFNESS MATRIX[\n\r][0-9a-zA-Z\.\+\-\n\r \t]*[\n\r][1-9][\n\r]', 'match');
stiffStr = regexprep(stiffStr{1}, '\r\n\r\n', '\n');
stiffStr = textscan(stiffStr, '%s%n%n%n%n%n%n', 'HeaderLines', 5);
numbers = cell2mat(stiffStr(2:end));
numbers(end,:) = [];
stiffness = numbers;
clear numbers

% get the mass and inertia matrix.
inertiaStr = regexp(fileStr, 'INERTIA MATRIX[0-9E\+\.\t \-\n\r]*\r\n\r\n', 'match');
temp = regexprep(inertiaStr{1}, '   *','\n');
temp = regexprep(temp, '[\n\r][\n\r]*', '\n');
numbers = textscan(temp, '%n', 'HeaderLines', 1);
inertia = reshape(numbers{1}, 3, 3)';

massStr = regexp(fileStr,'T O T A L[0-9a-zA-Z \t\+\-\.]*[\n\r]', 'match');
massNumbers = textscan(regexprep(massStr{1},'T O T A L[]*', ''), '%n');
mass = massNumbers{1}(2);

inertiaData = [diag([mass mass mass]) zeros(3) ; zeros(3)  inertia];
% Add to output structure

data.freq = freq;
data.inertia = inertiaData;
data.K = stiffness;
data.AddMass = A;
data.RadDamp = C;
data.diffMag = diffMag;
data.diffPhase = diffPhase .* (pi/180);
data.fkMag = fkMag;
data.dfPhase = fkPhase .* (pi/180);
data.dfkMag = dfkMag;
data.dfkPhase = dfkPhase .* (pi/180);



