function [hydParms] = read_hyd_file(filename)
% reads a .hyd file and outputs its contents
%
% [hydParms] = read_hyd_file(filename)

fid = fopen(filename,'r');
headers = textscan(fid,'%s',9,'delimiter',',');
data = textscan(fid,'%f%f%f%f%f%f%f%f%f','delimiter',',');
fclose(fid);


% Only take results less than 2 rads/s
%idx = find(data{1} <= 2);
idx = 1:length(data{1});

hydParms.freq      = data{1}(idx);  % rad/s
hydParms.diffMag   = data{2}(idx);
hydParms.diffPhase = data{3}(idx) .* (pi/180); % deg to rad
hydParms.fkMag     = data{4}(idx);
hydParms.fkPhase   = data{5}(idx) .* (pi/180); % deg to rad
hydParms.dfkMag    = data{6}(idx);
hydParms.dfkPhase  = data{7}(idx) .* (pi/180); % deg to rad
hydParms.radDamp   = data{8}(idx);
hydParms.addMass   = data{9}(idx);






