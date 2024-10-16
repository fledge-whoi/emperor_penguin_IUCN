function [year, temp]=read_temp_nc(fname);

% output: SIC for 66 colony until 2100
% input: name of the ncdf file

year=ncread(fname, 'year');

temp=ncread(fname, 'TREFHT_global'); %month, year, ens