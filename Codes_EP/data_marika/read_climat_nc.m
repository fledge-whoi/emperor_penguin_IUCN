function [year, num, sic]=read_climat_nc(fname)

% output: SIC for 66 colony until 2100
% input: nmae of the ncdf file

year=ncread(fname, 'year');

num=ncread(fname, 'new_num'); %colony id

sic=ncread(fname, 'sic_colony'); %month, year, ens, col
