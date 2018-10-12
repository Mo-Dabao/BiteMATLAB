% 将.nc文件ncpath中的var写入.tif的文件tifpath中

% ncpath - .nc文件的路径
% varname - .nc中要写入的变量名，.nc文件中要有对应的lon和lon
% tifpath - 目标.tif文件的路径，要求其上一级路径存在

% 墨大宝
% 2018/10/12

function netcdf2geotif(ncpath, varname, tifpath)
    disp(['开始处理 ', ncpath])
    lon = double(ncread(ncpath, 'lon'));
    lat = double(ncread(ncpath, 'lat'));
    variable = ncread(ncpath, varname);
    georef = georefcells(lat([1, end]), lon([1, end]), size(variable));
    geotiffwrite(tifpath, variable, georef)
    disp(['已经生成 ', tifpath])
end
