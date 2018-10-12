% ��.nc�ļ�ncpath�е�varд��.tif���ļ�tifpath��

% ncpath - .nc�ļ���·��
% varname - .nc��Ҫд��ı�������.nc�ļ���Ҫ�ж�Ӧ��lon��lon
% tifpath - Ŀ��.tif�ļ���·����Ҫ������һ��·������

% ī��
% 2018/10/12

function netcdf2geotif(ncpath, varname, tifpath)
    disp(['��ʼ���� ', ncpath])
    lon = double(ncread(ncpath, 'lon'));
    lat = double(ncread(ncpath, 'lat'));
    variable = ncread(ncpath, varname);
    georef = georefcells(lat([1, end]), lon([1, end]), size(variable));
    geotiffwrite(tifpath, variable, georef)
    disp(['�Ѿ����� ', tifpath])
end
