function initialize_data
try
    Eros = load('InData/ErosFit.mat');
catch
    b = 1
    cd('InData');
    MakeSpline;
    cd('..');
    close all;
    Eros = load('InData/ErosFit.mat');
end
end
