function [ramsize] = mysizeRAMmemory(nelements)

if strcmp(classtype,'single')
    bytes = 4;
elseif strcmp(classtype,'double')
    bytes = 8;
end
bytes2megas = 1d-6; %[megas * bytes-1]
bytes2gigas = 1d-9; %[gigas * bytes-1]

ramsizeMegas = nelements*bytes) * bytes2megas;
ramsizeGigas = nelements*bytes) * bytes2gigas;

ramsize = ramsizeMegas;
