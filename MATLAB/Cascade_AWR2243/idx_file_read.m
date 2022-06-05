idxFile = fopen(adcIdxFileName,'r');

indexInfo = fread(idxFile, 'uint64');

timeStampFrame = indexInfo(8:6:end); % data size for the effective number of frames

fclose(idxFile);