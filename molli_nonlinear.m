% Calculate T1 maps from MOLLI data
% 
% Daniel Bulte, IBME, University of Oxford, July 2019
%%
% this version uses the absolute value form of the model to fit the data,
% and uses a value of 300 in the 11th volume as a noise threshold (any
% voxel <300 in volume 11 is set to zero in the T1 map)

% Edited by E Bluemke July 2020 (questions? emma.bluemke@new.ox.ac.uk)

clear all
close all

%%%%%%%%%%%%%%%%%%%
disp('Select MOLLI Folder set')
dirName = uigetdir(); 
options = struct('recursive', true, 'verbose', true, 'loadCache', false);
[partitions, meta] = readDicomSeries(dirName, options);
 % Return values:
%   imagePartitions: Array of structs containing all partitions found
%   metaFilenames: Cell array of dicom filenames that contain no images

% Read image by partition index
% readDicomSeriesImage reads a dicom image (and optional dicominfo) from a
% dicom series partition that was found with the readDicomSeries function.
%
% This function can be used in two ways:
% [image, info] = readDicomSeriesImage(directory, partition):
% Reads a partition specified by the caller. partition should be one
% element of the imagePartitions structure returned by readDicomSeries.
% directory should be the same directory that was provided to
% readDicomSeries.
%
% The image return value will contain only the frames specified in the
% partition, typically in a 3D matrix. The type is the same as returned
% from dicomread (usually int16).
%
% The info return value is either a dicominfo structure in case of an
% enhanced dicom file, or a cell array containing dicominfo structures in
% case of a series of classic dicom files.
[image1, info1] = readDicomSeriesImage(dirName, partitions(1));

nbrow = size(image1,1);
nbcol = size(image1,2);
nbslice = size(image1, 3);
nbvols = length(partitions);
nbvoxels = nbrow*nbcol*nbslice;
nbseries = length(partitions);

% disp('Select MOLLI Folder set')
% mollidir = uigetdir(); 
% molli_orig=dicomreadVolume(mollidir);

%[stat,struc]=fileattrib('*'); % gives me the names of all of the files
%nbvols=size(struc); % gives me the number of files, and thus TI's

nbti = 11; % 11 TI's for MOLLI

%metadata = zeros(nbti,1);
tinv_acq = zeros(nbti,1);
ttrig_acq = zeros(nbti,1);


for i=1:nbti
    [image, info] = readDicomSeriesImage(dirName, partitions(i));
    metadata(i) = info{1,1}; % to get the dicom headers for every file (TI)
    
    if isfield(metadata,'InversionTime')
        tinv_acq(i)=metadata(i).InversionTime;  % builds a vector of all of the TI's
    else
        ttrig_acq(i) = metadata(i).TriggerTime;  % if no inversion time then use trigger time (will be close but wrong)
        tinv_acq(i) = ttrig_acq(i);  
    end
end

% % need to reorder the TI's in tinv(i), 11 TI's in MOLLI
[tinv,new_order] = sort(tinv_acq);

% tinv(1) = tinv_acq(1);
% tinv(2) = tinv_acq(4);
% tinv(3) = tinv_acq(7);
% tinv(4) = tinv_acq(2);
% tinv(5) = tinv_acq(5);
% tinv(6) = tinv_acq(8);
% tinv(7) = tinv_acq(3);
% tinv(8) = tinv_acq(6);
% tinv(9) = tinv_acq(9);
% tinv(10) = tinv_acq(10);
% tinv(11) = tinv_acq(11);

%cd(currpath)


for k = 1:nbseries
    [image, info] = readDicomSeriesImage(dirName, partitions(k));
	dataTmp = image;
	dataTmp = double(squeeze(dataTmp));	
	for ss = 1:nbslice 
		data(:,:,ss,k) = dataTmp(:,:,ss); 
    end
end 
size(data)

% for k = 1:nbseries
% 	dataTmp = d(k).imData;
%         data(:,:,:,1) = dataTmp(:,:,:,1);
%         data(:,:,:,2) = dataTmp(:,:,:,4);
%         data(:,:,:,3) = dataTmp(:,:,:,7);
%         data(:,:,:,4) = dataTmp(:,:,:,2);
%         data(:,:,:,5) = dataTmp(:,:,:,5);
%         data(:,:,:,6) = dataTmp(:,:,:,8);
%         data(:,:,:,7) = dataTmp(:,:,:,3);  
%         data(:,:,:,8) = dataTmp(:,:,:,6);
%         data(:,:,:,9) = dataTmp(:,:,:,9);
%         data(:,:,:,10) = dataTmp(:,:,:,10);
%         data(:,:,:,11) = dataTmp(:,:,:,11);
% end 
%         

% create a mask to speed up calc, thresholds data 300 in final volume
mask=data(:,:,:,11);
mask(le(mask,300))=0; % changed from 100 to 300 20190801 DB
mask(ge(mask,300))=1;


% initialise matrices 
t1vec = zeros(1,nbvoxels,'single');
slope = zeros(nbti,nbvoxels); % there are 11 TI's in the MOLLI sequence

%% Calculate T1

indechs = 1;

% create a 2D array with TI's as the 2nd dimension
for z=1:nbslice
    for y=1:nbcol
        for x=1:nbrow  
        if (mask(x,y,z)==1)
            slope(:,indechs) = data(x,y,z,:);
        end
        indechs = indechs + 1;
        end
    end 
end


fo = fitoptions('Method','NonlinearLeastSquares','Lower',[0,0,0],'Upper',[6000,12000,5000],'StartPoint',[1000,2000,1000]);

molli = fittype('abs(Axy - Bxy * exp(-tinv/tonestar))','dependent',{'y'},...
    'independent',{'tinv'},'coefficients',{'Axy','Bxy','tonestar'},'options',fo);


for i=1:nbvoxels
    recover = slope(:,i);
    i
    if recover(1)~=0
                f = fit(tinv,recover,molli); 
                coeffvals = coeffvalues(f);
                Tonestar =  coeffvals(3);
                t1vec(i)= Tonestar*(coeffvals(2)/coeffvals(1)-1); % LL correction

            if (isnan(t1vec(i)) || t1vec(i)<0 || isinf(t1vec(i) || t1vec(i)>8000)) % remove rubbish values, limit to 5sec max
                    t1vec(i)=0;
            end
    end
end

t1map = reshape(t1vec,nbrow,nbcol,nbslice);

figure
imagesc(t1map);
title(sprintf('%s-%d.png',dirName,i))
caxis([0,3500]);
colorbar
axis off
daspect([1 1 1])
saveas(gcf,sprintf('%s.png',dirName))


dicomt1map = uint16(reshape(t1map,[nbrow nbcol 1 nbslice])); % reshape undoes the squeeze, which removed the colour dimension

output_dcm = [dirName, '_MOLLI_T1_map.dcm'];
% output_nii = [loadpath '_MOLLI_T1_map.nii.gz'];
%cd(dirname)
dicomwrite(dicomt1map, output_dcm, metadata(1), 'CreateMode', 'Copy'); % save as a dicom, gets metadata from another dicom file
% niftiwrite(dicomt1map, output_nii, 'Compressed', true);

beep
%% end

