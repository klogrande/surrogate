%Simplified Parsing of Combustion Reaction Data Sets
%Author: Kevin LoGrande
%Description: Parses data to workable size by setting a minimum step size
              %in norm to maintain somewhat uniform data spacing. Then,
              %the data is normalized by the max column number to be
              %prepared for use in the diffusion maps algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%path to folder where files exist%
filePath='C:\Users\klogr\Documents\ARL\ARLTarballs\ExtractedFiles\2.0e_2.2dens\';

%file numbering, time step from the DPD simulation, and step minimum
firstFileNum=0;
lastFileNum=500; %Leave some files for testing data
totFiles=lastFileNum-firstFileNum+1;
DPDtimeStep=0.005;
minStep=0.005;

%initialize data matrices
parsedData=zeros(1,9);
nextData=zeros(1,10);

%determine data to read in
for fileNum=firstFileNum:lastFileNum
    
    %create string with name of file to read
    fileName=sprintf('2.0e_2.2dens_run_%i.csv', fileNum);
    fileVar=strcat(filePath,fileName);
    
    %read in full data from individual file as a matrix
    dataFile=readmatrix(fileVar,'Range','C2:K20002');
    
    %step through line by line and decide which data to keep
    row=1;
    numSteps=1;
    while row+1 <= size(dataFile,1)
        
        %if not a sizable step, delete the data
        if norm(dataFile(row,:)-dataFile(row+1,:))<minStep
            dataFile(row+1,:)=[];
            numSteps=numSteps+1;
        
        %if step is accepted, store appropriate data    
        else
            nextData2=[nextData; dataFile(row+1,:) numSteps*DPDtimeStep];
            clear nextData
            nextData=nextData2;
            clear nextData2
            row=row+1;
            numSteps=1;
        end 
    end
    
    %delete last row of dataFile since we have no derivative info for it
    dataFile(size(dataFile,1),:)=[];
    
    %store parsed data file in larger data matrix
    parsedData2=[parsedData;dataFile];
    clear parsedData
    parsedData=parsedData2;
    clear parsedData2
    clear dataFile
end
%}

%eliminate extraneous first rows
parsedData(1,:)=[];
nextData(1,:)=[];

%normalize data by its max value
[m,n]=size(parsedData);
maxVal=max(parsedData);
multiplier=1./maxVal;
multiplier=diag(multiplier);
normData=parsedData*multiplier;

%calculate discrete derivatives
deriv=zeros(m,n);
for row=1:m
    deriv(row,:)=(nextData(row,1:9)-parsedData(row,:))/(nextData(row,10));
end
%normalize derivative values
normDataDer=deriv*multiplier;