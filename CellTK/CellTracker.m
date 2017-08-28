%CellTracker.m takes a set of images and finds the trajectories of the
%paths of the cells.

function [distances] = CellTracker(fname, N)
    %tif_file contains layers of the cell movement images
    %N is sensitivity

    %Creates a struct with info like fnames, types, etc.
        %tif_file contains layers of the cell movement images
    %N is sensitivity

    %Creates a struct with info like fnames, types, etc.
    fileinfo = imfinfo(fname);
    timepoint = length(fileinfo);
    biggest=0; %Set up for matching labels later.

    %cell picture is every other layer in tif file.
    for i = 1:2:timepoint
        centers=CellCountTracker(fname,i,N);
        if length(centers)>biggest
            biggest=length(centers);
        end
        %To organize centers by image. Label, x, y, and checked.
        imcenters((i+1)/2).x=[];
        imcenters((i+1)/2).y=[];        
        imcenters((i+1)/2).label=1:length(centers);
        imcenters((i+1)/2).checked=zeros(1,length(centers));
        for j=1:length(centers)
            imcenters((i+1)/2).x=[imcenters((i+1)/2).x centers(j).Centroid(1)];
            imcenters((i+1)/2).y=[imcenters((i+1)/2).y centers(j).Centroid(2)];
        end
    end
    
    %Limit for linking cells from one layer to the next
    limit=10*N;

    for i=2:length(imcenters)
        match=-1;
        %The two frames to compare
        oldFrame=[imcenters(i-1).x; imcenters(i-1).y];
        newFrame=[imcenters(i).x; imcenters(i).y];
        for j=1:length(oldFrame) %j iterates through every cell center in old frame
            oldXY=oldFrame(:,j);
            smallest = 1e9; %find new point closest to old point
            for k=1:length(newFrame) %k iterates through not matched cell center in new frame
                if imcenters(i).checked(k)==0
                    newXY=newFrame(:,k);
                    dist=sqrt((oldXY(1)-newXY(1))^2+(oldXY(2)-newXY(2))^2);
                    if dist<smallest
                        smallest=dist;
                        match=k; %The matched label
                    end
                end
            end
            if smallest<limit %There's a match!
                imcenters(i).label(match)=imcenters(i-1).label(j);
                imcenters(i).checked(match)=1;
            else %Assign new label to unmatched cell.
                imcenters(i).label(match)=(biggest+1);
                biggest=biggest+1;
            end
        end
    end
    
    cellmvmt=cell(biggest,4); %x(list), y(list), displacement, distance

    for i=1:biggest %For each label
        for j=1:length(imcenters) %for each frame
            for k=1:length(imcenters(j).x) %for each cell in frame
                if imcenters(j).label(k)==i
                    cellmvmt{i,1}=[cellmvmt{i,1} imcenters(j).x(k)];
                    cellmvmt{i,2}=[cellmvmt{i,2} imcenters(j).y(k)]; 
                    break;
                end
            end
        end
        dist=0;
        lastPos=length(cellmvmt{i,1});
        if lastPos~=0
            %Displacement
            cellmvmt{i,3}=sqrt((cellmvmt{i,1}(1)-cellmvmt{i,1}(lastPos))^2+(cellmvmt{i,2}(1)-cellmvmt{i,2}(lastPos))^2);
            %Distance
            for m=2:length(cellmvmt{i,1})
                dist=dist+sqrt((cellmvmt{i,1}(m)-cellmvmt{i,1}(m-1))^2+(cellmvmt{i,2}(m)-cellmvmt{i,2}(m-1))^2);
            end
            cellmvmt{i,4}=dist;
        else
            cellmvmt{i,3}=0;
            cellmvmt{i,4}=0;
        end
    end
    
    cellMatrix=reshape(cellmvmt,[length(cellmvmt) 4]);
    filter=[cellMatrix{:,4}]==0 | isempty([cellMatrix{:,1}]);
    cellMatrix(filter,:)=[];
    
    %Should we make a filter for too few frames, too?
    cc=hot(length(cellMatrix));  
    dispSort = sortrows(cellMatrix, 3);
    
    distSort = sortrows(cellMatrix,4);
    
    image=imread(fname,1);
    figure;
    imshow(image); hold on;
    
    for i=1:length(dispSort)
        hold all;
        for j=length(dispSort{i,1}):-1:1
            plot(dispSort{i,1}(j), dispSort{i,2}(j),'o','color', cc(i,:));
        end
    end
    
    colormap(cc);
    c=colorbar;
    c.Ticks = [min([dispSort{:,3}]) median([dispSort{:,3}]) max([dispSort{:,3}])];
end