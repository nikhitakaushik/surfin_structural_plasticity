function [template, t, shifted, h] = make_stable_slice_average(im,varargin)
    I1=1:floor(size(im,3)/2);
    I2=floor((size(im,3)/2))+1:size(im,3);
    
    template=mean(im(:,:,I2),3);
    shifted=zeros(size(im));
    
    limits = [size(im,2)/8,size(im,1)/8]; 
    
    for i=I1
        [t(i,:), h(:,:,i)]=cvMotionCorrect(im(:,:,i),template,varargin{:});
        t(i,t(i,:)==limits) = 0;
        shifted(:,:,i)=BilinearImageRegistrator.shift(im(:,:,i),t(i,:));
    end
    
    template=nanmean(shifted(:,:,I1),3);
    for i=I2
        [t(i,:), h(:,:,i)]=cvMotionCorrect(im(:,:,i),template,varargin{:});
        t(i,t(i,:)==limits) = 0;
        shifted(:,:,i)=BilinearImageRegistrator.shift(im(:,:,i),t(i,:));
    end
    template=nanmean(shifted(:,:,I2),3);
    for i=I1
        [t(i,:), h(:,:,i)]=cvMotionCorrect(im(:,:,i),template,varargin{:});
        t(i,t(i,:)==limits) = 0;
        shifted(:,:,i)=BilinearImageRegistrator.shift(im(:,:,i),t(i,:));
    end
    template=nanmean(shifted,3);







