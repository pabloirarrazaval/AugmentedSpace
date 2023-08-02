function mkimshowsubplot(type,aux,rangetype,name,varargin)
% mkimshowsubplot(type,aux,rangetype,name)
%
% Displays three images in same plot
% type: 'complex' - Shows abs, real and imaginary of same image
%       'compare' - Shows three different images (real only) NxMx3
% aux: the data
% rangetype: 'common' - Shows the three images in the same scale
%            'individual' - Each image has its own scale
%            [m M] - Use this scaling for the three images
% name: to put in the titles (celss for comparison)

if strcmp(type,'complex') % Shows abs, real and imaginary of same image

    if ischar(rangetype)
        if strcmp(rangetype,'common')
            MM = max([abs(aux(:)); real(aux(:)); imag(aux(:))]);
            mm = min([abs(aux(:)); real(aux(:)); imag(aux(:))]);
            range = [mm MM];
        elseif strcmp(rangetype,'individual')
            range = [];
        else
            error('range must be "common", "individual" or range in [].')
        end
    else
        range = rangetype;
    end

    subplot('Position',[2.5 2.5 30 95]/100);
    imshow(abs(aux),range,varargin{:});
    set(gca,'YDir','Normal');
    title([name ' (abs)']);

    subplot('Position',[35 2.5 30 95]/100)
    imshow(real(aux),range,varargin{:});
    set(gca,'YDir','Normal');
    title([name ' (real)']);

    subplot('Position',[67.5 2.5 30 95]/100)
    imshow(imag(aux),range,varargin{:});
    set(gca,'YDir','Normal');
    title([name ' (imag)']);

elseif strcmp(type,'compare') % Shows three images (real only) (3 channels)

    if ischar(rangetype)
        if strcmp(rangetype,'common')
            MM = max(aux(:));
            mm = min(aux(:));
            range = [mm MM];
        elseif strcmp(rangetype,'individual')
            range = [];
        end
    else
        range = rangetype;
    end

    subplot('Position',[2.5 2.5 30 95]/100);
    imshow(aux(:,:,1),range,varargin{:});
    set(gca,'YDir','Normal');
    title(name{1});

    subplot('Position',[35 2.5 30 95]/100)
    imshow(aux(:,:,2),range,varargin{:});
    set(gca,'YDir','Normal');
    title(name{2});

    subplot('Position',[67.5 2.5 30 95]/100)
    imshow(aux(:,:,3),range,varargin{:});
    set(gca,'YDir','Normal');
    title(name{3});
end

end