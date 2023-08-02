function imshowzoomed(x,zoomarea)

% x = abs(m_dft_csm);
% zoomarea = [42 59 30 30]; % [xpos ypos width height]
zoomfactor = 2.5;
flag = false;
if (zoomarea(1)+zoomarea(3)-1)>size(x,1)
    zoomarea(3) = size(x,1)-zoomarea(1); flag = true;
end
if (zoomarea(2)+zoomarea(4)-1)>size(x,2)
    zoomarea(4) = size(x,2)-zoomarea(2); flag = true;
end
if flag
    fprintf('WARNING: zoom area reduced to %d x %d.\n',zoomarea(3),zoomarea(4));
end
imshow(x,[],'Init',300);set(gca,'YDir','Normal');
ax = gca;
inset = x(zoomarea(1)+(0:zoomarea(3)-1),zoomarea(2)+(0:zoomarea(4)-1));
inset = inset/max(inset(:))*ax.CLim(2);
hold on
xdata = [size(x,1)-round(zoomfactor*zoomarea(3)+1) size(x,1)];
ydata = [size(x,2)-round(zoomfactor*zoomarea(4)+1) size(x,2)];
imshow(inset,ax.CLim,'Init',300,'XData',xdata,'YData',ydata,'Interpolation','bilinear');
set(gca,'YDir','Normal');
hold off

end

