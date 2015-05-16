if ishandle(1), set(0, 'CurrentFigure', 1); else figure(1); end
clf
set(gcf,'DefaultLineLineWidth',2,'DefaultTextFontSize',12,...
    'DefaultTextFontWeight','bold','DefaultAxesFontSize',12,...
    'DefaultAxesFontWeight','bold');

plottype=1;

if plottype==1
    % Plot eta and density(z-eta)
    subplot(1,2,1); imagesc(x,z,eta); title('eta (m)');
    subplot(1,2,2); imagesc(x,z,density); title('density');
    for ii=1:2
        subplot(1,2,ii); set(gca,'ydir','normal');
        axis([0 L -H 0]); colorbar; xlabel('x (m)'); ylabel('z (m)');
    end
elseif plottype==2
    % Plot eight fields
    subplot(4,2,1); imagesc(x,z,eta);       title('eta (m)');
    subplot(4,2,2); imagesc(x,z,density);   title('density');
    subplot(4,2,3); imagesc(x,z,uwave);     title('u (wave) (m/s)');
    caxis([-1 1]*max(abs(uwave(:))));
    subplot(4,2,4); imagesc(x,z,w);         title('w (m/s)');
    caxis([-1 1]*max(abs(w(:))));
    subplot(4,2,5); imagesc(x,z,kewave);    title('kewave (m^2/s^2)');
    caxis([0 1]*max(abs(kewave(:))));
    subplot(4,2,6); imagesc(x,z,apedens);   title('ape (m^2/s^2)');
    caxis([0 1]*max(abs(apedens(:))));
    subplot(4,2,7); imagesc(x,z,ri);        title('Ri');
    caxis([0.2 1]);
    subplot(4,2,8); imagesc(x,z,vorticity); title('vorticity (1/s)');
    caxis([-1 1]*max(abs(vorticity(:))));
    for ii=1:8
        subplot(4,2,ii); set(gca,'ydir','normal');
        axis([0 L -H 0]); colorbar; xlabel('x (m)'); ylabel('z (m)');
    end
end
drawnow
