
function simpleChanTopo(sensors, chanlabels, data, cmap, comment, fontcolor, axiscol)
% simpleChanTopo(sensors, chanlabels, data, cmap, comment, fontcolor, axiscol)
        cfg = [];

        cfg.marker = 'labels'; % 'on', 'off', 'labels', 'numbers';
        cfg.markersymbol = '.';    
        cfg.markersize = 4;
        cfg.labelcolor = fontcolor;
        cfg.elec = sensors;
        cfg.elec.pnt = cfg.elec.pnt(:,[1 2 3]);

        cfg.layout = ft_prepare_layout(cfg);       %Layout of the electrodes

        cfg.xparam = 'time';
        cfg.zparam = 'avg';
        if ~isempty(cmap)
        cfg.colormap = cmap;
        end
        cfg.colorbar = 'yes';
        cfg.hotkeys = 1;

        data_eeg.time = 1;
        data_eeg.dimord = 'chan_time';
        data_eeg.label = chanlabels;
        data_eeg.avg = data;
        
        if ~isempty(axiscol)
            cfg.zlim = axiscol;
        end
        cfg.interactive = 'no';
        cfg.interpolation = 'v4'; % 'linear','cubic','nearest','v4' (default = 'v4')
        cfg.style =  'straight';
        cfg.interplimits = 'head';
        cfg.gridscale = 80; % resolution of figure (default = 67)

        cfg.commentpos = 'title';
        cfg.comment = comment;
        
        ft_topoplotER(cfg, data_eeg);
        
        title(cfg.comment)
end