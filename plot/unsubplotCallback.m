function unsubplotCallback(varargin)
    
    hh = varargin{3};  % Get structure.
    
    % create new image
    f = figure;
    
    % copy curret subplot
    c = copyobj(hh,get(hh,'Parent'));
    set(c,'Parent',f);
    set(c,'Position','default');
    
    % resize stuff
    ha=gca;
    hay=get(ha,'Ylabel');
    hax=get(ha,'Xlabel');
    NewFontSize=10;
    set(hay,'Fontsize',NewFontSize);
    set(hax,'Fontsize',NewFontSize);
end
