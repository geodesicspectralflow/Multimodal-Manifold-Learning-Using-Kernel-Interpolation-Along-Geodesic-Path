function [] = SaveFig(f,folder,filename,pos)
if nargin==3
    set(gcf, 'Position', get(0, 'Screensize'));
    set(gcf, 'Position',[1,1,1920,1007]);
else
    if not(pos==-1)
        set(gcf,'Position',pos);
    end
end
saveas(f,fullfile(folder,filename+".jpg"));

end

