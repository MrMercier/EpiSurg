printEm=0;

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='DK';
cfg.title='PT001: DK Atlas'; 
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasDKnoLabels'); 
end

%%
parcOut=elec2Parc('PT001','DK');

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='D';
cfg.title='PT001: Destrieux Atlas';
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasDleftLat'); 
end

%%
parcOut=elec2Parc('PT001','D');

%%
createIndivYeoMapping('PT001');

%%
cfg=[];
cfg.view='l';
cfg.overlayParcellation='Y17';
cfg.title='PT001: Yeo 17-Area';
cfgOut=plotPialSurf('PT001',cfg);
if printEm
   print(gcf,'-djpeg','pt001atlasY17'); 
end


%%
parcOut=elec2Parc('PT001','Y7');

%%
disp('Script testAtlases.m completed successfully.')