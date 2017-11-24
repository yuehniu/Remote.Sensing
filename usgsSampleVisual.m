%% read binary usgs file
fileHandle= fopen('./data/cup95RefEm.sli', 'rb');
cupLib = fread(fileHandle, [50, 38], 'single');
cupSpectrum  = [1.990800, 2.000900, 2.010900, 2.020900, 2.030900, 2.040900, 2.050900, ...
 2.060900, 2.071000, 2.081000, 2.091000, 2.101000, 2.111000, 2.121000, ...
 2.130900, 2.140900, 2.150900, 2.160900, 2.170900, 2.180900, 2.190800, ...
 2.200800, 2.210800, 2.220800, 2.230700, 2.240700, 2.250600, 2.260600, ...
 2.270600, 2.280500, 2.290400, 2.300400, 2.310400, 2.320300, 2.330200, ...
 2.340200, 2.350100, 2.360000, 2.370000, 2.379900, 2.389800, 2.399700, ...
 2.409600, 2.419600, 2.429500, 2.439400, 2.449300, 2.459200, 2.469100, ...
 2.479000];

figure;
hold on
xlim([cupSpectrum(1), cupSpectrum(50)])
xlabel('Wavelength(\mum)')
ylabel('Reflectance')

emRef = cupLib(:,5);
strid_ = 2;
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'b-^', 'LineWidth', 2);

emRef = cupLib(:,11);
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'k-+', 'LineWidth', 2);

emRef = cupLib(:,2);
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'r-d', 'LineWidth', 2);

emRef = cupLib(:,1);
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'g-*', 'LineWidth', 2);

emRef = cupLib(:,4);
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'm-s', 'LineWidth', 2);

emRef = cupLib(:,6);
plot(cupSpectrum(1:strid_:50), emRef(1:strid_:50), 'c-o', 'LineWidth', 2);

legend('Nontronite', 'Chalcedony','Kaolinite','Alunite','Andradite','Muscovite')