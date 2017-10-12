%% read binary usgs file
fileHandle= fopen('../data/Cup95RefEm.sli', 'rb');
cupLib = fread(fileHandle, [50, 38], 'single');
cupSpectrum  = 1000 * [1.990800, 2.000900, 2.010900, 2.020900, 2.030900, 2.040900, 2.050900, ...
 2.060900, 2.071000, 2.081000, 2.091000, 2.101000, 2.111000, 2.121000, ...
 2.130900, 2.140900, 2.150900, 2.160900, 2.170900, 2.180900, 2.190800, ...
 2.200800, 2.210800, 2.220800, 2.230700, 2.240700, 2.250600, 2.260600, ...
 2.270600, 2.280500, 2.290400, 2.300400, 2.310400, 2.320300, 2.330200, ...
 2.340200, 2.350100, 2.360000, 2.370000, 2.379900, 2.389800, 2.399700, ...
 2.409600, 2.419600, 2.429500, 2.439400, 2.449300, 2.459200, 2.469100, ...
 2.479000];
% extract neccesary band
figure(1)

subplot(231)
emMdc = H_test_orig(:,7)/1000;
emNfindr = H_I_orig(:,7)/1000;
emRef = cupLib(:,5);
plot(cupSpectrum,emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr , 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Nontronite', 'FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)

subplot(232)
emNfindr = H_I_orig(:,3)/1200;
emRef = cupLib(:,11);
emMdc = H_test_orig(:,3)/1400;
plot(cupSpectrum, emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr , 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Chalcedony','FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)

subplot(233)
emMdc = H_test_orig(:,4)/500;
emNfindr = H_I_orig(:,4)/600;
emRef = cupLib(:,2);
plot(cupSpectrum, emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr , 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Kaolinite', 'FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)

subplot(234)
emMdc = H_test_orig(:,5)/600;
emNfindr = H_I_orig(:,5)/600;
emRef = cupLib(:,1);
plot(cupSpectrum, emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr, 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Alunite','FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)

subplot(235)
emMdc = H_test_orig(:,8)/400;
emNfindr = H_I_orig(:,8)/400;
emRef = cupLib(:,4);
plot(cupSpectrum, emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr, 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Andradite', 'FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)

subplot(236)
emMdc = H_test_orig(:,9)/400;
emNfindr = H_I_orig(:,9)/400;
emRef = cupLib(:,6);
plot(cupSpectrum, emRef, 'k--', 'LineWidth', 2)
ylim([0, 1])
hold on
plot(cupSpectrum, emMdc, 'r', 'LineWidth', 2)
plot(cupSpectrum, emNfindr, 'g', 'LineWidth', 2)
ylabel('reflectance(%)')
xlabel('wavelength (nm)')
title('Muscovite', 'FontSize', 20)
legend('Lib', 'MDC', 'NFINDR')
sad(emRef, emMdc)
sad(emRef, emNfindr)
% plot(H_test_orig(:,3)/200, '--', 'LineWidth', 2)
fclose(fileHandle);