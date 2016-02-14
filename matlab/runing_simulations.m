noiseFactor=0.1; 
for dipi=1:5
    disp(num2str(dipi));
    marikVirtual29(dipi, noiseFactor, 100);
    close all
end


for noiseFactor=[0.1 0.3] 
    for dipi=1:5
        disp(num2str(dipi));
        marikVirtual29(dipi, noiseFactor, 1000);
        close all
    end
end


for noiseFactor=[0.1 0.3] 
    for dipi=1:5
        disp(num2str(dipi));
        marikVirtual29(dipi, noiseFactor, 1);
        close all
    end
end

% sensitivity 
% specificity 
% how much false postive??
% how much false negative??

% Intensity

% adding additional dipole in less confidence