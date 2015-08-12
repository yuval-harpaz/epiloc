function [ind,Vmodel,R]=dipoleFitAnat(Vdata,gain,pos,gainX)
% input
% Vdata is the measured field, N channels by 1. 
% gain matrix - N channels by M sources.
% gainX is the gain of mirror dipoles if you want to consider some activity
% on the other side

% output
% ind is the index of the source (from M sources)
% Vmodel is the reconstructed field from source estimation
% R is square correlation between measured and estimated fields
% 
%% consider xhemi dipole (WORKED WELL)




increment=1;
radius=sqrt(100./pi); %radius for area = 10cm^2
disp(['scanning through ',num2str(size(gain,2)),'sources  :  '])
for pnti=1:size(gain,2)
    % FIXME why are they (gainX and R) in different scales?
    mom=pinv([gainX,gain])*Vdata;;
    reconLR=[gainX,gain]*mom;
    R(pnti)=corr(reconLR,Vdata).^2;
    if R(pnti)==max(R)
        ind=pnti;
        maxR=R(pnti);
        Vmodel=reconLR;
    end
    prog(pnti)
end



