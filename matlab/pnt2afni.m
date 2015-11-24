function pnt2afni(pnt,mri)
% creates an AFNI file with the digitization points
% 'small' voxSize creates a hs+orig file with small digitization
% points
% 'big' voxSize (default) creates HS+orig file with 5x5x5mm cubes around the points
%requires hs_file and ortho+orig or warped+orig
% mri, file name instead of ortho/warped

%voxSize='small';

PNT=reshape(pnt',size(pnt,1)*3,1);
if exist('pntTxt','file')
    !rm pntTxt
end
txtFileName = 'pntTxt';
fid = fopen(txtFileName, 'w');
fprintf(fid,'%f\t%f\t%f\n',PNT);
fclose(fid);
if exist('pnt+orig.BRIK','file')
    !rm pnt+orig*
end

if exist('mri','var')
    eval(['!~/abin/3dUndump -orient PRI -xyz -dval 1 -master ',mri,' -prefix pnt pntTxt']);
else
    if exist ('ortho+orig.BRIK','file');
        !~/abin/3dUndump -orient PRI -xyz -dval 1 -master ortho+orig -prefix pnt pntTxt
    elseif exist ('warped+orig.BRIK','file');
        !~/abin/3dUndump -orient PRI -xyz -dval 1 -master warped+orig -prefix pnt pntTxt
    else
        display('cannot find ortho or warped+orig file')
        return
    end
end
% if strcmp(voxSize,'small')
%     return
% end
% !~/abin/3dresample -dxyz 5 5 5 -prefix hsT -inset hds+orig -rmode Cu
% !~/abin/3dfractionize -template hsT+orig -input hds+orig -prefix HS
% !rm hds+orig*
% !rm hsT+orig*

%!~/abin/afni -dset ortho+orig
end
