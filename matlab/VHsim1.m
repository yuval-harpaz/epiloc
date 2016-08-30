cd /home/yuval/Data/marik/som2/1

hdr=ft_read_header('c,rfhp0.1Hz');
hs=ft_read_headshape('hs_file');

xyz=hdr.grad.coilpos

plot3pnt(hdr.grad.coilpos,'.k')
hold on
plot3pnt(hdr.grad.coilori,'.r')


theta=pi/180*10;
Rx=[1,0,0;0,cos(theta),-sin(theta);0,-sin(theta),cos(theta)];
coilpos=hdr.grad.coilpos(1:248,:);
coilposNew=[Rx*coilpos']';
figure;
plot3pnt(coilpos,'.k')
hold on
plot3pnt(coilposNew,'.g')

coilori=hdr.grad.coilori(1:248,:)./100;
coiloriNew=[Rx*coilori']';
figure;
plot3pnt(coilposNew,'.k');
hold on
plot3pnt(coilposNew+coiloriNew,'.m');