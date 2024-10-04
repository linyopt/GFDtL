figure(1)
set(gcf, 'Color', [1,1,1]);
subplot(3,1,1);
plotCP(cpG,cpI,y,0);
title('Simulated data and changepoints detected with IFGL (red) and GFGL (blue)');
subplot(3,1,2);
coefCompare(edges(1:2,:),ThetaG,sigmainv);
title('Example of edge estimation for GFGL (estimates dashed ground-truth solid)');
xlabel('time');
ylabel('\Theta_{i,j}');
subplot(3,1,3);
coefCompare(edges(1:2,:),ThetaI,sigmainv);
title('Edge estimation for IFGL (same as above)');
xlabel('time');
ylabel('\Theta_{i,j}');

figure(2)
set(gcf, 'Color', [1,1,1]);
subplot(3,3,1);
plotGraph(squeeze(abs(sigmainv(:,:,3))),10,0);
title('Ground-truth graph (t=3)')
subplot(3,3,2);
plotGraph(squeeze(abs(sigmainv(:,:,7))),10,0);
title('(t=7)')
subplot(3,3,3);
plotGraph(squeeze(abs(sigmainv(:,:,12))),10,0);
title('(t=12)');

subplot(3,3,4);
plotGraph(squeeze(abs(ZG(:,:,3))),10,0);
title('GFGL estimate');
subplot(3,3,5);
plotGraph(squeeze(abs(ZG(:,:,7))),10,0);
subplot(3,3,6);
plotGraph(squeeze(abs(ZG(:,:,12))),10,0);

subplot(3,3,7);
plotGraph(squeeze(abs(ZI(:,:,3))),10,0);
title('IFGL estimate');
subplot(3,3,8);
plotGraph(squeeze(abs(ZI(:,:,7))),10,0);
subplot(3,3,9);
plotGraph(squeeze(abs(ZI(:,:,12))),10,1);