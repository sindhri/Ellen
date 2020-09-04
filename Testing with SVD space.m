% I am modifying Woods et al code to make some assessments of our PCA
% space.

%First, the PCA data is (this is slightly different than Schoppik's use of single value decomposition, because I have not mean centerered the data to be zero:

Fingerprint.coefs
Fingerprint.scores
Fingerprint.variances

Fingerprint.sr %is the Fingerprint.impute data divided by the standard dev.
Fingerprint.impute %is the imputed data of Fingerprint.data

% First I want to try David's code once more by mean centering the result:

for i = 1:120, 
  logdata(i,:) = Fingerprint.impute(i,:)-mean(Fingerprint.impute(i,:)); % make sure there is no DC offset)
  normdata(i,:) = logdata(i,:)./std(logdata(i,:)); 
end

% a figure describing the raw data
figure
imagesc(normdata,[-3 3])
set(gca,'position',[0 0 1 1],'visible','off')
set(gcf,'position',[100 100 800 350],'paperposition',[1 1 8 3.5])
colormap(map)
print(gcf,'-depsc2','woods1.eps')

% better to annotavte/reshape this in illustrator

% decompose the data matrix (i.e. do principal components analysis
[u fulls fullv] = svd(normdata);

% see how much of the variance we explain (i.e. a scree plot)
tmp = 100*(diag(fulls)./sum(diag(fulls)));

% and get the projections-- but I will try the first 10
fullproj = normdata*fullv(:,1:10);

% and the vector lengths
fulllen = sqrt(fullproj(:,1).^2 + fullproj(:,2).^2 + fullproj(:,3).^2+ fullproj(:,4).^2 + fullproj(:,5).^2 + fullproj(:,6).^2+ fullproj(:,7).^2 + fullproj(:,8).^2 + fullproj(:,9).^2+ fullproj(:,10).^2);


figure
hold on
plot(1,tmp(1),'r.','markersize',18)
plot(2,tmp(2),'g.','markersize',18)
plot(3,tmp(3),'b.','markersize',18)
plot(4:length(tmp),tmp(4:end),'k.','markersize',18)
plot([0 10],[10 10],'k:')
axis square,box off
set(gca,'xlim',[0 10],'ylim',[0 40],'ytick',[0 20 40],'xtick',[],'fontsize',12)
xlabel('Eigenvector')
ylabel('% Variance')
stevify
print(gcf,'-depsc2','woods2.eps')

% Does NOT explain that much of the data


% now, how did the first three do?
sprintf('The first three components account for %g of the variance',sum(tmp(1:3)))

% ok, now let's plot the first three eigenvectors
figure
hold on
plot(fullv(order,10),'r','linewidth',2);
plot(fullv(order,11),'g','linewidth',2);
plot(fullv(order,12),'b','linewidth',2);
plot([0 26],[0 0],'k:')
plot([8 8],[-.8 .8],'k:')
axis square, box off
set(gca,'xlim',[0 26],'ylim',[-.8 .8],'visible','off')
stevify
print(gcf,'-depsc2','woods3.eps')

% probably good to annotate the interesting peaks with text detailing what they are
% but that is easier to do in illustrator

% now, let's get to testing the generality of this space
% if we leave one peptide out, and redo the analysis
% then the spaces should be similar.
% the way we'll test that is to see 
% a) first, are the spaces the same (i.e. are the vectors similar?)
% b) is each peptide similarly far from the origin?
% c) is each peptide similarly far from its neighbor
% we'll use the Euclidian distance (i.e. sqrt of the sum of squares) as a
% distance metric.

for i = 1:120
 dex = circshift(1:120,[0 i]); % choose a unique set for each run
 [u s v] = svd(normdata(dex(1:119),:)); % take the first six
 
 % how much variance is accounted for?
 tmp = 100*(diag(s)./sum(diag(s)));
 varacct(i) = sum(tmp(1:10)); 
 
 % how similar are the vectors that define the space?
 eigvecproj(i,:) = diag(v(:,1:10)'*fullv(:,1:10));
 
 % save the indiviudal projections
 skiponev(i,:,:) = v;
 
 % get the new projections
 proj(:,:,i) = normdata*v(:,1:10);
 
 % now, how far is each point from the origin in this new space?
 veclen(i,:) = sqrt(proj(:,1,i).^2 + proj(:,2,i).^2 + proj(:,3,i).^2+proj(:,4,i).^2 + proj(:,5,i).^2 + proj(:,6,i).^2+proj(:,7,i).^2 + proj(:,8,i).^2 + proj(:,9,i).^2+proj(:,10,i).^2);
 
 % and how far is each point from each others?
 spacing(i,:) = pdist(proj(:,:,i));
end

% on the fourth run, eigenvectors 2 and 3 end up switched. so we'll
% fix that manually
if size(data,2) == 16
eigvecproj(4,2) = skiponev(4,:,2)*fullv(:,3);
eigvecproj(4,3) = skiponev(4,:,3)*fullv(:,2);
end

if size(data,2) == 17
eigvecproj(5,2) = skiponev(5,:,2)*fullv(:,3);
eigvecproj(5,3) = skiponev(5,:,3)*fullv(:,2);
end

% now, how did the first three do?
sprintf('The first three components in the leave one out analysis\naccount for %g of the variance, SD = %g',mean(varacct),std(varacct))


numruns = 10000;
randveclen = zeros(numruns,120);
randspacing = zeros(numruns,7140);
randeigvecproj = zeros(numruns,10);

for i = 1:numruns
  randnormdata = reshape(normdata(randperm(numel(normdata))),[size(data,1) size(data,2)]);
  [u s v] = svd(randnormdata);
  randeigvecproj(i,:) = diag(v(:,1:10)'*fullv(:,1:10));
  randproj = normdata*v(:,1:10);
  randveclen(i,:) = sqrt(randproj(:,1).^2 + randproj(:,2).^2 + randproj(:,3).^2+randproj(:,4).^2 + randproj(:,5).^2 + randproj(:,6).^2+randproj(:,7).^2 + randproj(:,8).^2 + randproj(:,9).^2+randproj(:,10).^2);
  % and how far is each point from each others?
  randspacing(i,:) = pdist(randproj);
end

% For random data, we'll just shuffle the original data so
% that we use the same values, but without any structure.

% first, the projections
bins = 0:.01:1;
P = hist(abs(eigvecproj(:)),bins);
randP = hist(abs(randeigvecproj(:)),bins);
P = P./sum(P);
randP = randP./sum(randP);

figure
hold on
plot(bins,cumsum(P),'k','linewidth',2);
plot(bins,cumsum(randP),'color',[.6 .6 .6],'linewidth',2)
axis square
set(gca,'xlim',[0 1],'ylim',[0 1], 'fontsize',12,'xtick',[0 .5 1],'ytick',[0 .5 1])
xlabel('Similarity of Eigenvectors')
ylabel('Cumulative Probability')
stevify 
print(gcf,'-depsc2','woods4.eps')

% that seems pretty clear; the vectors are largely the same even when we've
% left out one peptide. now, is each peptide about the same distance from
% the origin as it would have been in the original space?
bins = 0:.1:5;
randy = hist(randveclen(:),bins);
y = hist(veclen(:),bins);

figure
hold on
plot(bins,cumsum(randy./sum(randy)),'color',[.6 .6 .6],'linewidth',2)
plot(bins,cumsum(y./sum(y)),'color','k','linewidth',2)
axis square
box off
set(gca,'xlim',[0 5],'xtick',[0  5],'ylim',[0 1],'ytick',[0 .5 1],'fontsize',12)
xlabel('Length')
ylabel('Cumulative Probability')
stevify
print(gcf,'-depsc2','woods5.eps')


% looks pretty good too; now, what about the inter-point distance?
% here, we want to see how each value compares to the actual distance in
% the full space

% first, the random values
randx = repmat(pdist(fullproj),[numruns 1]);
x = repmat(pdist(fullproj),[120 1]);

% there's really no reason to plot all 10000 points
dex = randperm(numel(randx));
dex = dex(1:1000);

figure
hold on
plot(randx(dex),randspacing(dex),'.','color',[.6 .6 .6])
plot(x(:),spacing(:),'k.')
plot([0 10],[0 10],'k:')
axis square
box off
set(gca,'xlim',[0 10],'ylim',[0 8],'xtick',[0 10],'ytick',[0 8],'fontsize',12)
xlabel('True Distance')
ylabel('Distance')
stevify
print(gcf,'-depsc2','woods6.eps')

% OK, running this dataset gives pretty good graphic that are clean-- the
% key was to make sure for SVD, the space was centered with a mean zero.
% However, in the scree plot, it does not assess as much of the variance as
% the princomp analysis does.  Thus, I now try to redo this again, using my
% original PCA space as the testing set and the randomizing set:


% This time, I want to work in the PCA space I already generated, so, first
% off, the normdata should be equal to Fingerprint.sr

normdata=Fingerprint.sr;

%for i = 1:120, 
%  logdata(i,:) = Fingerprint.impute(i,:)-mean(Fingerprint.impute(i,:)); % make sure there is no DC offset)
%  normdata(i,:) = logdata(i,:)./std(logdata(i,:)); 
%end

% a figure describing the raw data
figure
imagesc(normdata,[-3 3])
set(gca,'position',[0 0 1 1],'visible','off')
set(gcf,'position',[100 100 800 350],'paperposition',[1 1 8 3.5])
colormap(map)
print(gcf,'-depsc2','woods1.eps')

% better to annotavte/reshape this in illustrator

% decompose the data matrix (i.e. do principal components analysis
% This is already done, where Fingerprint.variances is the variance,
% Fingerprint.score is the data mapped into the new PCA space, and
% Fingerprint.coefs is the coeffecients

%[u fulls fullv] = svd(normdata);

% see how much of the variance we explain (i.e. a scree plot)
tmp = 100*((Fingerprint.variances)./sum(Fingerprint.variances));

% and get the projections-- but I will try the first 10
fullproj = Fingerprint.scores(:,1:10);

% and the vector lengths
fulllen = sqrt(fullproj(:,1).^2 + fullproj(:,2).^2 + fullproj(:,3).^2+ fullproj(:,4).^2 + fullproj(:,5).^2 + fullproj(:,6).^2+ fullproj(:,7).^2 + fullproj(:,8).^2 + fullproj(:,9).^2+ fullproj(:,10).^2);


figure
hold on
plot(1,tmp(1),'r.','markersize',18)
plot(2,tmp(2),'g.','markersize',18)
plot(3,tmp(3),'b.','markersize',18)
plot(4:length(tmp),tmp(4:end),'k.','markersize',18)
plot([0 10],[10 10],'k:')
axis square,box off
set(gca,'xlim',[0 10],'ylim',[0 40],'ytick',[0 20 40],'xtick',[],'fontsize',12)
xlabel('Eigenvector')
ylabel('% Variance')
stevify
print(gcf,'-depsc2','woods2.eps')

% Does NOT explain that much of the data


% now, how did the first three do?
sprintf('The first three components account for %g of the variance',sum(tmp(1:3)))

% ok, now let's plot the first three eigenvectors
figure
hold on
plot(Fingerprint.coefs(:,1),'r','linewidth',2);
plot(Fingerprint.coefs(:,2),'g','linewidth',2);
plot(Fingerprint.coefs(:,3),'b','linewidth',2);
plot([0 26],[0 0],'k:')
plot([8 8],[-.8 .8],'k:')
axis square, box off
set(gca,'xlim',[0 26],'ylim',[-.8 .8],'visible','off')
stevify
print(gcf,'-depsc2','woods3.eps')

% probably good to annotate the interesting peaks with text detailing what they are
% but that is easier to do in illustrator

% now, let's get to testing the generality of this space
% if we leave one peptide out, and redo the analysis
% then the spaces should be similar.
% the way we'll test that is to see 
% a) first, are the spaces the same (i.e. are the vectors similar?)
% b) is each peptide similarly far from the origin?
% c) is each peptide similarly far from its neighbor
% we'll use the Euclidian distance (i.e. sqrt of the sum of squares) as a
% distance metric.

%OK, now I want to make my random data be in the princomp space I made:
for i = 1:120
 dex = circshift(1:120,[0 i]); % choose a unique set for each run
 [coefs scores variances] = princomp(normdata(dex(1:119),:)); % take the first six
 
 % how much variance is accounted for?
 tmp = 100*((variances)./sum(variances));
 varacct(i) = sum(tmp(1:10)); 
 
 % how similar are the vectors that define the space?
 eigvecproj(i,:) = diag(coefs(:,1:10)'*Fingerprint.coefs(:,1:10));
 
 % save the indiviudal projections
 skiponev(i,:,:) = coefs;
 
 % get the new projections
 proj(:,:,i) = scores(:,1:10);
 
 % now, how far is each point from the origin in this new space?
 veclen(i,:) = sqrt(proj(:,1,i).^2 + proj(:,2,i).^2 + proj(:,3,i).^2+proj(:,4,i).^2 + proj(:,5,i).^2 + proj(:,6,i).^2+proj(:,7,i).^2 + proj(:,8,i).^2 + proj(:,9,i).^2+proj(:,10,i).^2);
 
 % and how far is each point from each others?
 spacing(i,:) = pdist(proj(:,:,i));
end

% on the fourth run, eigenvectors 2 and 3 end up switched. so we'll
% fix that manually
%if size(data,2) == 16
%%eigvecproj(4,2) = skiponev(4,:,2)*fullv(:,3);
%eigvecproj(4,3) = skiponev(4,:,3)*fullv(:,2);
%end

%if size(data,2) == 17
%eigvecproj(5,2) = skiponev(5,:,2)*fullv(:,3);
%eigvecproj(5,3) = skiponev(5,:,3)*fullv(:,2);
%end

% now, how did the first three do?
sprintf('The first three components in the leave one out analysis\naccount for %g of the variance, SD = %g',mean(varacct),std(varacct))


numruns = 10000;
randveclen = zeros(numruns,120);
randspacing = zeros(numruns,7140);
randeigvecproj = zeros(numruns,10);

for i = 1:numruns
  randnormdata = reshape(normdata(randperm(numel(normdata))),[size(data,1) size(data,2)]);
  [coefs scores variances] = princomp(randnormdata);
  randeigvecproj(i,:) = diag(coefs(:,1:10)'*Fingerprint.coefs(:,1:10));
  randproj = scores(:,1:10);
  randveclen(i,:) = sqrt(randproj(:,1).^2 + randproj(:,2).^2 + randproj(:,3).^2+randproj(:,4).^2 + randproj(:,5).^2 + randproj(:,6).^2+randproj(:,7).^2 + randproj(:,8).^2 + randproj(:,9).^2+randproj(:,10).^2);
  % and how far is each point from each others?
  randspacing(i,:) = pdist(randproj);
end

% For random data, we'll just shuffle the original data so
% that we use the same values, but without any structure.

% first, the projections
bins = 0:.01:1;
P = hist(abs(eigvecproj(:)),bins);
randP = hist(abs(randeigvecproj(:)),bins);
P = P./sum(P);
randP = randP./sum(randP);

figure
hold on
plot(bins,cumsum(P),'k','linewidth',2);
plot(bins,cumsum(randP),'color',[.6 .6 .6],'linewidth',2)
axis square
set(gca,'xlim',[0 1],'ylim',[0 1], 'fontsize',12,'xtick',[0 .5 1],'ytick',[0 .5 1])
xlabel('Similarity of Eigenvectors')
ylabel('Cumulative Probability')
stevify 
print(gcf,'-depsc2','woods4.eps')

% that seems pretty clear; the vectors are largely the same even when we've
% left out one peptide. now, is each peptide about the same distance from
% the origin as it would have been in the original space?
bins = 0:.1:5;
randy = hist(randveclen(:),bins);
y = hist(veclen(:),bins);

figure
hold on
plot(bins,cumsum(randy./sum(randy)),'color',[.6 .6 .6],'linewidth',2)
plot(bins,cumsum(y./sum(y)),'color','k','linewidth',2)
axis square
box off
set(gca,'xlim',[0 5],'xtick',[0  5],'ylim',[0 1],'ytick',[0 .5 1],'fontsize',12)
xlabel('Length')
ylabel('Cumulative Probability')
stevify
print(gcf,'-depsc2','woods5.eps')


% looks pretty good too; now, what about the inter-point distance?
% here, we want to see how each value compares to the actual distance in
% the full space

% first, the random values
randx = repmat(pdist(fullproj),[numruns 1]);
x = repmat(pdist(fullproj),[120 1]);

% there's really no reason to plot all 10000 points
dex = randperm(numel(randx));
dex = dex(1:1000);

figure
hold on
plot(randx(dex),randspacing(dex),'.','color',[.6 .6 .6])
plot(x(:,1:7021),spacing(:,:),'k.')
plot([0 10],[0 10],'k:')
axis square
box off
set(gca,'xlim',[0 10],'ylim',[0 8],'xtick',[0 10],'ytick',[0 8],'fontsize',12)
xlabel('True Distance')
ylabel('Distance')
stevify
print(gcf,'-depsc2','woods6.eps')

% Plotting in the SVD defined space:
for i=1:20
DKOPair2(i,:)=[fullproj(i,1) fullproj(20+i,1) fullproj(40+i,1)];
DKOPair3(i,:)=[fullproj(i,4) fullproj(20+i,4) fullproj(40+i,4)];
WTPair3(i,:)=[fullproj(60+i,4) fullproj(80+i,4) fullproj(100+i,4)];
WTPair2(i,:)=[fullproj(60+i,1) fullproj(80+i,1) fullproj(100+i,1)];
end
x=[-5 4]
x2=[0 0]
y=[-3 4]
y2=[0 0]

figure; hold on
plot(fullproj(61:120,1),fullproj(61:120,4),'.','Color','blue')
plot(fullproj(1:60,1),fullproj(1:60,4),'.','Color','red')
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
plot(DKOPair2(20,:),DKOPair3(20,:),'Marker','o','Color','m')
plot(WTPair2(20,:),WTPair3(20,:),'Marker','o','Color','black')
plot(x,x2,'LineStyle','-','Color','black')
plot(y2,y,'LineStyle','-','Color','black')
axis([-11 12 -5 14])


