% glucose	50 mM NaCl	pH 6	42C	anaerobic	stationary 1 day	stationary 3 days	chemostat mu=0.12	chemostat mu=0.5	galactose	acetate	glycerol	LB	glucosamine	fumarate	succinate	pyruvate	chemostat mu=0.20	chemostat mu=0.35
clear all;
close all;

datmat = csvread('ecoli_copynumber.csv'); % we start with 2039 proteins

% remove stationary phase data
datmat = datmat(:,[1:5 8:19]);

% the first column is glucose -- everything needs to be normalised to that.
% First, find all 'below LOQ' entries those have value minus 1 -- remove those rows
whichLOQ = [];
for growthCond = 1:size(datmat,2)
	thisProt = datmat(:,growthCond);
	LOQThis = find(thisProt < 0);
	whichLOQ = [whichLOQ; LOQThis];
end
whichLOQ = unique(whichLOQ);
remainingCols = setdiff(1:size(datmat,1),whichLOQ);
datmat = datmat(remainingCols,:); % now we have 1968 proteins

glucProt = datmat(:,1);
galacProt = datmat(:,8);
acetateProt = datmat(:,9);
% First try just a scatter
figure;
plot(galacProt./glucProt, acetateProt./glucProt, 'ko','MarkerSize',4,'MarkerFaceColor','k')
makePretty
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'FontSize',22)
pbaspect([1 1.5 1])
axis([.001 100 .001 100])
xticks([.01 1 100])
box off
saveas(gcf,'FoldChange_GluctoGalac_GluctoAcetate.svg')

[coeff,score,latent,tsquared,explained,mu] = pca(datmat');
% randmat = rand(size(datmat));
randmat = datmat(:);
randmat = randmat(randperm(length(randmat)));
randmat = reshape(randmat, size(datmat,1), size(datmat,2));

[coeffRand,scoreRand,latentRand,tsquaredRand,explainedRand,muRand] = pca(randmat');

figure;
hold on
plot(cumsum(explainedRand),'ko','MarkerSize',9, 'MarkerFaceColor', [150 150 150]/255)
plot(cumsum(explained),'ko','MarkerSize',9, 'MarkerFaceColor', [217 95 2]/255)
hold off;
makePretty
axis([0 17 0 105])
pbaspect([1.5 1 1])
box on
saveas(gcf,'ecoli_VarianceExplained.svg')

figure;
plot(score(:,1),score(:,2),'ko','MarkerSize',8,'MarkerFaceColor','k')
makePretty
pbaspect([2 1 1])
% axis([])

function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',16)
end