clear all;
close all;


allFuncs = {'meanSquaredDisplacement', 'meanSquaredDisplacementWeighted', 'meanSquaredDisplacementWeighted_1'};

marky = {'^','x','o'};
markycols = {0.3*[1 1 1], 0*[1 1 1], 0.4*[1 1 1]}
markySize = [5 5 5];

for f = 1:length(allFuncs)
	funcName = allFuncs{f};

	matfiles = dir(['EpiMats/' funcName '/*.mat']);

	gapvec{f} = [];
	rsq{f} = [];
	rsqShuffle{f} = [];
	singularRadius{f} = [];
	for m = 1:length(matfiles)
		matname = matfiles(m).name;
		dat = load(['EpiMats/' funcName '/' matname]);

		matrixOfOnes = ones(size(dat.realSecondOrder));
		diagmat = diag(diag(matrixOfOnes));
		diagel = find(diagmat>0);
		matrixOfOnes(diagel) = 0;

		uppermat = triu(matrixOfOnes);
		uppervals = find(uppermat > 0);

		gapvec{f} = [gapvec{f} 10^dat.thisgap];

		[matGuess singularRadius{f}(m)] = rank1Reconst(dat.realSecondOrder);

		gval = matGuess(uppervals);
		rval = dat.realSecondOrder(uppervals);

		rsq{f} = [rsq{f} errorFunc(rval,gval)];


		% Randomise

		epiMatRandom = zeros(size(dat.realSecondOrder));
		epiMatRandom(uppervals) = rval(randperm(length(uppervals)));
		epiMatRandom = epiMatRandom + epiMatRandom';
		matGuessRandom = rank1Reconst(epiMatRandom);
		rsqShuffle{f}(m) = errorFunc(epiMatRandom(uppervals),matGuessRandom(uppervals));

		disp(num2str(m))
	end
end

singradRand = [];
for randy = 1:100
	randmat = normrnd(0, 1,length(dat.realSecondOrder), length(dat.realSecondOrder));
	randmat = randmat + randmat';
	randmat(diagel) = 0;
	s = svd(randmat);
	singradRand(randy) = s(1)/sum(s);
end


mkdir('Panels')

figure;
hold on
for f = 1:length(allFuncs)
	% plot(gapvec{f},rsq{f},marky{f},'MarkerEdgeColor',markycols{f},'MarkerFaceColor',markycols{f},'MarkerSize',markySize(f),'LineWidth',1)
	plot(gapvec{f},rsq{f},marky{f},'MarkerEdgeColor',markycols{f},'MarkerSize',markySize(f),'LineWidth',1)
end
hold off
set(gca,'XScale','log')
pbaspect([3 1 1])
axis([.8 130 0 0.65])
makePretty
xticks([1 10 100])
yticks([0 0.25 0.5])
saveas(gcf,'Panels/punchline.png')
saveas(gcf,'Panels/svg_punchline.svg')
% close



function [matGuess singularRad] = rank1Reconst(epiMat)
	[U S V] = svd(epiMat);
	S = diag(S);
	singularRad = S(1)/sum(S);
	S(2:end) = 0;
	S = diag(S);
	matGuess = U*S*V';
end

function errVal = errorFunc(rval,gval)
	err = rval - gval;
	errVal = median(abs(err))/mean(abs(rval));
	% errVal = median(abs(err));
	% errVal = mean((rval - mean(rval)).*(gval - mean(gval)))/(std(rval)*std(gval));
end

function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',12)
	box on
end