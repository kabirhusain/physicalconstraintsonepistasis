clear all;
close all;

% First, hard-code in the WT values
WT_inp = 1759616;
WT_sel = 3041819;


%%%% First, build a table of just the biggest effect

%% Go through each position's single mutants. Extract the mutant with the biggest effect, as well as average effect
allthesingles = readtable('sglmuts.csv');
allthesingles.Properties.VariableNames = {'WT','Residue','Mut','Inp','Sel'};

resnum = allthesingles.Residue;
residueLabels = unique(resnum);

singlemuts_maxfitW = zeros(size(residueLabels));
singlemuts_avfitW = zeros(size(residueLabels));
for i = 1:length(residueLabels)
	% extract subtable of mutations
	whichrows = find(resnum == residueLabels(i));
	thesemuts = allthesingles(whichrows,:);

	thesemuts_inp = thesemuts.Inp;
	thesemuts_sel = thesemuts.Sel;

	thesemuts_fitW = (thesemuts_sel./thesemuts_inp)/(WT_sel/WT_inp);
	thesemuts_fitW(find(thesemuts_fitW < 0.01)) = 0.01;

	singlemuts_avfitW(i) = exp(mean(log(thesemuts_fitW)));
	% keyboard
	[maxval maxpos] = max((log(thesemuts_fitW)));

	% maxpos = find(thesemuts_fitW == median(thesemuts_fitW));
	% maxpos = 7;
	singlemuts_maxfitW(i) = thesemuts_fitW(maxpos);
	craziestAA{i} = thesemuts(maxpos,:).Mut{1};
end

% let's save some memory
clear allthesingles

% Build the matrix of double mutants. Hopefully it's all there -- if its not, whoops

allthedoubles = readtable('dblmuts.csv');
allthedoubles.Properties.VariableNames = {'WT1', 'Res1', 'Mut1', 'WT2', 'Res2', 'Mut2', 'Inp', 'Sel','W1','W2'};

averageEpistasis = zeros(length(residueLabels),length(residueLabels));
bigEpistasis = zeros(length(residueLabels),length(residueLabels));
doublemut_avfitW = zeros(length(residueLabels),length(residueLabels));
res1list = allthedoubles.Res1;
res2list = allthedoubles.Res2;
for i = 1:(length(residueLabels)-1)
	for j = i+1:length(residueLabels)
		% Need to figure out where this combination is. Either 1 is i and 2 is j or the other way around
		whichrows = [intersect(find(res1list == residueLabels(i)), find(res2list == residueLabels(j))); intersect(find(res2list == residueLabels(i)), find(res1list == residueLabels(j)))];

		thesemuts = allthedoubles(whichrows,:);

		thesemuts_inp = thesemuts.Inp;
		thesemuts_sel = thesemuts.Sel;

		W1vals = thesemuts.W1; W1vals(find(W1vals < 0.01)) = 0.01;
		W2vals = thesemuts.W2; W2vals(find(W2vals < 0.01)) = 0.01;

		thismuts_additiveModelW = log(W1vals) + log(W2vals);

		thesemuts_fitW = (thesemuts_sel./thesemuts_inp)/(WT_sel/WT_inp);
		thesemuts_fitW(find(thesemuts_fitW < 0.01)) = 0.01;

		doublemut_avfitW(i,j) = exp(mean(log(thesemuts_fitW)));

		averageEpistasis(i,j) = mean(log(thesemuts_fitW) - thismuts_additiveModelW);


		% Now find the entry with the max amino acid for i and j
		thosewith1i = intersect(find(res1list == residueLabels(i)), find(strcmp(allthedoubles.Mut1, craziestAA{i})) );
		thosewith2j = intersect(find(res2list == residueLabels(j)), find(strcmp(allthedoubles.Mut2, craziestAA{j})) );

		those1i2j = intersect(thosewith1i, thosewith2j);

		thosewith1j = intersect(find(res1list == residueLabels(j)), find(strcmp(allthedoubles.Mut1, craziestAA{j})) );
		thosewith2i = intersect(find(res2list == residueLabels(i)), find(strcmp(allthedoubles.Mut2, craziestAA{i})) );

		those1j2i = intersect(thosewith1j, thosewith2i);

		therowIwant = allthedoubles([those1j2i those1i2j],:);
		theDoubleW = max((therowIwant.Sel/therowIwant.Inp)/(WT_sel/WT_inp),0.01);

		bigEpistasis(i,j) = log(theDoubleW) - log(singlemuts_maxfitW(i)) - log(singlemuts_maxfitW(j));
	end
end

save('computedepistasis.mat');
load('computedepistasis.mat');

epiMat = averageEpistasis;
epiMat(find(isnan(epiMat))) = 0;
epiMat = epiMat + epiMat';

maxval = .25*max(abs(epiMat(:)));
% maxval = 1
minval = -1*maxval;
cmap = gentleRedToBlue(100);

dirname = 'WithFilt_AvEpi_Panels'

mkdir(dirname)

figure('visible','off');
imagesc(residueLabels,residueLabels,epiMat)
colorbar
colormap(cmap)
caxis([minval maxval])
pbaspect([1 1 1])
makePretty
set(gca,'XTick',[],'YTick',[])
saveas(gcf,[dirname '/epi.svg'])
saveas(gcf,[dirname '/epi.png'])
close

% figure('visible','off'); 
% plot(sum(epiMat.^2,1),'ko')
% makePretty

%%%%% Rank 3 Decomposition
for k = 2:4
	[U,S,V] = svd(epiMat);
	S = diag(S);
	svec = S;
	S(k:end) = 0;
	S = diag(S);
	matGuess = U*S*V';

	diagmat = diag(ones(length(epiMat),1));
	diagel = find(diagmat == 1);
	matGuess(diagel) = 0;

	figure('visible','off');
	imagesc(residueLabels,residueLabels,matGuess)
	colorbar
	colormap(cmap)
	caxis([minval maxval])
	makePretty
	% set(gca,'XTick',[],'YTick',[])
	pbaspect([1 1 1])
	saveas(gcf,[dirname '/epi_rank' num2str(k-1) '.svg'])
	saveas(gcf,[dirname '/epi_rank' num2str(k-1) '.png'])
	close
end

matrixOfOnes = ones(size(epiMat));

diagmat = diag(diag(matrixOfOnes));
diagel = find(diagmat>0);



uppervals = setdiff(find(triu(matrixOfOnes) > 0), diagel);
err = epiMat(uppervals) - matGuess(uppervals);
rsq = 1 - sum(err.^2)/sum( (matGuess(uppervals) - mean(matGuess(uppervals))).^2 );

matGuess(diagel) = 0;


%%%%%%%%%% Compare against structure
structMat = epiMat - matGuess;
% whichhighlight = [27 31 43];
% whichhighlight = [27 31 43 39 40 41 54];
whichhighlight = [];
maxval = .2*max(abs(structMat(:)));
minval = -1*maxval;

cmapgreys = invGreys(100);

% pdb = prepare_pdb(pdbread('1pga.pdb'));


pdb = pdbread(['1pga.pdb']);

alltheatoms = squeeze(struct2cell(pdb.Model.Atom));
atomNames = alltheatoms(2,:);
atomChains = alltheatoms(5,:);

whichChain = strmatch('A',atomChains, 'exact');

% Make a list of Calphas and their positions
calphas=strmatch(['CA'],atomNames,'exact');

calphas = intersect(calphas,whichChain);

xyz_alphas = [];
resNum = [];
for atnum = 1:length(calphas)
	thisatom = alltheatoms(:,calphas(atnum));

	resNum(atnum) = thisatom{6};
	xyz_alphas = [xyz_alphas; thisatom{8} thisatom{9} thisatom{10}];
end

% Fucking altLocs
[resNum uniqOnes] = unique(resNum);
xyz_alphas = xyz_alphas(uniqOnes, :);

xyz = [xyz_alphas];
springy.springMatrix = zeros(size(xyz,1));

% Now, find all C-alpha C-alpha distances. What's a quick way to do this?
aaDistMat = zeros(length(resNum));
for i = 1:length(resNum)
	for j = 1:length(resNum)
		aaDistMat(i,j) = norm(xyz_alphas(i,:) - xyz_alphas(j,:));
	end
end

contactMap = aaDistMat(2:end,2:end)>12;

figure('visible','off');
imagesc(residueLabels,residueLabels, contactMap)
makePretty
% colorbar
colormap(cmapgreys)
pbaspect([1 1 1])
caxis([0 1])

% saveas(gcf,'Panels/structureOnly.svg')
% close

whichpercentile = 50

hold on
scalefactor = 7/max(abs(structMat(:)));
for i = 1:length(structMat)
	for j = i+1:length(structMat)
		if abs(structMat(i,j)) > prctile(abs(structMat(:)),whichpercentile);
			if length( find(whichhighlight == residueLabels(i)) ) > 0 && length( find(whichhighlight == residueLabels(j)) ) > 0
				plot(residueLabels(j),residueLabels(i), 'o','MarkerSize', scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [217 95 2]/255)
				plot(residueLabels(i),residueLabels(j), 'o','MarkerSize', scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', [217 95 2]/255)
			else
				plot(residueLabels(j),residueLabels(i), 'o','MarkerSize', scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor',[0 0 0])
				plot(residueLabels(i),residueLabels(j), 'o','MarkerSize', scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor',[0 0 0])

				% scatter(residueLabels(j),residueLabels(i), 3*scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor',[0 0 0])
				% scatter(residueLabels(i),residueLabels(j), 3*scalefactor*abs(structMat(i,j)), 'MarkerFaceColor', 'k', 'MarkerEdgeColor',[0 0 0])
			end
		end
	end
end
hold off
set(gcf,'Renderer','painters')
saveas(gcf, [dirname '/contactWithPred.svg'])
saveas(gcf, [dirname '/contactWithPred.png'])
close



% How many are contacts?
percentil = 80;
pervec = 0.25:.25:99;

whichupper = find(triu(ones(length(epiMat)),1) > 0);
epivec = epiMat(whichupper);
contactOrNot = 1-contactMap(whichupper);
[dumdum sortind] = sort(abs(structMat(whichupper)));


contFrac = [];
contFracbelow = [];
for percentil = pervec
	lowerval = ceil(length(sortind)*percentil/100);
	whichindabove = sortind(lowerval:end);
	whichindbelow = sortind(1:(lowerval-1));
	contFrac = [contFrac mean(contactOrNot(whichindabove))];
	contFracbelow = [contFracbelow mean(contactOrNot(whichindbelow))];
end

figure('visible','off');
plot(pervec, contFrac, '-','LineWidth',2,'Color',[27,158,119]/255)
hold on;
plot([pervec(1) pervec(end)], [ mean(contactOrNot), mean(contactOrNot) ], '--', 'Color',[.5, .5,.5],'LineWidth',2)
hold off
pbaspect([1.5 1 1])
axis([pervec(1) pervec(end) 0.35 0.85])
makePretty
saveas(gcf,[dirname '/perectileContacts.png'])
saveas(gcf,[dirname '/perectileContacts.svg'])
close






function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',12)
	box on
	% set(gca,'FontWeight','bold')
end

function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end

function cmap = invGreys(spac)
	cmap = [27,158,119;255,255,255];

	cmap = cmap/255;

	c1 = cmap(:,1);
	c2 = cmap(:,2);
	c3 = cmap(:,3);
	lenvec = [0:(size(cmap,1)-1)]/(size(cmap,1)-1);

	len2 = linspace(0,1,spac);

	cmap = [interp1(lenvec,c1,len2)' interp1(lenvec,c2,len2)' interp1(lenvec,c3,len2)'];
end

function cmap = gentleRedToBlue(spac)
	cmap = [103,0,31
	178,24,43
	214,96,77
	244,165,130
	253,219,199
	247,247,247
	209,229,240
	146,197,222
	67,147,195
	33,102,172
	5,48,97];

	cmap = cmap/255;

	c1 = cmap(:,1);
	c2 = cmap(:,2);
	c3 = cmap(:,3);
	lenvec = [0:(size(cmap,1)-1)]/(size(cmap,1)-1);

	len2 = linspace(0,1,spac);

	cmap = [interp1(lenvec,c1,len2)' interp1(lenvec,c2,len2)' interp1(lenvec,c3,len2)'];
end









%%% Physics
function dynMatrix = dynamicalMatrix3D(springy, equilibriumConfiguration)
	springMatrix = springy.springMatrix;
	% restLengthMatrix = springy.restLengthMatrix;

	% Compute dynamical matrix. Ingredients needed: unit vectors
	N = size(springMatrix,2);

	% The dynamical matrix is going to be 3N x 3N. the ij'th entries are 
	dynMatrix = zeros(3*N,3*N);
	for i = 1:N
		summedmatrices = zeros(3,3);
		whichj = find(springMatrix(i,:)>0);
		for j = whichj
			nx = equilibriumConfiguration(i,1) - equilibriumConfiguration(j,1);
			ny = equilibriumConfiguration(i,2) - equilibriumConfiguration(j,2);
			nz = equilibriumConfiguration(i,3) - equilibriumConfiguration(j,3);
			% keyboard
			nmag = sqrt(nx^2 + ny^2 + nz^2);
			nx = nx/nmag;
			ny = ny/nmag;
			nz = nz/nmag;

			thismat = -springMatrix(i,j)*[nx^2 nx*ny nx*nz; nx*ny ny^2 ny*nz; nz*nx ny*nz nz^2];

			% if length(find(isnan(thismat))) > 0
			% 	keyboard
			% end

			summedmatrices = summedmatrices + thismat;
			dynMatrix((3*(i-1) + 1):3*i,(3*(j-1) + 1):3*j) = thismat;
		end
		dynMatrix((3*(i-1) + 1):3*i,(3*(i-1) + 1):3*i) = -summedmatrices;
	end
end


function [forcex forcey forcez] = computeForce3D(springy,xy)
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	xvals = xy(:,1);
	yvals = xy(:,2);
	zvals = xy(:,3);

	% compute distancematrix(i,j) = |r_i - r_j|
	deltax = xvals - xvals';
	deltay = yvals - yvals';
	deltaz = zvals - zvals';

	distanceMatrix = sqrt(deltax.*deltax + deltay.*deltay + deltaz.*deltaz);
	forcePrefactor = springMatrix .* (distanceMatrix - restLengthMatrix)./(distanceMatrix);
	forcePrefactor(isnan(forcePrefactor)) = 0;
	forceMatx = -forcePrefactor .* deltax ;
	forceMaty = -forcePrefactor .* deltay ;
	forceMatz = -forcePrefactor .* deltaz ;

	forcex = sum(forceMatx,2);
	forcey = sum(forceMaty,2);
	forcez = sum(forceMatz,2);
end

function [springy xyz] = createENM(pdb, enmParam)
	cutoff = enmParam.cutoff;
	backboneSpring = enmParam.backboneSpring;
	aaSpringK = enmParam.aaSpringK;

	alltheatoms = squeeze(struct2cell(pdb.Model.Atom));
	atomNames = alltheatoms(2,:);
	atomChains = alltheatoms(5,:);

	whichChain = strmatch('A',atomChains, 'exact');

	% Make a list of Calphas and their positions
	calphas=strmatch('CA',atomNames,'exact');

	calphas = intersect(calphas,whichChain);

	xyz_alphas = [];
	resNum = [];
	for atnum = 1:length(calphas)
		thisatom = alltheatoms(:,calphas(atnum));

		resNum(atnum) = thisatom{6};
		xyz_alphas = [xyz_alphas; thisatom{8} thisatom{9} thisatom{10}];
	end

	% Fucking altLocs
	[resNum uniqOnes] = unique(resNum);
	xyz_alphas = xyz_alphas(uniqOnes, :);

	xyz = [xyz_alphas];
	springy.springMatrix = zeros(size(xyz,1));

	% Now, find all C-alpha C-alpha distances. What's a quick way to do this?
	aaDistMat = zeros(length(resNum));
	for i = 1:length(resNum)
		for j = 1:length(resNum)
			aaDistMat(i,j) = norm(xyz_alphas(i,:) - xyz_alphas(j,:));
		end
	end
	aaSprings = zeros(size(aaDistMat));
	aaSprings(find(aaDistMat <= cutoff)) = aaSpringK;
	for i = 1:(length(aaSprings)-1)
		aaSprings(i,i+1) = backboneSpring;
		aaSprings(i+1,i) = backboneSpring;
	end

	% Build the spring system
	springy.springMatrix = [aaSprings];
	% Get rid of self-interactions
	diagmat = diag(ones(length(springy.springMatrix),1));
	diagel = find(diagmat == 1);
	springy.springMatrix(diagel) = 0;

	% Compute restlength matrix
	x = xyz(:,1);
	y = xyz(:,2);
	z = xyz(:,3);

	DistMat = sqrt((x - x').^2 +(y - y').^2 + (z - z').^2);
	whichzero = find(springy.springMatrix == 0);
	springy.restLengthMatrix = DistMat;
	springy.restLengthMatrix(whichzero) = 0;
end


%%%% Epistasis
% Going to compute the epistatic `matrix' from deformations. How? Do all single and double mutants (for simplicity -- make AA bigger), compute some function (hmm?)
function epimat = computeEpistasis(springy, xyWT)
	% How many residues?
	resNum = length(springy.springMatrix);
	geno.MTeffectvec = 1.1*ones(resNum,1);

	% Compute inverse of hessian
	dynMatrix = dynamicalMatrix3D(springy, xyWT);
	[v d] = eig(dynMatrix);
	d = diag(d);
	% d(7:9) = .01;
	% d([7 8 10:end]) = 0;
	d = diag(d);
	dynMatrix = v*d*(v');
	invy = pinv(dynMatrix);

	% keyboard

	xy = xyWT;
	fitWT = distanceMatrixFitness(xy, xyWT);
	forceWT = computeForce3D(springy, xyWT);

	singlemut = [];
	for i = 1:resNum
		% Mutate residue i, compute fitness effect of mutant
		springyMT = springy;
		[springyMT.springMatrix springyMT.restLengthMatrix] = mutationToSprings(geno, springyMT, i);

		[fx fy fz] = computeForce3D(springyMT, xyWT);
		% keyboard
		f = zeros(3*length(fx),1);
		f(1:3:end) = fx;
		f(2:3:end) = fy;
		f(3:3:end) = fz;

		% keyboard

		dispi = invy*f;
		dispi = [dispi(1:3:end) dispi(2:3:end) dispi(3:3:end)];

		xy = xyWT + dispi;
		singlemut(i) = distanceMatrixFitness(xy, xyWT);
	end

	% keyboard

	% Compute epistasis
	epimat = zeros(resNum,resNum);
	for i = 1:resNum
		for j = (i+1):resNum
			% Compute fitness effect of double mutant i,j -- then epistasis
			springyMT = springy;
			[springyMT.springMatrix springyMT.restLengthMatrix] = mutationToSprings(geno, springyMT, [i j]);

			[fx fy fz] = computeForce3D(springyMT, xyWT);
			% keyboard
			f = zeros(3*length(fx),1);
			f(1:3:end) = fx;
			f(2:3:end) = fy;
			f(3:3:end) = fz;

			dispi = invy*f;
			dispi = [dispi(1:3:end) dispi(2:3:end) dispi(3:3:end)];

			xy = xyWT + dispi;
			doublemut = distanceMatrixFitness(xy, xyWT);

			epimat(i,j) = doublemut - ( - fitWT + singlemut(i) + singlemut(j) );
		end
	end

	epimat = epimat + epimat';
end


function m = distanceMatrixFitness(xy, xy0)
	% Compute the distance matrix
	x = xy(:,1);
	y = xy(:,2);
	z = xy(:,3);
	x0 = xy0(:,1);
	y0 = xy0(:,2);
	z0 = xy0(:,3);


	DistMat = sqrt((x - x').^2 +(y - y').^2 + (z - z').^2);
	DistMat0 = sqrt((x0 - x0').^2 +(y0 - y0').^2 + (z0 - z0').^2);

	diagmat = diag(ones(length(DistMat0),1));
	DistMat0Norm = DistMat0;
	DistMat0Norm(find(diagmat > 0)) = 1;

	distMatDiff = (DistMat-DistMat0)./DistMat0Norm;
	m = -1*mean(distMatDiff(:).*distMatDiff(:));
end

function [springMatrix restLengthMatrix] = mutationToSprings(geno, springyWT, sitestomut)
	% WT spring matrices
	springMatrix = springyWT.springMatrix;
	restLengthMatrix = springyWT.restLengthMatrix;
	
	sitestomut = sort(sitestomut); % need it sorted to avoid double counting

	% keyboard
	muteffect = @(i,j) sum(double(ismember(sitestomut,i))) + sum(double(ismember(sitestomut,j)));
	ismut = @(i) sum(double(ismember(sitestomut,i)));
	for i = sitestomut
		whichj = find(springMatrix(i,:)>0);
		for j = whichj
			% remember that A is 1, B is 2
			% restLengthMatrix(i,j) = (geno.MTeffect^muteffect(i,j))*springyWT.restLengthMatrix(i,j);
			restLengthMatrix(i,j) = (geno.MTeffectvec(i)^ismut(i))*((geno.MTeffectvec(j)^ismut(j)))*springyWT.restLengthMatrix(i,j);
			restLengthMatrix(j,i) = restLengthMatrix(i,j);
		end
	end
end
