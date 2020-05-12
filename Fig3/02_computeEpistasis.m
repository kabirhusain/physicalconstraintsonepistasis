% Specifically for my random target shape

clear all;
close all;




cmap = gentleRedToBlue(100);


% funcName = 'meanSquaredDisplacement';
funcName = 'meanSquaredDisplacementWeighted_1';


dirnames = dir('RawDat/Scan*');


mkdir('EpiMats')
mkdir(['EpiMats/' funcName])

mkdir(['Analysis'])
mkdir(['Analysis/' funcName])

mkdir('EpiPlots')
mkdir(['EpiPlots/' funcName])

%%%% Stuff for fitness functions
% Want a symmetric random matrix -- for weighting distances
matrixOfPain = 100*exp(3*rand(81));
matrixOfPain = matrixOfPain + matrixOfPain';
vectorOfPain = rand(81,1);
% Want a random set of positions to compare against
lat.wid = 9;
lat.len = 9;
lat.spacing = 1;
lat.rando = 0.3;
[xyRand lat.adjacencyMatrix nonContactList] = perturbedHexagonalLattice(lat);

gappy = [];
rsq = [];
k = 1;
for dirno = 1:length(dirnames)
	gapnames = dir(['RawDat/' dirnames(dirno).name '/Gap*']);

	for gp = 1:length(gapnames)
		dirname = ['RawDat/' dirnames(dirno).name '/' gapnames(gp).name];
		if length(dir([dirname '/*.mat'])) > 0
			dat = load([dirname '/data.mat']);

			xyWT = dat.xyEnd;
			if strcmp(funcName, 'distanceMatrixFitness')
				funcFunc = @(xy) distanceMatrixFitness(xy,xyWT);
			elseif strcmp(funcName, 'radiusOfGyration')
				funcFunc = @(xy) radiusOfGyration(xy);
			elseif strcmp(funcName, 'distanceMatrixFitnessWeighted')
				funcFunc = @(xy) distanceMatrixFitnessWeighted(xy, xyWT, matrixOfPain);
			elseif strcmp(funcName, 'distanceMatrixFitnessRandoWeighted')
				funcFunc = @(xy) distanceMatrixFitnessWeighted(xy, xyRand, matrixOfPain);
			elseif strcmp(funcName, 'distanceMatrixFitnessRando')
				funcFunc = @(xy) distanceMatrixFitness(xy, xyRand);
			elseif strcmp(funcName, 'onlyTwoLowestModes')
				dynMatrix = dynamicalMatrix(dat.springyEnd, xyWT);

				[modes modeenergies] = eig(dynMatrix);
				[modeenergies sorty] = sort(diag(modeenergies));
				% Shuffle so its [x y]
				lmVec = modes(:,sorty(4)); % lowest mode
				nlmVec = modes(:,sorty(5)); % next lowest mode
				lmx = lmVec(1:2:end);
				lmy = lmVec(2:2:end);
				nlmx = nlmVec(1:2:end);
				nlmy = nlmVec(2:2:end);
				lmVec = [lmx lmy];
				nlmVec = [nlmx nlmy];

				funcFunc = @(xy) onlyTwoLowestModes(xy, xyWT, lmVec, nlmVec);
			elseif strcmp(funcName, 'meanSquaredDisplacement')
				funcFunc = @(xy) meanSquaredDisplacement(xy, xyWT);
			elseif strcmp(funcName, 'meanSquaredDisplacementRando')
					funcFunc = @(xy) meanSquaredDisplacement(xy, xyRand);
			elseif length(regexp(funcName, 'meanSquaredDisplacementWeighted')) > 0 
				funcFunc = @(xy) meanSquaredDisplacementWeighted(xy, xyWT, vectorOfPain);
			end
			
			plotty = dat.plotty;

			gappy = [gappy dat.gapvecEvol(end)];

			%%%%% Predict Epistatic Matrices
			[uppertri realSecondOrder realFirstOrder predLow predLowvec] = secondOrderEpistasis(dat.springyEnd,dat.xyEnd,dat.genoEnd,dat.lat,funcFunc);
			% keyboard

			% maxval = .5*max(abs(realSecondOrder(:)));
			maxval = 2;
			minval = -maxval;


			%%%%% Rank 1 Decomposition
			% realSecondOrder = predLow;
			[U,S,V] = svd(realSecondOrder);
			S = diag(S);
			svec{k} = S;
			S(2:end) = 0;
			S = diag(S);
			matGuess = U*S*V';

			[vec eigval] = eig(realSecondOrder);
			eigval = diag(eigval);
			eigvalvec{k} = eigval;

			matrixOfOnes = ones(size(realSecondOrder));

			diagmat = diag(diag(matrixOfOnes));
			diagel = find(diagmat>0);

			matrixOfOnes(diagel) = 0;
			uppermat = triu(matrixOfOnes);
			lowermat = tril(matrixOfOnes);
			uppervals = find(uppermat > 0);

			matGuess(diagel) = 0;

			% keyboard

			thisgap = dat.gapvecEvol(end);
			save(['EpiMats/' funcName '/' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '.mat'],'realSecondOrder','thisgap','realFirstOrder');

			%%%% Plot
			epiMat = realSecondOrder;
			figure('visible','off');
			imagesc(epiMat)
			colorbar
			colormap(cmap)
			caxis([minval maxval]);
			makePretty
			xticks([]);
			yticks([]);
			% title(['Gap ' num2str(dat.gapvecEvol(end))])
			pbaspect([1,1,1])
			saveas(gcf,['EpiPlots/' funcName '/' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '.png']);
			saveas(gcf,['EpiPlots/' funcName '/svg_' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '.svg']);
			close

			% %%%%% Plot
			epiMat = matGuess;
			figure('visible','off');
			imagesc(epiMat)
			colorbar
			colormap(cmap)
			caxis([minval maxval]);
			makePretty
			xticks([]);
			yticks([]);
			% title(['Gap ' num2str(dat.gapvecEvol(end))])
			pbaspect([1,1,1])
			saveas(gcf,['EpiPlots/' funcName '/' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '_Reconst.png']);
			saveas(gcf,['EpiPlots/' funcName '/svg_' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '_Reconst.svg']);
			close

			figure('visible','off');
			imagesc(realSecondOrder - matGuess)
			colorbar
			colormap(cmap)
			caxis([minval maxval]);
			makePretty
			xticks([]);
			yticks([]);
			% title(['Gap ' num2str(dat.gapvecEvol(end))])
			pbaspect([1,1,1])
			saveas(gcf,['EpiPlots/' funcName '/' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '_Difference.png']);
			saveas(gcf,['EpiPlots/' funcName '/svg_' num2str(dat.gapvecEvol(end),'%.3f') '_' dirnames(dirno).name '_Difference.svg']);
			close
			
			disp([howfar({[gp , length(gapnames)], [dirno, length(dirnames)]})]);
			k = k+1;
		end
	end
end

save(['Analysis/' funcName '/data.mat'])

% figure;
% plot(10.^gappy, rsq,'ko','MarkerSize',8)
% makePretty
% set(gca,'XScale','log')



















function [xy adjacencyMatrix chosenList] = perturbedHexagonalLattice(lat)
	spacing = lat.spacing;
	wid = lat.wid;
	len = lat.len;
	rando = lat.rando;

	xvals = [];
	yvals = [];
	pxvals = [];
	pyvals = [];
	for i = 1:wid
		for j = 1:len
			xvals = [xvals; rem(j,2)*spacing/2 + (i*spacing) + 2*rando*(rand-0.5)];
			yvals = [yvals; j*spacing*(sin(acos(0.5))) + 2*rando*(rand-0.5)];

			pxvals = [pxvals rem(j,2)*spacing/2 + (i*spacing)];
			pyvals = [pyvals j*spacing*(sin(acos(0.5)))];
		end
	end

	% Find neighbours - brute force
	adjacencyMatrix = zeros(wid*len);
	neib = [];
	boringNodes = [];
	for node = 1:wid*len
		x0 = pxvals(node);
		y0 = pyvals(node);
		disty = sqrt( (pxvals - x0).^2 + (pyvals - y0).^2);
		whichones = find(disty < 1.1*spacing); % nearest and next-nearest
		if length(whichones) < 4
			boringNodes = [boringNodes node];
		end
		adjacencyMatrix(node,whichones) = 1;
		neib = [neib length(whichones)-1];
	end
	adjacencyMatrix = adjacencyMatrix - eye(length(adjacencyMatrix));

	chosenlist = [];
	for i=1:2:lat.wid
		chosenlist = [chosenlist lat.len*(i-1)+[1:2:lat.len]];
	end

	% keyboard;

	chosenlist = setdiff(chosenlist, boringNodes);
	chosenList = chosenlist;

	xy = zeros(length(xvals),2);
	xy(:,1) = xvals;
	xy(:,2) = yvals;
end


function [springMatrix restLengthMatrix] = adjacencyToSprings(lat,springy,xy)
	adjacencyMatrix = lat.adjacencyMatrix;

	% Generate a uniform spring constant matrix
	springMatrix = springy.springk*adjacencyMatrix;
	% and add some noise:
	% springMatrix = springMatrix.*(0.8 + .4*rand(size(springMatrix)));
	% springMatrix = (springMatrix + springMatrix')/2;

	restLengthMatrix = zeros(size(springMatrix));
	N = length(restLengthMatrix);
	for i = 1:N
		xi = xy(i,1);
		yi = xy(i,2);
		whichj = intersect(find(adjacencyMatrix(i,:)>0),[i:N]); % the intersect is to avoid double counting
		for j = whichj
			xj = xy(j,1);
			yj = xy(j,2);
			% remember that A is 1, B is 2
			restLengthMatrix(i,j) = sqrt((xi-xj)^2 + (yi-yj)^2);
			restLengthMatrix(j,i) = restLengthMatrix(i,j);
		end
	end
end


function [uppertri realSecondOrder realfirstOrder predLow predFirstOrder] = secondOrderEpistasis(springy,xy,geno,lat,funcFunc)
	dynMatrix = dynamicalMatrix(springy, xy);
	[v d] = eig(dynMatrix);
	[d sorty] = sort(diag(d));
	lm = v(:,sorty(4));
	% keyboard
	invy = pinv(dynMatrix);
	forces = zeros(length(dynMatrix),geno.numberOfResidues);
	predictedLow = zeros(geno.numberOfResidues,1);
	for i = 1:geno.numberOfResidues
		springyMT = springy;
		[springyMT.springMatrix springyMT.restLengthMatrix] = mutationToSprings(geno, springyMT, geno.residuesToMutate(i));
		[fx fy] = computeForce(springyMT, xy);
		% keyboard
		f = zeros(size(forces,1),1);
		f(1:2:end) = fx;
		f(2:2:end) = fy;
		forces(:,i) = f;

		% keyboard
		actualdisp = geno.ModeDisplacement{1}(num2str([i]));
		actualdisp = [actualdisp(1:end/2)  actualdisp(1+end/2:end)];
		a = zeros(length(geno.ModeDisplacement{1}(num2str([i]))),1);
		% keyboard
		a(1:2:end) = actualdisp(:,1);
		a(2:2:end) = actualdisp(:,2);

		predictedLow(i) = (lm'*f)/d(4);
	end
	lowMode = [lm(1:2:end) lm(2:2:end)];

	realSecondOrder = zeros(geno.numberOfResidues);
	realfirstOrder = zeros(geno.numberOfResidues,1);
	predLow = zeros(geno.numberOfResidues);
	uppertri = zeros(geno.numberOfResidues);
	for i = 1:geno.numberOfResidues-1
		fi = forces(:,i);
		dispi = invy*fi;
		dispi = [dispi(1:2:end) dispi(2:2:end)];
		actualdispi = geno.ModeDisplacement{1}(num2str([i]));
		actualdispi = [actualdispi(1:end/2)  actualdispi(1+end/2:end)];

		realfirstOrder(i) = funcFunc(xy+dispi) - funcFunc(xy);

		for j = i+1:geno.numberOfResidues
			fj = forces(:,j);

			dispj = invy*fj;
			dispj = [dispj(1:2:end) dispj(2:2:end)];
			dispij = dispi+dispj;

			actualdispj = geno.ModeDisplacement{1}(num2str([j]));
			actualdispj = [actualdispj(1:end/2)  actualdispj(1+end/2:end)];
			

			actualdispij = geno.ModeDisplacement{2}(num2str([i j]));
			actualdispij = [actualdispij(1:end/2)  actualdispij(1+end/2:end)];

			newpos = xy+actualdispj;
			predpos = xy+dispj;
			

			% realSecondOrder(i,j) = funcFunc(xy + dispi + dispj) - ( -funcFunc(xy) + funcFunc(xy+dispi) + funcFunc(xy+dispj) );
			realSecondOrder(i,j) = funcFunc(xy + actualdispij) - ( - funcFunc(xy) + funcFunc(xy+actualdispi) + funcFunc(xy+actualdispj) );
			predLow(i,j) = funcFunc(xy + dispij) - ( - funcFunc(xy) + funcFunc(xy+dispi) + funcFunc(xy+dispj) );
			
			% figure;
			% hold on
			% % plot(dispij(:)-actualdispij(:),actualdispij(:),'x');
			% plot(dispj(:),actualdispj(:),'x');
			% plot(dispi(:),actualdispi(:),'x');
			% % plot([min(dispij(:)) max(dispij(:))] , [min(dispij(:)) max(dispij(:))])
			% hold off;
			% makePretty

			% keyboard

			% err = (dispij(:) - actualdispij(:));

			% figure;
			% hist(abs(err))
			% makePretty

			% keyboard

			% keyboard;
			predLowSumma = funcFunc(xy + (dispi + dispj) ) - ( - funcFunc(xy) + funcFunc(xy + dispi) + funcFunc(xy + dispj) );
			% keyboard
			uppertri(i,j) = 1;
		end
	end
	% keyboard
	normFactor = mean(abs(realfirstOrder));
	% keyboard
	realSecondOrder = (realSecondOrder + realSecondOrder')/normFactor;
	predLow = (predLow + predLow')/normFactor;

	predFirstOrder = predictedLow;
	% keyboard
end



function dynMatrix = dynamicalMatrix(springy, equilibriumConfiguration)
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	% Compute dynamical matrix. Ingredients needed: unit vectors, resting strain for each bond
	N = size(springMatrix,2);

	% The dynamical matrix is going to be 2N x 2N. the ij'th entries are 
	dynMatrix = zeros(2*N,2*N);
	for i = 1:N
		summedmatrices = zeros(2,2);
		whichj = find(springMatrix(i,:)>0);
		for j = whichj
			nx = equilibriumConfiguration(i,1) - equilibriumConfiguration(j,1);
			ny = equilibriumConfiguration(i,2) - equilibriumConfiguration(j,2);
			nmag = sqrt(nx^2 + ny^2);
			nx = nx/nmag;
			ny = ny/nmag;

			restStrain = (nmag - restLengthMatrix(i,j))/restLengthMatrix(i,j);

			thismat = -0.5*springMatrix(i,j)*[nx^2 nx*ny; nx*ny ny^2]/(restStrain+1);
			thismat = thismat - 0.5*springMatrix(i,j)*restStrain*eye(2)/(restStrain+1);

			summedmatrices = summedmatrices + thismat;
			dynMatrix(2*(i-1) + 1:2*i,2*(j-1) + 1:2*j) = thismat;
		end
		dynMatrix(2*(i-1) + 1:2*i,2*(i-1) + 1:2*i) = -summedmatrices;
	end
	dynMatrix = 2*dynMatrix;
end


function energy = computeEnergy(springy,xy);
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	xvals = xy(:,1);
	yvals = xy(:,2);

	% compute deltax and deltay
	deltax = xvals' - xvals;
	deltay = yvals' - yvals;

	distanceMatrix = sqrt(deltax.*deltax + deltay.*deltay);

	energymat = 0.5*springMatrix.*((distanceMatrix - restLengthMatrix).^2);

	energy = sum(energymat(:))/2; % divide by 2 to get rid of overcounting
end

function [forcex forcey] = computeForce(springy,xy)
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	xvals = xy(:,1);
	yvals = xy(:,2);

	% compute distancematrix(i,j) = |r_i - r_j|
	deltax = xvals - xvals';
	deltay = yvals - yvals';

	distanceMatrix = sqrt(deltax.*deltax + deltay.*deltay);
	forcePrefactor = springMatrix .* (distanceMatrix - restLengthMatrix)./(distanceMatrix);
	forcePrefactor(isnan(forcePrefactor)) = 0;
	forceMatx = -forcePrefactor .* deltax ;
	forceMaty = -forcePrefactor .* deltay ;

	forcex = sum(forceMatx,2);
	forcey = sum(forceMaty,2);
end

function [xy] = relaxSprings(springy,xy)
	dt = springy.dt;
	maxforce = springy.maxforce;
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	[forcex forcey] = computeForce(springy, xy);

	% en = computeEnergy(springy,xy);
	% disp(['Starting Energy: ' num2str(en) ' Magnitude of largest force: ' num2str(max(sqrt(forcex.^2 + forcey.^2)))]);

	maxforcemag = sum(sqrt(forcex.^2 + forcey.^2));
	
	while maxforcemag > maxforce
		xy = xy + dt*[forcex forcey];
		[forcex forcey] = computeForce(springy, xy);
		maxforcemag = sum(sqrt(forcex.^2 + forcey.^2));
	end
	xy = xy + dt*[forcex forcey];

	% en = computeEnergy(springy,xy);
	% disp(['Final Energy: ' num2str(en) ' Magnitude of largest force: ' num2str(max(sqrt(forcex.^2 + forcey.^2)))]);
end


function [springMatrix restLengthMatrix] = mutationToSprings(geno, springyWT, sitestomut)
	% WT spring matrices
	springMatrix = springyWT.springMatrix;
	restLengthMatrix = springyWT.restLengthMatrix;
	residuesToMutate = geno.residuesToMutate;

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






























%%% Plot the network
function f=plotSpringNetwork(geno,xy,springy,plotty)
	N = size(xy,1);
	f = figure('visible','off');
	axis off
	hold on

	for i = 1:N
		whichj = intersect(find(springy.springMatrix(i,:) > 0),[1:N]);
		for j = whichj
			resty = springy.restLengthMatrix(i,j);
			disty = sqrt((xy(i,1) - xy(j,1))^2 + (xy(i,2) - xy(j,2))^2);
			strain = disty - resty;
			isstrained = double(abs(strain/resty) >= plotty.minStrain);
			straincolor = double(strain > 0)*plotty.colourStretched + double(strain < 0)*plotty.colourCompressed;
			straincolor = (1-isstrained)*[0 0 0] + isstrained*straincolor;

			plot([xy(i,1) xy(j,1)], [xy(i,2) xy(j,2)], '-','LineWidth',2,'Color',straincolor);
		end
	end

	% plot(xy(aNodes,1),xy(aNodes,2),'o','MarkerSize',10,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0],'LineWidth',1);
	% plot(xy(bNodes,1),xy(bNodes,2),'o','MarkerSize',10,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],'LineWidth',1);

	% % Which ones do you mutate?
	% plot(xy(geno.residuesToMutate,1),xy(geno.residuesToMutate,2),'o','MarkerEdgeColor',[1 0 0],'LineWidth',3,'MarkerSize',10);
end


function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',15)
	% set(gca,'FontWeight','bold')
end

function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
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




function cmap = gentleBlues(spac)
	cmap = [247,252,240
	224,243,219
	204,235,197
	168,221,181
	123,204,196
	78,179,211
	43,140,190
	8,104,172
	8,64,129];

	cmap = cmap/255;

	c1 = cmap(:,1);
	c2 = cmap(:,2);
	c3 = cmap(:,3);
	lenvec = [0:(size(cmap,1)-1)]/(size(cmap,1)-1);

	len2 = linspace(0,1,spac);

	cmap = [interp1(lenvec,c1,len2)' interp1(lenvec,c2,len2)' interp1(lenvec,c3,len2)'];
end



































function m = radiusOfGyration(xy)
	meanx = mean(xy(:,1));
	meany = mean(xy(:,2));

	x = xy(:,1) - meanx;
	y = xy(:,2) - meany;

	m = mean(x.^2 + y.^2);
end

function m = findConvexHull(xy)
	[bah m] = convhull(xy(:,1),xy(:,2));
end

function m = findConvexHullSubset(xy,subset)
	x = xy(subset,1);
	y = xy(subset,2);
	[bah m] = convhull(x,y);
end

function m = distanceFromTarget(xy,targetxy)
	xvals = targetxy(:,1);
	yvals = targetxy(:,2);
	deltax = xvals' - xvals;
	deltay = yvals' - yvals;
	distanceMatrixTarget = sqrt(deltax.*deltax + deltay.*deltay);

	xvals = xy(:,1);
	yvals = xy(:,2);
	deltax = xvals' - xvals;
	deltay = yvals' - yvals;
	distanceMatrixActual = sqrt(deltax.*deltax + deltay.*deltay);

	m = (distanceMatrixActual - distanceMatrixTarget).^2;
	m = sqrt(sum(m(:)))/2;
	% keyboard
end

function m = projectionOntoTarget(xy,targetxy)
	xvals = targetxy(:,1);
	yvals = targetxy(:,2);
	deltaxTarget = xvals' - xvals;
	deltayTarget = yvals' - yvals;
	distanceMatrixTarget = sqrt(deltaxTarget.*deltaxTarget + deltayTarget.*deltayTarget);
	% deltarall
	deltaxTarget = deltaxTarget./(eps+distanceMatrixTarget);
	deltayTarget = deltayTarget./(eps+distanceMatrixTarget);

	xvals = xy(:,1);
	yvals = xy(:,2);
	deltax = xvals' - xvals;
	deltay = yvals' - yvals;
	distanceMatrixActual = sqrt(deltax.*deltax + deltay.*deltay);

	deltax = deltax./(eps+distanceMatrixActual);
	deltay = deltay./(eps+distanceMatrixActual);

	projyx = deltaxTarget.*deltax;
	projyt = deltayTarget.*deltay;

	projtot = projx + projy;

	nonz = find(abs(deltax) > 0);

	m = mean(projtot(nonz));
	% keyboard
end

function m = distTwo(xy,whichTwo)
	x1 = xy(whichTwo(1),1);
	y1 = xy(whichTwo(1),2);
	x2 = xy(whichTwo(2),1);
	y2 = xy(whichTwo(2),2);
	% [x1 y1] = [xy(whichTwo(1),1) xy(whichTwo(1),2)];
	% [x2 y2] = [xy(whichTwo(2),1) xy(whichTwo(2),2)];

	m = sqrt((x1-x2).^2 + (y1-y2).^2);
end

function m = distanceMatrixFitness(xy, xy0)
	% Compute the distance matrix
	x = xy(:,1);
	y = xy(:,2);
	x0 = xy0(:,1);
	y0 = xy0(:,2);


	DistMat = sqrt((x - x').^2 +(y - y').^2);
	DistMat0 = sqrt((x0 - x0').^2 +(y0 - y0').^2);

	diagmat = diag(ones(length(DistMat0),1));
	DistMat0Norm = DistMat0;
	DistMat0Norm(find(diagmat > 0)) = 1;

	distMatDiff = (DistMat-DistMat0)./DistMat0Norm;
	m = -1*mean(distMatDiff(:).*distMatDiff(:));
end

function m = distanceMatrixFitnessWeighted(xy, xy0, weightMatrix)
	% Compute the distance matrix
	x = xy(:,1);
	y = xy(:,2);
	x0 = xy0(:,1);
	y0 = xy0(:,2);


	DistMat = sqrt((x - x').^2 +(y - y').^2);
	DistMat0 = sqrt((x0 - x0').^2 +(y0 - y0').^2);

	diagmat = diag(ones(length(DistMat0),1));
	DistMat0Norm = DistMat0;
	DistMat0Norm(find(diagmat > 0)) = 1;

	distMatDiff = weightMatrix.*(DistMat-DistMat0)./DistMat0Norm;
	m = (-1*mean(distMatDiff(:).*distMatDiff(:)));
end

function m = meanSquaredDisplacement(xy, xy0)
	m = -1*mean(sum((xy-xy0).^2, 2));
end

function m = meanSquaredDisplacementWeighted(xy, xy0, weightVec)
	m = -1*mean(weightVec.*sum((xy-xy0).^2, 2));
end

function m = onlyTwoLowestModes(xy, xy0, mode1, mode2)
	% Find the deformation vector
	defo = xy - xy0;
	% Project onto the two lowest modes
	defo1 = sum(sum(defo.*mode1));
	defo2 = sum(sum(defo.*mode2));

	m = -1*defo1^2-5*defo2^2;
end