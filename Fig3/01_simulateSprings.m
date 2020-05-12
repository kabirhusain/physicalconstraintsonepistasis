clear all;
close all;

%%%%% Parameters

rng('shuffle')

% Lattice parameters
lat.wid = 9;
lat.len = 9;
lat.spacing = 1; 
lat.rando = 0.3; % how much the initial lattice should be randomised

% Physics parameters
springy.springk = 1;
springy.maxforce = 0.0001; % Relaxation Settings -- total force on the system
springy.dt = 0.1;

% Plot Settings
plotty.colourStretched = [27,158,119]/255;
plotty.colourCompressed = [217,95,2]/255;
plotty.minStrain = .01; % anything more than 1% strain, and it's plotted as strained

% Settings for `evolving' a low mode
evol.stepsForLowModeEvolution = 60000;
evol.desiredModeGap = 1.2; % units of log 10
evol.desiredIPR = 2/(lat.wid*lat.len);
evol.shift = .025;

% `Genetic' settings
geno.MTeffect = 0.1; % fold-change in bond-lengths upon mutation.
geno.epistaticOrder = 2; % how many residues to mutate 

%%%%% The Real Deal

masterdirname = ['RawDat/Scan_Gap_' datestr(now)]
mkdir(masterdirname);
% Make a spring system
[lat.xy0 lat.adjacencyMatrix nonContactList] = perturbedHexagonalLattice(lat);
[springy.springMatrix springy.restLengthMatrix] = adjacencyToSprings(lat,springy,lat.xy0);


% geno.numberOfResidues = length(nonContactList);
% geno.residuesToMutate = nonContactList;

geno.numberOfResidues = lat.wid*lat.len;
geno.residuesToMutate = 1:geno.numberOfResidues;

save([masterdirname '/initCond.mat'])


gapvec = [0:.2:2];
for gp = 1:length(gapvec)
	clearvars -except gp masterdirname gapvec
	load([masterdirname '/initCond.mat']);
	evol.desiredModeGap = gapvec(gp) + 2*0.05*(rand-0.5);

	disp(['Desired Gap: ' num2str(evol.desiredModeGap)]);


	% Evolve le low mode
	[xyGeneric xyEnd gapvecEvol] = evolveLowModeByGeometry(lat,springy,evol);
	% springyGeneric = springy;
	springyEnd = springy;
	% [springyGeneric.springMatrix springyGeneric.restLengthMatrix] = adjacencyToSprings(lat,springyGeneric,xyGeneric);
	[springyEnd.springMatrix springyEnd.restLengthMatrix] = adjacencyToSprings(lat,springyEnd,xyEnd);

	% Compute their epistatic terms
	disp('Computing Epistatic Terms for (hopefully) Low Mode Network...')
	lat.xy0 = xyEnd;
	genoEnd = exhaustiveMutagenesis(lat,springyEnd,geno);


	dirname = [masterdirname '/Gap_' num2str(gapvec(gp))];
	mkdir(dirname);
	copyfile('linearreaponse_springymise.m',[dirname '/'])
	save([dirname '/data.mat']);
end








%%%%%%% Physics

%% Energy and force functions
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

%% Modes
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

function [energy minusforces] = energyAndForceToMin(xysplice,springy)
	xvals = xysplice(1:end/2);
	yvals = xysplice(1+end/2:end);

	energy = computeEnergy(springy,[xvals,yvals]);

	[forcex forcey] = computeForce(springy,[xvals,yvals]);
	minusforces = [-forcex; -forcey];
	minusforces = minusforces(:);
end

function [xy] = relaxSpringsComplex(springy,xy)
	dt = springy.dt;
	maxforce = springy.maxforce;
	springMatrix = springy.springMatrix;
	restLengthMatrix = springy.restLengthMatrix;

	[forcex forcey] = computeForce(springy, xy);

	thisEnergyAndForce = @(xysplice) energyAndForceToMin(xysplice,springy);

	xy0 = xy(:);
	% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off','StepTolerance',1e-14,'MaxFunctionEvaluations',1000000,'OptimalityTolerance',1e-8);
	options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'Display','off','StepTolerance',1e-10,'MaxFunctionEvaluations',1000000,'OptimalityTolerance',1e-9,'FunctionTolerance',1e-9);
	xy1 = fminunc(thisEnergyAndForce,xy0,options);
	x0 = xy1(1:end/2);
	y0 = xy1(1+end/2:end);

	xy = [x0 y0];
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
	% keyboard;
	while maxforcemag > maxforce
		xy = xy + dt*[forcex forcey];
		[forcex forcey] = computeForce(springy, xy);
		maxforcemag = sum(sqrt(forcex.^2 + forcey.^2));
		% disp([num2str(computeEnergy(springMatrix, restLengthMatrix,xy))])
		% keyboard
	end
	xy = xy + dt*[forcex forcey];

	% en = computeEnergy(springy,xy);
	% disp(['Final Energy: ' num2str(en) ' Magnitude of largest force: ' num2str(max(sqrt(forcex.^2 + forcey.^2)))]);
end













%%%%%%% Initial Conditions
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





%%%%%%% Protocols


% Do upto epistaticOrder of every residue
function geno = exhaustiveMutagenesis(lat,springy,geno)
	N = length(lat.adjacencyMatrix);
	epistaticOrder = geno.epistaticOrder;
	numberOfResidues = geno.numberOfResidues;

	% Relax the wild type.
	disp('Relaxing Wild Type...')
	% [springy.springMatrix springy.restLengthMatrix] = adjacencyToSprings(lat,springy,lat.xy0);
	xy_wt = lat.xy0;
	springyWT = springy;
	% Compute the lowest mode of the wt
	dynMatrix = dynamicalMatrix(springy, xy_wt);
	[modes modeenergies] = eig(dynMatrix);
	[modeenergies sorty] = sort(diag(modeenergies));
	lowModeVec = modes(:,sorty(4));
	invy = pinv(dynMatrix);

	% Now go through epistatic order and sort out what you need
	disp(['Exhaustively Mutagenising ' num2str(numberOfResidues) ' Bonds upto Order ' num2str(epistaticOrder)])
	

	geno.MTeffectvec = 1 +  geno.MTeffect*(1 - 2*(rand(length(springy.springMatrix),1) > 0.5));

	n = 0;
	for ord = 1:epistaticOrder
		% Do all ord-mutants
		allmut{ord} = nchoosek([1:numberOfResidues],ord);
		ModeDisplacement{ord} = containers.Map;
		energies{ord} = zeros(size(allmut{ord},1),1);
		% keyboard
		for mut = 1:size(allmut{ord},1) %length(allmut{ord})
			whichmut = allmut{ord}(mut,:); % this is a row vector of sites to mutate
			sitestomut = geno.residuesToMutate(whichmut);

			springyMT = springyWT;

			% mutate and compute energy
			[springyMT.springMatrix springyMT.restLengthMatrix] = mutationToSprings(geno, springyWT, sitestomut);
			
			[fx fy] = computeForce(springyMT, xy_wt);
			
			f = zeros(length(fx)*2,1);
			f(1:2:end) = fx;
			f(2:2:end) = fy;
			xyguess = invy*f;
			xyguess = [xyguess(1:2:end) xyguess(2:2:end)];
			xyguess = xy_wt + xyguess;
			
			% xy = relaxSprings(springyMT,xyguess);
			xy = xyguess;
			
			% Compute displacement from WT
			dispy = (xy(:) - xy_wt(:));

			ModeDisplacement{ord}(num2str(sort(whichmut))) = xy(:) - xy_wt(:);

			fprintf(repmat('\b',1,n));
			msg = [ 'Finished Mutant ' howfar({[mut,size(allmut{ord},1)]}) ' at Order ' howfar({[ord,epistaticOrder]})];
			fprintf(msg);
			n = numel(msg);
		end
	end
	fprintf('\n');

	geno.ModeDisplacement = ModeDisplacement;
end











%%%%%% Plotting Functions

function f = plotTowerOfModes(springy,xy)
	% Dynamical Matrix
	dynMatrix = dynamicalMatrix(springy,xy);

	[modes modeenergies] = eig(dynMatrix);
	modeenergies = diag(modeenergies);

	modeenergies = sort(modeenergies);

	modeenergies = modeenergies(4:end);

	f=figure('visible','off');
	hold on;
	for i = 1:length(modeenergies)
		plot([1 2],[modeenergies(i) modeenergies(i)],'-','LineWidth',2,'Color',[117,112,179]/255)
		% 217,95,2
	end
	hold off;
	set(gca,'YScale','log')
	axis([0 3 min(modeenergies)/2 2*max(modeenergies)])
	ax1 = gca;
	ax1.XAxis.Visible = 'off'; 
	makePretty
	% keyboard
end

%%% Plot the network
function f=plotSpringNetwork(xy,springy,plotty)
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
	set(gca,'FontSize',12)
	set(gca,'FontWeight','bold')
end

function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end






% Can you evolve a low mode?
function [xyStart xyEnd gapvec] = evolveLowModeByGeometry(lat,springy,evol)
	xy = lat.xy0;

	% Now let's try to affect the mode structure
	dynMatrix = dynamicalMatrix(springy, xy);
	[modes modeenergies] = eig(dynMatrix);
	[modeenergies sorty] = sort(diag(modeenergies));
	modegap = modeenergies(5) - modeenergies(4);
	moderat = log10(modeenergies(5)/modeenergies(4));
	
	lmVec = modes(:,sorty(4));

	lmx = lmVec(1:2:end);
	lmy = lmVec(2:2:end);

	lm = sqrt(lmx.^2 + lmy.^2);
	lm = lm/sum(lm.^2); % normalise L2
	invpartratio = sum(lm.^4)/sum(lm.^2);
	costfun = -modegap + .1*invpartratio;
	
	disp(['Starting gap: ' num2str(moderat) ' and Part. Ratio: ' num2str(invpartratio)]);

	% save starting
	xyStart = xy;
	springyStart = springy;

	N = lat.wid*lat.len;
	numberofsteps = evol.stepsForLowModeEvolution;
	n = 0;
	gapvec = []; %zeros(numberofsteps,1);
	for i = 1:numberofsteps
		% pick a residue, move it
		whichresid = randi([1 N]);
		xyNew = xy;
		xyNew(whichresid,:) = xyNew(whichresid,:) + 2*evol.shift*(rand(1,2)-0.5);
		springyNew = springy;

		[springyNew.springMatrix springyNew.restLengthMatrix] = adjacencyToSprings(lat,springyNew,xyNew);
		
		dynMatrix = dynamicalMatrix(springyNew, xyNew);
		[modes modeenergies] = eig(dynMatrix);
		[modeenergies sorty] = sort(diag(modeenergies));
		
		lmVec = modes(:,sorty(4));
		lmx = lmVec(1:2:end);
		lmy = lmVec(2:2:end);
		lm = sqrt(lmx.^2 + lmy.^2);
		lm = lm/sum(lm.^2); % normalise L2
		invpartratioNew = sum(lm.^4)/sum(lm.^2);

		if modeenergies(1) < -.001 % this means we're not really at a minimum
			keyboard
		end
		modegapNew = modeenergies(5) - modeenergies(4);
		moderatNew = log10(modeenergies(5)/modeenergies(4));

		% costfun = -modegap + invpartratio;
		% costfunNew = -modegapNew + invpartratioNew;

		if (moderatNew > moderat && (invpartratioNew <= invpartratio || invpartratioNew < evol.desiredIPR))
			modegap = modegapNew;
			moderat = moderatNew;
			invpartratio = invpartratioNew;
			xy = xyNew;
			springy = springyNew;
		end

		% gapvec(i) = modegap;
		gapvec = [gapvec moderat];
		% disp([howfar({[i numberofsteps]}) ': Mode gap is ' num2str(modegap)]);
		fprintf(repmat('\b',1,n));
		msg = [howfar({[i numberofsteps]}) ': Mode gap is ' num2str(moderat) ': and Inv Part Ratio is ' num2str(invpartratio)];
		fprintf(msg);
		n = numel(msg);

		if(moderat >= evol.desiredModeGap)
			break;
		end
	end
	fprintf('\n');

	xyEnd = xy;
end

