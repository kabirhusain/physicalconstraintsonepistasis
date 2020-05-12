% Analyse networks simulated by makeSignallingNetworks.m
% Computes epistasis for three fitness functions as defined in SI Sec. IIIB.

clear all;
close all;


masterdatdir = 'Data';

% Specify the random matrices used to compute the fitness function. Let the first just be the identity matrix
matrixOfPain{1} = eye(30);
for i = 2:3
	randmat = diag(exp(normrnd(0,1,30,1)));
	matrixOfPain{i} = randmat;
end

dirlist = dir([masterdatdir '/*Networks_*']);

for i = 1:3
	gapvec{i} = [];
	errorVec{i} = [];
end
for dirno = 1:length(dirlist)
	matlist = dir([masterdatdir '/' dirlist(dirno).name '/*ModeGap*.mat']);

	for matno = 1:length(matlist)
		dat = load([masterdatdir '/' dirlist(dirno).name '/' matlist(matno).name]);

		epistatLow = dat.epistatLow;
		matrixOfOnes = ones(size(dat.epistatLow.epimat));
		diagmat = diag(diag(matrixOfOnes));
		diagel = find(diagmat>0);
		matrixOfOnes(diagel) = 0;
		uppermat = triu(matrixOfOnes);
		uppervals = find(uppermat > 0);

		for i = 1:3
			epistatLow = dat.epistatLow;
			epistatLow.matrixOfPain = matrixOfPain{i};
			epistatLow = redoTheEpiEdges(dat.signalLow, epistatLow, dat.n0Low);
			epiMat = epistatLow.epimat;

			matGuess = rank1Reconst(epiMat);
			gvalLow = matGuess(uppervals);
			rvalLow = epiMat(uppervals);

			% Error
			gapvec{i} = [gapvec{i} dat.signalLow.optimisationTimeSeries(end)];
			errorVec{i} = [errorVec{i} errorFunc(rvalLow, gvalLow)];
		end

		disp(howfar({[matno, length(matlist)], [dirno, length(dirlist)]}));
	end
end



save('tempdata.mat');


load('tempdata.mat');

marky = {'^','x','o'};
markycols = {0.2*[1 1 1], 0*[1 1 1], 0.4*[1 1 1]}
markysize = [5 5 5];



mkdir('EpiPanels')

figure
hold on;
for i = 1:3
	thosewanted{1} = 1:length(gapvec{i});

	for f = 1:length(thosewanted)
		plot(gapvec{i}(thosewanted{f}), errorVec{i}(thosewanted{f}), marky{i}, 'MarkerSize', markysize(f), 'MarkerEdgeColor', markycols{i}, 'LineWidth', 1);
	end
end
pbaspect([3 1 1])
axis([1.2 15 -0.025 0.45])
set(gca,'XScale','log')
makePretty
xticks([2:2:10])
% saveas(gcf,'EpiPanels/Punchline.svg')
% saveas(gcf,'EpiPanels/Punchline.png')


function matGuess = rank1Reconst(epiMat)
	[U S V] = svd(epiMat);
	S = diag(S);
	S(2:end) = 0;
	S = diag(S);
	matGuess = U*S*V';
end

function errVal = errorFunc(rval,gval)
	err = rval - gval;
	errVal = median(abs(err))/mean(abs(rval));
end


%%%%% Basic Functions - [F1]
function signal = generateKij(signal)
	% Decide whether links are regulatory or not
	signal.regulatoryLinks = full(sprand(signal.internalNetwork,signal.internalNetwork,signal.internalSparsity)) > 0;
	
	% No self regulation -- kill terms on diagonal
	diagmat = diag(ones(length(signal.regulatoryLinks),1));
	diagel = find(diagmat > 0);
	signal.regulatoryLinks(diagel) = 0;

	poslinks = find(signal.regulatoryLinks > 0);
	

	signal.kcatij = zeros(size(signal.regulatoryLinks));
	signal.kcatij(poslinks) = rand(size(poslinks));
	% signal.kcatij(poslinks) = 1;

	signal.Kij = zeros(size(signal.regulatoryLinks));
	signal.Kij(poslinks) = rand(size(poslinks));
	% signal.Kij(poslinks) = 1;	
end


function n = runToSteadyState(n,signal)
	t = 0;
	notdone = 1;
	while notdone
		thist = 0;
		lastno = n;
		while thist < 1
			[t n] = goForth(t,n,signal);
			thist = thist + signal.dt;
		end
		newno = n;
		if abs(lastno - newno)/lastno < signal.reqAccu
			notdone = 0;
		end

		if t > 100000
			keyboard
		end
	end
end

function [t n] = goForth(t,n,signal)
	dt = signal.dt;

	nmat = repmat(n,[1 signal.internalNetwork]);

	% fn = (1-n).*( sum(transpose(signal.kcatij.*nmat./(signal.Kij + nmat)), 2) ) - signal.decayRate*n;
	fn = ( sum(transpose(signal.kcatij.*nmat./(signal.Kij + nmat)), 2) ) - signal.decayRate*n;

	n = n + fn*dt;
	t = t + dt;
end

function jacob = generateJacobian(signal,n)
	jacob = zeros(signal.internalNetwork);

	nmat = repmat(n, [1 signal.internalNetwork]);
	jacob = transpose((signal.kcatij./(signal.Kij + nmat)).*( 1 - nmat./(signal.Kij + nmat) ) );

	diagmat = diag(ones(length(jacob),1));
	diagel = find(diagmat == 1);
	jacob(diagel) = - signal.decayRate;
end


%%%%% Double-check your jacobian calculation
function fn = forceFelt(n,signal)
	nmat = repmat(n,[1 signal.internalNetwork]);

	fn = ( sum(transpose(signal.kcatij.*nmat./(signal.Kij + nmat)), 2) ) - signal.decayRate*n;
end


function numJacob = numericallyGenerateJacobian(signal,n)
	numJacob = zeros(length(n));
	for i = 1:length(n)
		n2 = n;
		n2(i) = n(i) + 0.0001*n(i);

		fn = forceFelt(n2,signal);
		numJacob(:,i) = fn/(0.0001*n(i));
	end
end



function [n0 signal] = perturbKij(evosize,signal,n0)
	% Enzymes
	enzymeLabels = find(abs(signal.regulatoryLinks) > 0); % what are all of the enzymes?

	whichenz = enzymeLabels(randi(length(enzymeLabels)));

	signal.kcatij(whichenz) = signal.kcatij(whichenz) + 2*evosize*(rand - 0.5);
	signal.kcatij(whichenz) = max(min(signal.kcatij(whichenz),1),0);

	signal.Kij(whichenz) = signal.Kij(whichenz) + 2*evosize*(rand - 0.5);
	signal.Kij(whichenz) = max(min(signal.Kij(whichenz),1),0);

	n0 = runToSteadyState(n0,signal);
end





function f = plotTowerOfModes(jacob)
	[v d] = eig(jacob);
	d = abs(real(diag(d)));
	modeenergies = sort(d);

	f=figure('visible','off');
	hold on;
	plot([1 2],[modeenergies(1) modeenergies(1)],'-','LineWidth',3,'Color',[217,95,2]/255)
	for i = 2:length(modeenergies)
		plot([1 2],[modeenergies(i) modeenergies(i)],'-','LineWidth',2,'Color',[0,0,0]/255)
		% 217,95,2
	end
	hold off;
	set(gca,'YScale','log')
	axis([0 3 min(modeenergies)/2 2*max(modeenergies)])
	ax1 = gca;
	ax1.XAxis.Visible = 'off'; 
	makePretty
end














%%%%% Protocols

function [n0 signal] = gimmeALowMode(signal,lowMode)
	signalWT = signal;


	n0 = 0.5*ones(signal.internalNetwork,1);

	% Find the steady state and the relaxational mode
	n0 = runToSteadyState(n0, signal);
	jacob = generateJacobian(signal,n0);
	[v d] = eig(jacob);
	d = abs(real(diag(d)));
	ds = sort(d);
	moderat = ds(2)/ds(1);

	minone = find(d==min(d));
	minone = minone(1);
	lm = v(:,minone);
	invpartratio = sum((lm.*conj(lm)).^2,1);

	evosize = lowMode.evosize;
	desiredIPR = lowMode.desiredIPR;
	numberofsteps = lowMode.numberofsteps;
	n = 0; 
	moderatvec = zeros(numberofsteps,1);

	tBeta = 50;

	disp(['Starting Gap is: ' num2str(moderat)])
	for i = 1:numberofsteps
		% Generate a perturbed network
		signalNew = signal;

		[n0New signalNew] = perturbKij(evosize,signal,n0);
		
		jacob = generateJacobian(signalNew,n0New);
		[v d] = eig(jacob);
		% if max(real(diag(d))) > 0
		% 	disp('whoops')
		% 	keyboard
		% end
		d = abs(real(diag(d)));
		ds = sort(d);
		moderatNew = ds(2)/ds(1);
		minone = find(d==min(d));
		minone = minone(1);
		lm = v(:,minone);
		invpartratioNew = sum((lm.*conj(lm)).^2,1);


		% how's it going?
		if (invpartratioNew <= invpartratio || invpartratioNew < desiredIPR)
			if moderatNew > moderat
				signal = signalNew;
				n0 = n0New;
				moderat = moderatNew;
				invpartratio = invpartratioNew;
			end
		end

		
		% fprintf(repmat('\b',1,n));
		msg = [howfar({[i numberofsteps]}) ': Mode gap is ' num2str(moderat) ': and Inv Part Ratio is ' num2str(invpartratio) ];
		% fprintf(msg);
		n = numel(msg);
		moderatvec(i) = moderat;

		if moderat >= lowMode.desiredRat
			moderatvec = moderatvec(1:i);
			break;
		end
	end

	disp(['Final Gap is: ' num2str(moderat)])


	signal.optimisationTimeSeries = moderatvec;
end


% Mutate each edge
function epistat = doTheEpiEdges(signal,epistat, n0)
	% WT
	signalWT = signal;
	n0_WT = n0;

	% Enzymes
	allEnzymeLabels = find(signal.Kij > 0);
	if epistat.mutNumber == 0
		% mutate all
		enzymeLabels = allEnzymeLabels;
		enzymeNumber = length(enzymeLabels);
	else
		% find the biggest kcats
		nonzerocats = find(signal.kcatij > 0);

		nmat = repmat(n0,[1 signal.internalNetwork]);

		fn = signal.kcatij.*nmat./(signal.Kij + nmat);

		[a kcatid] = sort(fn(nonzerocats),'desc');

		enzymeLabels = nonzerocats(kcatid(1:epistat.mutNumber));
		enzymeNumber = length(enzymeLabels);
		% keyboard
	end

	mutEffect = 1 + epistat.mutationMagnitude*(-2*(rand(size(enzymeLabels))>0.5) + 1);
	mutEffect2 = 1 + epistat.mutationMagnitude*(-2*(rand(size(enzymeLabels))>0.5) + 1);

	fitval = exp(-0.5*transpose(n0_WT-n0_WT)*epistat.matrixOfPain*(n0_WT-n0_WT));
	fitvalWT = fitval;
	epivec = zeros(enzymeNumber,1);
	epimat = zeros(enzymeNumber);
	for i = 1:enzymeNumber
		% Mutate
		signal = signalWT;
		signal.kcatij(enzymeLabels(i)) = mutEffect(i)*signalWT.kcatij(enzymeLabels(i));
		signal.Kij(enzymeLabels(i)) = mutEffect2(i)*signalWT.Kij(enzymeLabels(i));
		
		n0 = runToSteadyState(n0, signal);

		epistat.singleMut{i} = n0-n0_WT;

		fitval = exp(-0.5*transpose(n0-n0_WT)*epistat.matrixOfPain*(n0-n0_WT));;

		epivec(i) = fitval - fitvalWT;
		% disp(howfar({[i enzymeNumber]}));
	end

	for i = 1:(enzymeNumber-1)
		for j = (i + 1):enzymeNumber
			% Mutate
			signal = signalWT;
			signal.kcatij(enzymeLabels(i)) = mutEffect(i)*signalWT.kcatij(enzymeLabels(i));
			signal.kcatij(enzymeLabels(j)) = mutEffect(j)*signalWT.kcatij(enzymeLabels(j));
			signal.Kij(enzymeLabels(i)) = mutEffect2(i)*signalWT.Kij(enzymeLabels(i));
			signal.Kij(enzymeLabels(j)) = mutEffect2(j)*signalWT.Kij(enzymeLabels(j));

			n0 = runToSteadyState(n0, signal);

			epistat.doubleMut{i,j} = n0-n0_WT;

			fitval = exp(-0.5*transpose(n0-n0_WT)*epistat.matrixOfPain*(n0-n0_WT));;

			epimat(i,j) = fitval - (fitvalWT + epivec(i) + epivec(j));
			% disp(howfar({[i enzymeNumber],[j enzymeNumber]}));
		end
	end
	epistat.epivec = epivec;
	epistat.epimat = (epimat + transpose(epimat))/mean(abs(epivec));
	epistat.enzymeLabels = enzymeLabels;
	epistat.enzymeNumber = enzymeNumber;
	epistat.mutEffect = mutEffect;
end


function epistat = redoTheEpiEdges(signal,epistat, n0)
	% WT
	n0_WT = n0;
	enzymeNumber = epistat.enzymeNumber;

	fitval = exp(-0.5*transpose(n0_WT-n0_WT)*epistat.matrixOfPain*(n0_WT-n0_WT));
	fitvalWT = fitval;
	epivec = zeros(enzymeNumber,1);
	epimat = zeros(enzymeNumber);
	for i = 1:epistat.enzymeNumber
		n0 = n0_WT + epistat.singleMut{i};

		fitval = exp(-0.5*transpose(n0-n0_WT)*epistat.matrixOfPain*(n0-n0_WT));
		% keyboard

		epivec(i) = fitval - fitvalWT;
	end

	for i = 1:(enzymeNumber-1)
		for j = (i + 1):enzymeNumber
			n0 = n0_WT + epistat.doubleMut{i,j};

			fitval = exp(-0.5*transpose(n0-n0_WT)*epistat.matrixOfPain*(n0-n0_WT));;

			epimat(i,j) = fitval - (fitvalWT + epivec(i) + epivec(j));
		end
	end
	epistat.epivec = epivec;
	epistat.epimat = (epimat + transpose(epimat))/mean(abs(epivec));
	epistat.enzymeNumber = enzymeNumber;
end













%%%%% The Favs


function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',12);
	box on
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
