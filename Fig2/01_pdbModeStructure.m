% Read the PDB, extract 3D locations of atoms. Find C-alpha, build a spring network, diagonalise Hessian
clear all;
close all;

cutoff = 10;
backboneSpring = 1;
aaSpringK = 1;

pdbnames = {'1i5s.pdb','1w7j.pdb' ,'1BFE.pdb', '1erk.pdb','105m.pdb', '2vkn.pdb', '1dkg.pdb'};
protName = {'Kif1A', 'MyosinV', 'PDZ3','Erk','Myoglobin', 'Sho1-SH3','Hsp70'};
chains = {'A', 'A' , 'A' , 'A' , 'A' , 'A' , 'A' ,'D'};
slowmodes = [1,2,2,1,1,2,2];

mkdir('ModeEnergies/')
mkdir('Projontolowest')

for nam = 1:length(pdbnames)
	pdbname = pdbnames{nam};
	pdb = pdbread(['pdbs/' pdbname]);

	alltheatoms = squeeze(struct2cell(pdb.Model.Atom));
	atomNames = alltheatoms(2,:);
	atomChains = alltheatoms(5,:);

	whichChain = strmatch(chains{nam},atomChains, 'exact');

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

	% springy.springMatrix = springy.springMatrix(1:2,1:2);
	hessian = dynamicalMatrix3D(springy,xyz);
	eigvals = sort(eig(hessian));
	eigvals = eigvals(7:end);

	[v,d] = eig(hessian);
	[bah, sortyeig] = sort(diag(d));
	lowestfellow = v(:,sortyeig(7:7+slowmodes(nam)-1));
	pinvhess = pinv(hessian);

	randEig = sort(rand(size(diag(d))));
	randEig(1:6) = 0;
	pinvrandHess = pinv(v*diag(randEig)*(v'));

	projontolowest = [];
	projontolowestRand = [];
	for rando = 1:5000
		randyvec = randn(length(lowestfellow),1);
		randyvec = randyvec/sqrt(sum(randyvec.^2));

		transofrandy = pinvhess*randyvec;
		transofrandy = transofrandy/sqrt(sum(transofrandy.^2));
		projontolowest = [projontolowest, sqrt(sum( (lowestfellow'*transofrandy).^2 ))];

		transofrandy = pinvrandHess*randyvec;
		transofrandy = transofrandy/sqrt(sum(transofrandy.^2));
		projontolowestRand = [projontolowestRand, sqrt(sum( (lowestfellow'*transofrandy).^2 ))];
	end

	figure('visible','off');
	hold on
	histogram(projontolowestRand,50, 'FaceColor',[0 0 0]/255,'Normalization','probability');
	histogram(projontolowest,50, 'FaceColor',[217 95 2]/255,'Normalization','probability');
	hold off
	pbaspect([0.6 1 1])
	axvals = axis;
	axis([-0.025 1.025 0 1.1*axvals(end)])
	makePretty
	set(gcf,'Renderer','painter')
	title(protName{nam})
	saveas(gcf,['Projontolowest/svg_' protName{nam} '_' pdbname(1:end - 4) '.svg'])
	saveas(gcf,['Projontolowest/' protName{nam} '_' pdbname(1:end - 4) '.png'])
	close

	f = plotTowerOfModes(hessian,3);
	% figure(f)
	axis([0.5 2.5 .5 2*max(eigvals)/min(eigvals)])
	set(gcf,'Renderer','painter')
	title(protName{nam})
	saveas(gcf,['ModeEnergies/svg_' protName{nam} '_' pdbname(1:end - 4) '.svg'])
	saveas(gcf,['ModeEnergies/' protName{nam} '_' pdbname(1:end - 4) '.png'])
	close

	disp(['Done with ' protName{nam}])

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

	% keyboard

	energy = sum(energymat(:))/2; % divide by 2 to get rid of overcounting
end


function makePretty()
	set(gca,'LineWidth',2)
	set(gca,'FontSize',12);
	box on
end

function f = plotTowerOfModes(jacob,whichcol)
	[v d] = eig(jacob);
	d = abs(real(diag(d)));
	modeenergies = sort(d);
	modeenergies = modeenergies(7:end);
	minen = modeenergies(1);

	f=figure('visible','off');
	hold on;
	% plot([1 2],[modeenergies(1) modeenergies(1)],'-','LineWidth',2,'Color',[217,95,2]/255)
	% for i = 1:whichcol
	% 	plot([1 2],[modeenergies(i) modeenergies(i)],'-','LineWidth',2,'Color',[217,95,2]/255)
	% end
	for i = 1:length(modeenergies)
		plot([1 2],[modeenergies(i) modeenergies(i)]/minen,'-','LineWidth',2,'Color',[0,0,0]/255)
		% 217,95,2
	end
	hold off;
	set(gca,'YScale','log')
	% axis([0.5 2.5 whichlim(1) whichlim(2)])
	pbaspect([0.4 1 1])
	box on
	% ax1 = gca;
	% ax1.XAxis.Visible = 'off'; 
	makePretty
end

function f = plotTowerOfModesWithIPRs(jacob)
	[v d] = eig(jacob);
	d = abs(real(diag(d)));
	[modeenergies sorty] = sort(d);
	modeenergies = modeenergies(7:end);

	v = v(:,sorty);
	v = v(:,7:end);

	f=figure('visible','off');
	hold on;
	% plot([1 2],[modeenergies(1) modeenergies(1)],'-','LineWidth',2,'Color',[217,95,2]/255)
	for i = 1:length(modeenergies)
		ipr = sum((v(:,i)).^4);
		plot([ipr],[modeenergies(i)],'ko','MarkerSize',7,'LineWidth',1,'MarkerFaceColor',[0,0,0]/255)
		% 217,95,2
	end
	hold off;
	set(gca,'YScale','log')
	set(gca,'XScale','log')
	% axis([0.5 2.5 whichlim(1) whichlim(2)])
	pbaspect([1 1 1])
	box on
	% ax1 = gca;
	% ax1.XAxis.Visible = 'off'; 
	makePretty
end

% Assume that equilibrium is unstrained
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