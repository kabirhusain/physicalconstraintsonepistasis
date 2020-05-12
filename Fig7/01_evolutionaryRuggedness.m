clear all;
close all;

epsvec = [0:.05:1];
repNum = 100;
landNum = 100;

% Generate a fitness landscape
sequenceL = 50;
epsilon = 0;
diagmat = diag(ones(sequenceL,1));
diagel = find(diagmat == 1);

matrixOfOnes = ones(size(diagmat));

matrixOfOnes(diagel) = 0;
uppermat = triu(matrixOfOnes);
uppervals = find(uppermat > 0);

for epsnum = 1:length(epsvec)
	epsilon = epsvec(epsnum);
	numberOfUnique{epsnum} = [];
	for land = 1:landNum
		hiOther = normrnd(0,1,sequenceL,1);
		thetai = normrnd(0,1,sequenceL,1);
		JijLow = thetai*thetai';
		JijFull = normrnd(0,1,sequenceL,sequenceL);
		JijFull(uppervals) = 0;
		JijFull = JijFull + JijFull';

		JijLow(diagel) = 0;
		JijFull(diagel) = 0;

		% do some adaptive walks.
		Jij = epsilon*JijLow + (1-epsilon)*JijFull;
		stepvec{epsnum} = [];
		spinmat = zeros(repNum, sequenceL);
		for rep = 1:repNum
			% spinvec = (randi(2, sequenceL, 1) - 1.5)*2;
			spinvec = randi(2, sequenceL, 1) - 1;
			fitness = hiOther'*spinvec + spinvec'*Jij*spinvec/2;
			notdone = 1;
			stepstaken = 0;
			envec = [fitness];
			while notdone
				maxfit = fitness;
				whichi = 0;
				for i = 1:sequenceL
					newspin = spinvec;
					newspin(i) = 1 - spinvec(i);
					% newspin(i) = -1*spinvec(i);
					newfitness = hiOther'*spinvec + newspin'*Jij*newspin/2;
					if newfitness > maxfit
						maxfit = newfitness;
						whichi = i;
					end
				end

				if whichi > 0
					stepstaken = stepstaken + 1;
					spinvec(whichi) = 1 - spinvec(whichi);
					% spinvec(whichi) = -1*spinvec(whichi);
					fitness = maxfit;
					envec = [envec fitness];
				else
					notdone = 0;
				end
			end
			spinmat(rep,:) = spinvec';
			stepvec{epsnum}(rep) = stepstaken;
		end
		numberOfUnique{epsnum}(land) = size(unique(spinmat,'rows'),1);
		disp([howfar({[ epsnum, length(epsvec) ], [land landNum]})])
	end
end


figure;
hold on;
maxval = 0;
for i = 1:length(epsvec)
	plot([epsvec(i) epsvec(i)], [mean(numberOfUnique{i}) - 1*std(numberOfUnique{i}) mean(numberOfUnique{i}) + std(numberOfUnique{i})], '-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5])
	plot(epsvec(i), mean(numberOfUnique{i}), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6)
	maxval = max(maxval,  mean(numberOfUnique{i}) + std(numberOfUnique{i}));
end
hold off
makePretty
axis([-.05 1.05 0 1.1*maxval]);
pbaspect([2.5 1 1])
saveas(gcf,'withRandomField.svg');
saveas(gcf,'withRandomField.fig');



function stringy = howfar(pairs)
	stringy = '';
	for i = 1:length(pairs)
		stringy = [stringy num2str(pairs{i}(1)) ' of ' num2str(pairs{i}(2)) ' '];
	end
end

function makePretty()
	set(gca, 'LineWidth',2)
	set(gca ,'FontSize', 12)
	box on
end