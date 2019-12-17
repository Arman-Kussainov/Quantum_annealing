%close all;
clear all; clc;
N=50;JJ=+1;k=1;
d_ata=ones(N,N);

dT=0.02;MCSteps=400;
MCSteps=MCSteps*(N*N);

for T=2.3:dT:2.5
	% thermalization production stage
	for m=1:round(0.2*MCSteps)
		x=round(rand*(N-1))+1;y=round(rand*(N-1))+1;
		% create alternative configuration d_ata2
		% expand it to meet periodic boundary conditions
		d_ata2=[          0,  d_ata(N,:), 0;...
			d_ata(:,N),    d_ata,    d_ata(:,1);...
			          0,  d_ata(1,:), 0];
		% need to account the changes in matrix size
		x=x+1;y=y+1;
		s_um=(d_ata2(x-1,y)+d_ata2(x+1,y)+d_ata2(x,y-1)+d_ata2(x,y+1));
		dE=2*JJ*d_ata2(x,y)*s_um;
		w=exp(-dE/(k*T));

		if w>rand
			d_ata(x-1,y-1)=-d_ata(x-1,y-1);
			figure(1);
			imagesc(d_ata);axis square;title(strcat('dE = ',num2str(dE), ', T=',num2str(T) ));drawnow;colorbar;
		end

	end

	% production stage
	M=0;
	for m=1:MCSteps
		x=round(rand*(N-1))+1;y=round(rand*(N-1))+1;
		% create alternative configuration d_ata2
		% expand it to meet periodic boundary conditions
		d_ata2=[0,d_ata(N,:),0;...
			d_ata(:,N),d_ata,d_ata(:,1);...
			0,d_ata(1,:),0];
		% need to account the changes in matrix size
		x=x+1;y=y+1;
		s_um=(d_ata2(x-1,y)+d_ata2(x+1,y)+d_ata2(x,y-1)+d_ata2(x,y+1));
		dE=2*JJ*d_ata2(x,y)*s_um;
		w=exp(-dE/(k*T));

		if w>rand
			d_ata(x-1,y-1)=-d_ata(x-1,y-1);
			%figure(2);
			%imagesc(d_ata);axis square;title(strcat('dE = ',num2str(dE), ', T=',num2str(T) ));drawnow;colorbar;
		end
	figure(1);
	M=M+sum(sum(d_ata));
%	plot(T,sum(sum(d_ata)/(N*N));hold on;title(m);drawnow;grid on;
	end

figure(1);
plot(T, M/(N*N)/MCSteps,'b.');title(strcat('T=',num2str(T)));drawnow;hold on;grid on;
xlabel('T, temperature'); ylabel('<M>, magnetization per spin');
end





