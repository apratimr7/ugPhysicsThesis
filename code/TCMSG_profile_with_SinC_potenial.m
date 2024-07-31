% This program is modified by Anjali Saini to determine the effect of potential on CS dynamics
% in 2D VCSEl incorporated with graphene saturable absorber and frequency selective feedback

clear all;
close all;

theta = 1.1; %detuning parameter
mu1 = 1.37; %pump parameter respectively for active  medium.   % range 1.39 - 1.507
alpha =2.7; %linewidth enhancement factor for active
gamma =0.55; %pump parameter respectively for passive medium.
beta = 1.4; %linewidth enhancement factor for passive
s = 10;
G = 4; %g1=g2=G=4
lambda = 0.6; %feed back bandwidth
sigma =0.4; %feedback strength
omega0 = 1.2;%resonance frequency
A=0.25;%amplitude of CS
alphaNS= 0.020;%amount of absorption which cannot be saturated
a1 = sigma*lambda^2/(lambda^2+omega0^2);
b1 = sigma*lambda*omega0/(lambda^2+omega0^2);
peako=A.^2;% peak intensity



%% Section to make transeverse space in X and Y dimesion %%%
largo=10;%10;%100;% % spatial pulse width
% fixing positioning window & maximum
fmax=40/(2*pi*largo); % maximum spectral frequency
X0 = largo*40; % spectral window
%initialization of space and frequency points
df=1/X0; % sampled spatial freq width
dx=1.0/(2.0*fmax); % sampled spectral width
nu=1;

while dx*(2^nu) < X0 % loop to get 2^n number of points to ease FFT
    nu=nu+1;
end

nx=2^nu; % no. of sample points in 2^n
dx=X0/nx; % adjusting the sampling space width
fmax=1.0/(2.0*dx); % adjusting the max frequency
dw=2*pi*df;

% Definition of spectral and frequency variables
n1=nx/2;
n2=n1-1;
x1=(-n1:n2)*dx;
x2=(-n1:n2)*dx;
[X,Y]= meshgrid(x1,x2);
size(X)
w1=(-n1:n2)*dw;
w2=(-n1:n2)*dw;
[Kx,Ky]= meshgrid(w1,w2);

%% %% calculation of L_D L_NL & time steps %%%%%%
Aireff=55; % effective core area in micrometre^2
wavelength=1330;  %*1e-9;    % emmission wavelength (in meters) of VCSEL (NRadwell Thesis)
coef_dt=1.e-2; % a fraction useful to set temporal steps
gamma1=2*pi*(2.7e-20)/(280e-9*Aireff*(1e-6)^2);
L_NL=1/(gamma1*peako);
L_D=4*pi^2*largo^2/wavelength;  %diffraction length in meter
% L_D=largo*largo/abs(gvd);5this is for fibre where dispersion will occur
Lmin=min([L_NL L_D]); % minimum length effect is more
% L=1*a*L_D; % total length of the fiber set to one dispersion length
L=wavelength;
tstep=coef_dt*Lmin; % step size in time direction
npt=round(L/tstep); % no. of steps in time direction
tstep=(L/npt); % re-adjusting the tstep for whole number of no. of steps

%% Input pulse form of field%%%
f=w1/(2*pi);
% u1=A*exp(-(X/largo).*(Y/largo));
% u1=A*sech((X)/largo).*sech((1.5*Y)/(largor));
% u1=A*sech((X)/largo).*sech((Y)/(largo)); %this is for original trapping

%u1=A*(exp(-(X).^2./largo^2-(Y).^2./largo^2));

delta = 5;
u1=A*sin(delta*(X.*Y./largo^2)).*(exp(-(X).^2./largo^2-(Y).^2./largo^2));

 %u1=A*(exp(-(X).^2./largo^2-(Y).^2./largo^2));
%u1 = A*(cosh(X)*sinh(Y));
UU = u1;%required for saving field
Vpt=fftshift(fft2(u1)); % Fourier transform
d = abs(u1); %input intensity
k=max(max(d)); %maximum intesity, same as peako
KK=k;
figure;
mesh(x1,x2,d);
ax=gca;
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
ax.ZAxis.FontSize=15;
ax.LineWidth=3;
ax.FontWeight = 'normal';
title( ["Input", num2str(k)] );
xlabel("x",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
ylabel("y",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
zlabel("|E|","Rotation",0,"FontSize",20,"FontWeight","bold")

figure
contour(x1,x2,d,'LineWidth', 3);
ax=gca;
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
ax.ZAxis.FontSize=15;
ax.LineWidth=3;
ax.FontWeight = 'normal';
xlabel("x",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
ylabel("y",'FontName','Gabriola','FontSize',20,"FontWeight","bold","Rotation",0)

%% Defining potential

% POT = 1*(0.0008).*(cos(0.031*X).*cos(0.031*Y));%trapped at centre without forming other CS
% POT = 1*(0.00085).*(cos(0.031*X).*cos(0.031*Y));%trapped at centre without forming other CS
%POT = 1*(0.000).*(cos(0.031*X).*cos(0.031*Y));%trapped at centre without forming other CS

r = sqrt(X.^2+Y.^2);
height_potential = 155;
period_potential = 20;
POT = height_potential.*(sinc(period_potential*r));

%delta = 15;
%POT = sin(delta*(X.*Y./largo^2)).*(exp(-(X).^2./largo^2-(Y).^2./largo^2));

%bessel_order = 2
%POT = height_potential.*besselj(2,period_potential*r)

OO= max(max(POT)); % peak amplitude of potential
figure;
mesh(X,Y,POT);
hold on;
SC=meshc(X,Y,POT);
%ax.ZLim(2)=0.001
SC(2).ZLocation="zmax";

ax=gca;
% ax.Zlim(2)=
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
ax.ZAxis.FontSize=15;
ax.LineWidth=3;
ax.FontWeight = 'normal';

xlabel("x",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
ylabel("y",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
zlabel("V(x,y)","Rotation",0,"FontSize",20,"FontWeight","bold")
drawnow;
figure;
contour(X,Y,POT,'LineWidth', 2);
ax=gca;
ax.XAxis.FontSize=15;
ax.YAxis.FontSize=15;
ax.ZAxis.FontSize=15;
ax.LineWidth=3;
ax.FontWeight = 'normal';
xlabel("x",'FontName','Gabriola','FontSize',20,"FontWeight","bold")
ylabel("y",'FontName','Gabriola','FontSize',20,"FontWeight","bold","Rotation",0)
drawnow;

%% Temporal evolution of input profile is sudied by itetration
t = 0; % initial time
li1 = 1; % running variable used for storing the data during propagation
li=1;%variable to draw trajectory
U = u1; % U for using in the SSFM loop

%%% loop for number of soliton periods %%%%%
iteration = 5; % number of soliton periods
niter=400; % value set to save the pulse field envelope
%%%%%%%%%%%%
f1 = figure;
%nm = char(floor(rand()*(123-97))+97);
nm = join([string(height_potential),string(period_potential)],'-')

%v = VideoWriter(plus("Sinc_",nm),"MPEG-4");

% v = VideoWriter(plus("Sinc_",nm),"MPEG-4");

v = VideoWriter(plus(nm,plus("TCMSG_",string(delta))),"MPEG-4");

open(v)

f2 = figure;

%v2 = VideoWriter(plus("Sinc_cont_",nm),"MPEG-4");

% v2 = VideoWriter(plus("Sinc_cont_",nm),"MPEG-4");

v2 = VideoWriter(plus(nm,plus("TCMSG_cont_",string(delta))),"MPEG-4");

open(v2)
for li2=1:iteration
    h=tstep;

    %%%% split-step Fourier transform method loop starts here %%%%
    diffraction=1*(-1+1i*theta+a1-1i*b1-1i*Kx.^2-1i*Ky.^2).*(h/2);
    for indice=1:npt
        % first-half diffraction
        V = fftshift(ifft2(U)); % Fourier transform
        V = V.*exp(diffraction); % calculation of diffraction
        U = (fft2(fftshift(V))); % Inverse Fourier transform
        %full nonlinearity
        P = abs(U).^2; % Intensity
        nonlinearity=1*(1i*POT+mu1*(1-1i*alpha)./(1+G*P)-gamma*(1-1i*beta)./(1+s*G*P)+alphaNS).*h;
        %  nonlinearity=1*(mu1*(1-1i*alpha)./sqrt((1+P))-gamma*(1-1i*beta)./sqrt((1+s*G*P))+alphaNS).*h;
        % nonlinearity=1*(mu1*(1-1i*alpha)./(1+G*P)-gamma*(1-1i*beta)./(1+s*P)+alphaNS).*h;                             % Active
        U = U.*exp(nonlinearity); % calculation of nonlinear
        %second-half diffraction
        V = fftshift(ifft2(U)); % Fourier transform
        V = V.*exp(diffraction); % calculation of diffraction
        U = (fft2(fftshift(V))); % Inverse Fourier transform
        % Incrementation of t
        t = t + h;
        %% condition to calculate the pulse parameters. if indice=niter then
        %%%%% pulse parameters are calculated for each tstep %%%%%%%%%%%%%

        if mod(indice,niter) == 0
            li1=li1+1;
            % saving the field
            UU =[UU; U];
            %%%% saving the temporal data %%%%%%%%%%%%%%%%%
            tt(li1)=t;
            D=abs((U));
            K=max(max(D));
            KK = [KK;K];


            % % % % plotting 3d plot of field after ssfm
            set(0, 'CurrentFigure', f1);
            set(gcf,'color','w');
            mesh(X,Y,D);
            ax=gca;
            ax.XAxis.FontSize=15;
            ax.YAxis.FontSize=15;
            ax.ZAxis.FontSize=15;
            ax.LineWidth=2;
            ax.FontWeight = 'normal';
            hold on
            % plot3(X, 200*ones(size(Y)), S); % project in x-z axis at y=200
            plot3(50*ones(size(X)), Y, max(D,[],2)); % project in y-z axis at x=200
            % plot3(X, Y,-0.1*ones(size(D))); % project in y-z axis at x=200

            title( [ "At time ",num2str(t)],FontSize=25 );
            xlabel("x",'FontName','Gabriola','FontSize',20,"FontWeight","bold",'Position',[59.558882937006274,-268.3584679325122,-0.031134065281436])
            ylabel("y",'FontName','Gabriola','FontSize',20,"FontWeight","bold",'Position',[-270.178730329013,49.55470237863301,-0.005902343738333])
            zlabel("|E|","Rotation",0,"FontSize",20,"FontWeight","bold",'Position',[-274.8251 213.229 0.25])
            %ax.XLim=[-50 50]; % comment or uncomment
            %ax.YLim=[-50 50];

            drawnow;
            video = getframe(gcf);
            writeVideo(v,video)
            hold off

            %plot *contour plots* of CS with contour plot of potenital
            set(0, 'CurrentFigure', f2)
            clf();
            set(gcf,'color','w');
            ax1 = axes;
            contour(ax1,X,Y,POT,'LineWidth', 1.5);
            title( [ 'At time',num2str(t)],FontSize=25 );
            xlabel(ax1,"x",'FontName','Gabriola','FontSize',20,"FontWeight","bold",'Position',[39.81781991419457 -230.7146732453936 1.000000000000014]);
            ylabel(ax1,"y",'FontName','Gabriola','FontSize',20,"FontWeight","bold",'Position',[-250.7480754539066 -46.88837655891389 1.000000000000014],Rotation=0.0);
            ax2 = axes;
            [CS,HI]=contour(ax2,X,Y,D,'LineWidth', 2);
            % linkaxes([ax1,ax2])
            hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];
            ax1.XAxis.FontSize=15;
            ax1.YAxis.FontSize=15;
            ax1.FontWeight = 'normal';
            ax1.LineWidth=2;

            %ax1.XLim=[-50 50]; % comment or uncomment
            %ax1.YLim=[-50 50];

            colormap(ax1,'parula')
            colormap(ax2,'jet')
            %%Then add colorbars and get everything lined up
            set([ax1,ax2],'Position',[0.211458333333333 0.1792656587473 0.643541666666666 0.745734341252686]);

            cb1 = colorbar(ax1,'Position',[0.08125,0.173416407061266,0.018750000000001,0.74558670820353],'FontSize',15,'FontWeight','bold');
            cb1.Label.String = 'V';
            cb1.Label.Position = [1.001851930220922,-0.2940983007688791,0]; % to change its position
            cb1.Label.Rotation = 0; % to rotate the text

            cb2 = colorbar(ax2,'Position',[0.913541666666665,0.176531671858775,0.018750000000001,0.74558670820353],'FontSize',15,'FontWeight','bold');
            cb2.Label.String = '|E|';
            pos2 = get(cb2,'Position');
            cb2.Label.Position = [1.09602271185981,0.041829903602393,0]; % to change its position
            cb2.Label.Rotation = 0; % to rotate the text
            %clim([0 max(max(D))]);
            drawnow;
            video2 = getframe(gcf);
            writeVideo(v2,video2)
            hold off

        end % if

    end % for indice
    %print("done")
    %close(vid2)
end % for iteration
close(v)

close(v2)


