% add description here
%Sanket Deshpande

clear all;
close all;
digits(32)

%% Import the data and convert the wavenumber to wavelength (in ascending order)

dataImport = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Data\7_9_21_TiO2SiO2_277_82_nm.0.dpt', 'FileType', 'text');
data = zeros(size(dataImport));

for i=1:size(dataImport,1)
    data(i,1)=1e7/dataImport(size(dataImport,1)-(i-1),1);
    data(i,2)=dataImport(size(dataImport,1)-(i-1),2);
end

%% Set up the simulation

%wavelength range. Use this if you're defining the refractive indices using Sellmeier relation
%lambda = linspace(data(1,1),data(size(data,1),1),size(data,1)/10);

%use this definition of wavelength if importing refractive indices of
%materials from CSV file
% Au_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\Au_n_k.csv', 'FileType', 'text');
% Ti_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\Ti_n_k.csv', 'FileType', 'text');
% lambda(1,:) = Au_ref(:,1).*1000;

TiO2_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\TiO2_n.csv', 'FileType', 'text');
SiO2_ref = readmatrix('G:\My Drive\Projects\Quantum Array Generator\Thin Film matlab codes\TF theory vs data code\SiO2_n.csv', 'FileType', 'text');
lambda(1,:) = TiO2_ref(:,1).*1000;

%% free-space refractive indices of various media

%n_Au = 0.15659+4.9908i; %gold refractive index at 810nm
% n_Au = zeros(length(lambda),1);
% for j=1:size(Au_ref,1)
%     n_Au(j,1) = Au_ref(j,2)+Au_ref(j,3)*1i;
% end

%n_Ti = 3.1745+4.01i; %titanium refractive index at 810 nm
% n_Ti = zeros(length(lambda),1);
% for j=1:size(Ti_ref,1)
%     n_Ti(j,1) = Ti_ref(j,2)+Ti_ref(j,3)*1i;
% end

%n_FSi = 1.4531; %fused silica refractive index at 810 nm
% n_FSi = sqrt(1 + 0.6961663.*((lambda./1000).^2)./((lambda./1000).^2 - 0.0684043^2) + 0.4079426.*((lambda./1000).^2)./((lambda./1000).^2 - 0.1162414^2) + 0.8974794.*((lambda./1000).^2)./((lambda./1000).^2 - 9.896161^2) ); %Sellmeier relation from https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson
n_FSi = zeros(length(lambda),1);
for j=1:size(SiO2_ref,1)
    n_FSi(j,1) = SiO2_ref(j,2);
end

n_air = ones(length(lambda),1); %air refractive index

%n_TiO2 = 2.5173; %Titanium dioxide at 810 nm
%n_TiO2 = vpa(sqrt(5.913 + 0.2441./(lambda.^2 - 0.0803))); %Sellmeier relation from https://refractiveindex.info/?shelf=main&book=TiO2&page=Devore-o

n_TiO2 = zeros(length(lambda),1);
for j=1:size(TiO2_ref,1)
    n_TiO2(j,1) = TiO2_ref(j,2);
end

%% variable declarations
n_mat = vpa(zeros(4,length(lambda))); %matrix holding refractive indices of all layers
l_mat = zeros(4,1); %matrix holding lengths of all layers. Length units are nanometers
t_mat = zeros(3,length(lambda)); %transmission coeff matrix for all interfaces
r_mat = zeros(3,length(lambda)); %reflection coeff matrix for all interfaces
phi_mat = zeros(4,length(lambda)); %phase shift matrix for all layers
T_mat = zeros(2,2,3,length(lambda)); %Transfer Matrix for interfaces
P_mat = zeros(2,2,4,length(lambda)); %Phase shift Matrix for layers
M = zeros(2,2,length(lambda)); %Total transfer matrix
for i=1:length(lambda)
    M(:,:,i) = [1 0;0 1];
end

%% defining layers in the structure
%if you intend to use less than 4 layers, you can define layers 2 and 3 to
%be made up of the same material. For example, 20 nm of Au can be defined
%as 10 nm Au in layer 1 and 10 nm Au in layer 2

%layer 1 is air with no length
n_mat(1,:) = n_air;
l_mat(1) = 0;

%layer 2 (thin film)
n_mat(2,:)=n_TiO2;
l_mat(2)=245;

%layer 3 (thin film)
n_mat(3,:)=n_FSi;
l_mat(3)=40;

%layer 4 (substrate)
n_mat(4,:)=n_FSi;
l_mat(4)=42;

%calculating the transfer matrix
for i=1:length(lambda)
    for a=1:3
        t_mat(a,i)= 2*n_mat(a,i)/(n_mat(a,i)+n_mat(a+1,i));
        r_mat(a,i)= (n_mat(a,i)-n_mat(a+1,i))/(n_mat(a,i)+n_mat(a+1,i));

        T_mat(:,:,a,i) = t_mat(a,i)^(-1) * [1 r_mat(a,i); r_mat(a,i) 1];
    end
end

%calculating the propagation matrix
for i=1:length(lambda)
    for b=1:3
        phi_mat(b,i) = (2*pi/lambda(i))*n_mat(b+1,i)*l_mat(b+1);
        P_mat(:,:,b,i) = [exp(-1i*phi_mat(b,i)) 0; 0 exp(1i*phi_mat(b,i))];
    end
end

%calculating the total transfer matrix
for i=1:length(lambda)
    for c=1:3
        M(:,:,i) = M(:,:,i) * T_mat(:,:,c,i) * P_mat(:,:,c,i);
    end
end

phase=zeros(length(lambda),1);
T=zeros(length(lambda),1);

for i=1:length(lambda)
    %phase shift due to all layers in degrees
    phase(i) = angle(1/M(1,1,i)).*180/pi;

    %Transmittance from all the layers
    T(i)=(n_mat(size(n_mat,1),i)/n_mat(1,i))*(abs(1/M(1,1,i)))^2;
end

%plot the transmission vs lambda curve for TiO2 on SiO2
figure(1)
plot(lambda,T.*100,'-','Color','b','Linewidth',3)
hold on
plot(data(:,1),data(:,2).*100,'*','Color','g','Linewidth',3)
ylabel('Transmittance [%]')
yyaxis right
plot(lambda,phase,'-','Color','r','Linewidth',3)
xline(810,'-',{'810 nm'},'LabelVerticalAlignment','bottom','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
set(gca, 'XColor','k', 'YColor','k')
ylabel('Phase shift [deg]', 'Color','r')
hold off
title('TiO2 (xxnm) on SiO2 (82nm)')
xlabel('Free Space Wavelength [nm]')

legend('T (theoretical)','T (experimental)','\phi (theoretical)')
ax = gca;
ax.FontSize = 45;

% %plot the transmission vs lambda curve for Au-Ti
% figure(1)
% plot(lambda,T.*100,'-','Color','b','Linewidth',3)
% hold on
% plot(data(:,1),data(:,2).*100,'*','Color','g','Linewidth',3)
% ylabel('Transmittance [%]')
% yyaxis right
% plot(lambda,phase,'-','Color','r','Linewidth',3)
% xline(810,'-',{'810 nm'},'LabelVerticalAlignment','bottom','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
% set(gca, 'XColor','k', 'YColor','k')
% ylabel('Phase shift [deg]', 'Color','r')
% hold off
% title('Au (25nm) on Ti (8nm) on SiO2 (42nm)')
% xlabel('Free Space Wavelength [nm]')
% legend('T (theoretical)','T (experimental)','\phi (theoretical)')
% ax = gca;
% ax.FontSize = 45;

%plot the refractive index of TiO2 and SiO2 vs lambda curve
figure(2)
plot(lambda,n_TiO2(:,1),'+','Color','b','Linewidth',2,'MarkerSize',20)
hold on
plot(lambda,n_FSi(:,1),'+','Color','r','Linewidth',2,'MarkerSize',20)
ylabel('Magnitude')
xline(810,'-',{'810 nm'},'LabelVerticalAlignment','middle','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
set(gca, 'XColor','k', 'YColor','k')
hold off
title('Refractive index of TiO2 and SiO2')
xlabel('Free Space Wavelength [nm]') 
legend('n: TiO_2','n: SiO_2')
ax = gca;
ax.FontSize = 45;


% %plot the Au and Ti n,k vs lambda curve
% figure(3)
% plot(lambda,Au_ref(:,2),'+','Color','b','Linewidth',2,'MarkerSize',20)
% hold on
% plot(lambda,Au_ref(:,3),'*','Color','b','Linewidth',2,'MarkerSize',20)
% plot(lambda,Ti_ref(:,2),'+','Color','r','Linewidth',2,'MarkerSize',20)
% plot(lambda,Ti_ref(:,3),'*','Color','r','Linewidth',2,'MarkerSize',20)
% ylabel('Magnitude')
% xline(810,'-',{'810 nm'},'LabelVerticalAlignment','top','LabelOrientation','horizontal','Linewidth',2,'FontSize',25);
% set(gca, 'XColor','k', 'YColor','k')
% hold off
% title('Refractive index of Au and Ti')
% xlabel('Free Space Wavelength [nm]') 
% 
% legend('n: Au','k: Au','n: Ti','k: Ti')
% ax = gca;
% ax.FontSize = 45;

