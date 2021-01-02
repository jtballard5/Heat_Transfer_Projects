function [T] = fire_implicit(t)
%6-nodes materials equal thickness

%Constants
time=0;
delta_x=0.018;
delta_t=0.01;
h=90;
index=(t/delta_t)+1;

%Material Properties
%Material 1
k=0.016;
p=100;
Cp=2.1;

%Material 2
k2=0.062;
Cp2=920;
p2=2115;

%Material 3
k3=0.032;
Cp3=72;
p3=1680;

%Material 4
k4=0.023;
Cp4=24;
p4=1600;

%Air Gap 
kAir=0.02588;
pAir=1.164;
CpAir=1007;


%Temperatures
T_ambient=800;
T_initial=27;

%Temperature Array
T=zeros(6,1);
T(:,:)=T_initial;
Tnew(:,:)=T(:,:);
nodeTemps=zeros(index,6);
nodeTemps(:,:)=T_initial;
i=2;r=1;

%Begin loop for 6-node Explicit Finite Difference equations
while time<t
   
    for j=1:100
    Tnew(1,1)= (( (2*delta_t)/(delta_x*p*Cp))*((h*(T_ambient-Tnew(1,1)))+((k/delta_x)*(Tnew(2,1)-Tnew(1,1)))))+T(1,1);
    
    Tnew(2,1) = ((2*delta_t/((delta_x^2)*(p*Cp+p2*Cp2)))*((k*(Tnew(1,1)-Tnew(2,1)))+(k2*(Tnew(3,1)-Tnew(2,1)))))+T(2,1);
   
    Tnew(3,1) = ((2*delta_t/((delta_x^2)*(p2*Cp2+pAir*CpAir)))*((k2*(Tnew(2,1)-Tnew(3,1)))+(kAir*(Tnew(4,1)-Tnew(3,1)))))+T(3,1);
    
    Tnew(4,1) = ((2*delta_t/((delta_x^2)*(pAir*CpAir+p3*Cp3)))*((kAir*(Tnew(3,1)-Tnew(4,1)))+(k3*(Tnew(5,1)-Tnew(4,1)))))+T(4,1);
   
    Tnew(5,1) = ((2*delta_t/((delta_x^2)*(p4*Cp4+p3*Cp3)))*((k3*(Tnew(4,1)-Tnew(5,1)))+(k4*(Tnew(6,1)-Tnew(5,1)))))+T(5,1);
       
    Tnew(6,1) = (((2*delta_t*k4)/((delta_x^2)*p4*Cp4))*(Tnew(5,1)-Tnew(6,1)))+T(6,1);
    
    end
   
    
    T(:,:)=Tnew(:,:);
    
    %Step Forward in Time
    time=time+delta_t;
    
    %Matrix for temperatures over time
    nodeTemps(i,:)=T(:,:);
    i=i+1;
    
end
    
    for j=1:9001
        temperatures(j,:)=nodeTemps(r,:);
        r=r+100;
    end


writematrix(temperatures, 'implicit_temperatures_6node.csv');

end

