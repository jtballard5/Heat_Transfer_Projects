function [T] = fire5(t)
%13-nodes materials equal thickness

%Constants
time=0;
delta_x=0.004;
delta_t=0.1;
h=90;

%Material Properties
%Material 1
k4=10.1;
p4=8940;
Cp4=427;


%Material 2
k2=0.032;
Cp2=72;
p2=1680;

%Material 3
k3=0.023;
Cp3=24;
p3=1600;

%Material 4
%k4=15.6;
%p4=7913;
%Cp4=456;
k=0.062;
Cp=920;
p=2115;

%Air Gap 
kAir=15.6;
pAir=7913;
CpAir=456;



%Temperatures
T_ambient=800;
T_initial=27;

%Temperature Array
T=zeros(13,1);
T(:,:)=T_initial;
Tnew(:,:)=T(:,:);

%Begin loop for 5-node Explicit Finite Difference equations
while time<t
    
    
    Tnew(1,1)= (((2*delta_t)/(delta_x*p*Cp))*((h*(T_ambient-T(1,1)))+((k/delta_x)*(T(2,1)-T(1,1)))))+T(1,1);
    
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
     
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
    
    Tnew(4,1) = ((2*delta_t/((delta_x^2)*(p*Cp+p2*Cp2)))*((k*(T(3,1)-T(4,1)))+((k2)*(T(5,1)-T(4,1)))))+T(4,1);
    
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
      
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
        
    Tnew(7,1) = ((2*delta_t/((delta_x^2)*(p2*Cp2+p3*Cp3)))*((k2*(T(6,1)-T(7,1)))+(k3*(T(8,1)-T(7,1)))))+T(7,1);
    
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
     
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
    
    Tnew(10,1) = ((2*delta_t/((delta_x^2)*(p3*Cp3+p4*Cp4)))*((k3*(T(9,1)-T(10,1)))+((k4)*(T(11,1)-T(10,1)))))+T(10,1);
    
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
     
    Tnew(8,1)= (((k4*delta_t)/(p4*Cp4*delta_x^2))*(T(7,1)+T(9,1)-2*T(8,1)))+T(8,1);
     
    Tnew(13,1) = (((2*delta_t*k4)/((delta_x^2)*p4*Cp4))*(T(12,1)-T(13,1)))+T(13,1);
    
    
    T(:,:)=Tnew(:,:);
    time=time+delta_t;
end
end

