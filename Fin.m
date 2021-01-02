function [T] = Fin(t,h)
%2D Rectangular Fin Heat Transfer Simulation

time=0;
%Constants
k=200;
p=2700;
Cp=890;
delta_x=0.0025;
delta_y=0.0025;
x=0.05;
y=0.01;
T_ambient=25;
T_initial=250;
thermal_diffusivity = 0.000084;
delta_t=0.01;
Ac=y*1;
perimeter=2*1+2*0.01;
m=sqrt((h*perimeter)/(k*Ac));

%Construct Matrix
rows=1+y/delta_y;
columns=1+x/delta_x;
T=zeros(rows,columns);
T(:,:)=250;
Tnew(:,:)=T(:,:);
xdist=linspace(0,5,21);
ydist=linspace(0,1,5);
Ttheory=xdist;
Ttheory2=xdist;
Ttheory3=xdist;

while time<t
    time=time+delta_t;
    %Nodes
    %Top Row
    for column=2:1:20
        row=1;
        Tnew(row,column)=T(row,column)+(((2*delta_t)/(p*Cp*delta_x*delta_y))*((0.5*k*(T(row,column-1)-T(row,column)))+(0.5*k*(T(row,column+1)-T(row,column)))+(k*(T(row+1,column)-T(row,column)))+(h*delta_x*(T_ambient-T(row,column)))));
    
    
    end
    
    %Bottom Row
    for column=2:1:20
        row=5;
        Tnew(row,column)=T(row,column)+(((2*delta_t)/(p*Cp*delta_x*delta_y))*((0.5*k*(T(row,column-1)-T(row,column)))+(0.5*k*(T(row,column+1)-T(row,column)))+(k*(T(row-1,column)-T(row,column)))+(h*delta_x*(T_ambient-T(row,column)))));
    
    
    end
    
    %Interior
    for row=2:1:4
        for column=2:1:20
        
         Tnew(row,column)=T(row,column)+(((delta_t)/(p*Cp*delta_x*delta_y))*((k*(T(row,column-1)-T(row,column)))+(k*(T(row,column+1)-T(row,column)))+(k*(T(row+1,column)-T(row,column)))+(k*(T(row-1,column)-T(row,column)))));
    
    
        end
    end
    
    %Corner
     Tnew(1,21)=T(1,21)+(((4*delta_t)/(p*Cp*delta_x*delta_y))*((0.5*k*(T(1,20)-T(1,21)))+(0.5*k*(T(2,21)-T(1,21)))+(0.5*h*delta_x*(T_ambient-T(1,21)))+(0.5*h*delta_y*(T_ambient-T(1,21)))));
     Tnew(5,21)=T(5,21)+(((4*delta_t)/(p*Cp*delta_x*delta_y))*((0.5*k*(T(5,20)-T(5,21)))+(0.5*k*(T(4,21)-T(5,21)))+(0.5*h*delta_x*(T_ambient-T(5,21)))+(0.5*h*delta_y*(T_ambient-T(5,21)))));
     
     %Right Side
     for row=2:1:4
         column=21;
         Tnew(row,column)=T(row,column)+(((2*delta_t)/(p*Cp*delta_x*delta_y))*((0.5*k*(T(row+1,column)-T(row,column)))+(0.5*k*(T(row-1,column)-T(row,column)))+(k*(T(row,column-1)-T(row,column)))+(h*delta_y*(T_ambient-T(row,column)))));
         
     end
     
     if Tnew-T==0
        disp(time);
        break
      end
     
     T(:,:)=Tnew(:,:);
end

%Infinitely Long Fin
for i=1:21
Ttheory(i)=T_ambient + exp(-m*xdist(i)/100)*(250-T_ambient);
end
Qdot2=sqrt(h*perimeter*k*Ac)*(250-T_ambient)

%
%Adiabatic Fin Tip
for i=1:21
Ttheory2(i)=T_ambient + (cosh(m*(x-(xdist(i)/100)))/cosh(m*x))*(250-T_ambient);
end
Qdot3=sqrt(h*perimeter*k*Ac)*(250-T_ambient)*tanh(m*x)

%Convective Tip
for i=1:21
Ttheory3(i)=T_ambient + ((cosh(m*(x-(xdist(i)/100)))+(h/(m*k))*sinh(m*(x-(xdist(i)/100))))/(cosh(m*x)+(h/(m*k))*sinh(m*x)))*(250-T_ambient);
end
Qdot4=sqrt(h*perimeter*k*Ac)*(250-T_ambient)*(((sinh(m*x))+(h/(m*k))*cosh(m*x))/(cosh(m*x)+(h/(m*k))*sinh(m*x)))

plot(xdist,T(3,:),'g');
title('Centerline Temperature')
ylabel('Temperature (degrees C)')
xlabel('distance (cm)')

hold on
plot(xdist,Ttheory,'--');
plot(xdist,Ttheory2,'*--');
plot(xdist,Ttheory3,'b--o');
legend({'2D Explicit','Infinitely Long Fin','Adiabatic Fin Tip','Convective Tip'})
hold off

Q=0;

for row=2:1:4
    for column=2:1:20 
        
        Q=Q+(p*delta_x*delta_y*Cp)*(250-T(row,column));
                
    end
end

for column=2:1:20
    Q=Q + ((p*delta_x*delta_y*Cp)/2)*(250-T(1,column));
    Q=Q + ((p*delta_x*delta_y*Cp)/2)*(250-T(5,column));
end


for row=2:1:4 
    Q=Q + ((p*delta_x*delta_y*Cp)/2)*(250-T(row,21));
end
    
   
     Q=Q+((p*delta_x*delta_y*Cp)/4)*(250-T(1,21)); 
     Q=Q+((p*delta_x*delta_y*Cp)/4)*(250-T(5,21)); 
     
Q %displays total heat loss across the fin

Qdot=0; 

Qdot=Qdot + (h*delta_x/2*(T(1,1)-T_ambient)); 
Qdot=Qdot + (h*delta_x/2*(T(5,1)-T_ambient)); 
Qdot=Qdot + (k/2*(T(1,1)-T(1,2))); 
Qdot=Qdot + (k*(T(2,1)-T(2,2))); 
Qdot=Qdot + (k*(T(3,1)-T(3,2))); 
Qdot=Qdot + (k*(T(4,1)-T(4,2))); 
Qdot=Qdot + (k/2*(T(5,1)-T(5,2))); 
Qdot %total heat transfer rate through fin 

%contourf(xdist,ydist,T);
%title('2D Temperature Distribution')
%xlabel('distance (cm)')
%ylabel('distance (cm)')
%c = colorbar;
%c.Label.String = 'Temperature (degrees C)';
%plot(xdist,T(3,:));
%title('Centerline Temperature')
%ylabel('Temperature (degrees C)')
%xlabel('distance (cm)')