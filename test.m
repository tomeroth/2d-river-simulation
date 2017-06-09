clear;
clc;

%parameters
length=10;%[m]
width=4;%[m]
depth=1;%[m]
meanFlow=0.3;%[m/s]
dispersionCoeff=0.02;%[m^2/s]
injectionPoint=[2,2];%[m]

Cevap=0.1; % Coefficient of evaporation
River_flow=[5000 4500 4000 3000 2500 2000 2000 3000 5000 5500 5000 5000]; % input from source
V=500; % initial volume of water
Dem=[100 400 400 200 200 100 50 50 50 100 150 200]; % Water from rain
Rain=[0 0 0 0 0 0 50 300 500 0 0 0]; % Initial Demand
t=0; % STart time
DT=1; % Step size - 1-month
Tsim=120; % Simulation time (100 Yrs.)
i1=1;

%numerical parameters
dx = 0.1;% spacital resolution
dy = 0.2;
DT = 1;
%boundary conditions; 
%left side - Dirichlet
left=0 ;%c(0,t)=0
% right side - von Neumann
right=0;% 

%xNodes=100;
xNodes=round(length/dx)+1;
yNodes=round(length/dy)+1;
maxTime=1000;
maxTimeSteps=round(maxTime/DT);
river=zeros(xNodes,yNodes);

 %initial conditions
 for x=1:length
     for y=1:length
        if(mod(x*dx,10) == 0 && mod(y*dy,10) == 0)
            Demand=Dem(i1)*exp(0.003*t);
            Vin= Demand + River_flow(i1);
            Asurface=0.01*V;
            Evaporation=Asurface*Cevap;
            Seepage = 0.2*V;
            Tloss = Seepage + Evaporation;
            V = V+Vin-Tloss-Demand;
            river(x, y)= V;
            i1=i1+1;
            if i1>12
                i1=1;
            end;
        else
            Demand=Dem(3)*exp(0.003*t); %default
            river(x, y) = Demand;
         end;
         t=t+DT;
     end
end;
%loop
run=1;
n=1;
mesure=zeros(maxTimeSteps,1);

Ca=(meanFlow * DT) / (dx*dy);
Cd=(dispersionCoeff * DT) / (dx)^2/ dy^2;

w1=Cd*(1-Ca)-(Ca/6)*(Ca^2-3*Ca+2);
w2=Cd*(2-3*Ca)-(Ca/2)*(Ca^2-2*Ca-1);
w3=Cd*(1-3*Ca)-(Ca/2)*(Ca^2-Ca-2);
w4=Cd*Ca+(Ca/6)*(Ca^2-1);
xx=0:dx:length;
yy=0:dy:length;

endMass=0;
while run
    n=n+1;
    newRiver=zeros(xNodes, yNodes);
    for x=1:xNodes
        for y=1:yNodes
            %boundary condition left side
            if(x==1 || x==2 || y==1 || y == 2)
                newRiver(x,y)=left;
            else
                %boundary condition right side
                 if(x==xNodes)
                    newRiver(x,y)=right;
                 else                  
                   newRiver(x,y)=river(x,y)+ w1*river(x+1,y)-w2*river(x,y)+w3*river(x-1,y-1)+w4*river(x-2,y-2);            
                  end;
            end;      
        end
    end;
    river=newRiver;
    % stop simulation 
    if(n*DT>=maxTimeSteps) 
        run=0;
    end;
  
    endMass=sum(mesure) * DT;
    figure(1)
    [X,Y] = meshgrid(xx,yy);
    size(xx)
    size(yy)
    size(river)
    mesh(yy,xx,river);
    title(['Simulation time= ' num2str(n*DT) ' (s)']);
    xlabel('x(m) ');
    ylabel('y(m) ');
    zlabel('water volume');  

pause(0.01); 

end;
