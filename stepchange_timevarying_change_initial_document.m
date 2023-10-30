
clear;

clc;

kappa=10;

N=16;

epsilon=2; 

varepsilon=50; 

constraint_gain=200;

localcons_gain=10;

box_gain=100;

%=========================================Graph parameters====================================

Line1=[1;1;2;2;3;3;4;4; ...
       5;5;6;6;7;7;8;8; ...
       9;9;10;10;11;11;12;12; ...
       13;13;14;14;15;15;16;16; ...
       3;3;8;8;9;9;14;14];

Line2=[2;4;1;3;2;4;1;3; ...
       6;8;5;7;6;8;5;7; ...
       10;12;9;11;10;12;9;11; ...
       14;16;13;15;14;16;13;15; ...
       8;14;3;9;8;14;3;9];

BN=size(Line1,1)/2;

GN=size(Line1,1); 

Adjacent=zeros(N,N);                                                       %adjacent matrix 

for i=1:GN
     
         Adjacent(Line1(i),Line2(i))=1;                                    %adjacent matrix 
         
         Adjacent(Line2(i),Line1(i))=1;                                    %adjacent matrix 
              
end    
    
Degree=diag(Adjacent*ones(N,1));                                           %Degree matrix

Lap=1000*(Degree-Adjacent);                                                 %Laplacian matrix

Lap_BN_plus_N=kron(Lap,eye(BN+N+1));                                       %Laplacian matrix kronect

Lap_BN=kron(Lap,eye(BN));                                                  %Laplacian matrix kronect

%-------------------------------------Construct incidence matrix------------------------------ 

Lineindex1=zeros(BN,2);                                                    %initial matrix

Lineindex2=zeros(BN,3);

kk=1;                                                                      %initial count number

for i=1:GN                                                                 %set line belongs index, i.e. \mathcal{C}_i

    if isempty(intersect([min(Line1(i),Line2(i)),max(Line1(i),Line2(i))],Lineindex1,'row'))

        Lineindex1(kk,:)=[min(Line1(i),Line2(i)),max(Line1(i),Line2(i))];     
        
        Lineindex2(kk,:)=[Lineindex1(kk,:),kk];                            %include line number

        kk=kk+1;

    end

end

Lineindex=zeros(BN,4);

kk=1;

for j=1:N                                                                 %initial count number
    
    for i=1:BN
        
        if Lineindex1(i,1)==j

            Lineindex(kk,:)=[j,Lineindex2(i,:)];                            %Lines are managed by agent j
            
            kk=kk+1;
            
        end

    end

end

NCLN=zeros(N,1);

for i=1:N                                                                  %How many lines are under control of each agent
    
    if find(Lineindex(:,1)==i)~=0
    
        NCLN(i,1)=size(find(Lineindex(:,1)==i),1);                               %node control's lines number
        
    else
        
         NCLN(i,1)=0;                                                      %No line is under control of this node
        
    end
    
end


B=zeros(N,BN);                                                             %initial incidence matrix

for i=1:BN

    B(Lineindex1(i,1),i)=-1;                                                %Set the initial node as 1   
    
    B(Lineindex1(i,2),i)=1;                                                %Set the end node as -1 
    
end

rho_V_g_d=300*ones(N,1);

rho_I_l_d=100*ones(BN,1);


%========================================initial condition====================================

% x_ini=[42.9885180803861	53.1083898608177	39.5499847011037	36.9105442723708 ...
%        376.975097530925	377.480024559161	377.062747253705	377.117393398486 ...
%        -7.21324326051656	8.34554610911788	-0.683076809760794	-2.37159779268566 ...
%        42.9885180803861	53.1083898608177	39.5499847011037	36.9105442723708 ...
%        376.975097530925	377.480024559161	377.062747253705	377.117393398486 ...
%        -7.21324326051656	8.34554610911788	-0.683076809760794	-2.37159779268566 ...
%        42.9885180803861	53.1083898608177	39.5499847011037	36.9105442723708 ...
%        376.975097530925	377.480024559161	377.062747253705	377.117393398486 ...
%        -7.21324326051656	8.34554610911788	-0.683076809760794	-2.37159779268566 ...
%        42.9885180803861	53.1083898608177	39.5499847011037	36.9105442723708 ...
%        376.975097530925	377.480024559161	377.062747253705	377.117393398486 ...
%        -7.21324326051656	8.34554610911788	-0.683076809760794	-2.37159779268566 ]';

x_ini=zeros(2*N+BN,1);

sigma_ini=zeros(N,1);

u_ini=[378.9289  380.4146  379.7396  379.2243  380.1955  381.1042 ...
       379.9247  379.8974  380.0252  380.4045  379.6763  379.8341 ...
       380.1194  380.5682  379.3881  379.9696]';                                                     %set u_g_d initial value as 0        

x_hat_ini=[14.515064754915	3.95935897697232	32.9787419273847	14.9071228299889 ...
           37.6257931691129	61.1095883901259	41.4721611959256	33.6364102939300 ...
           44.9867290171943	52.8132621486412	28.1355733315040	33.0707668633162 ...
           36.6024320622732	49.3936276267038	2.38033527987224	36.4476162672594 ...
           378.638597580330	380.343321138855	379.211948658584	379.000752058610 ...
           379.367771584999	379.943120819626	378.846459424323	378.955535765946 ...
           379.260398628211	379.401077323141	379.001050112745	379.073436746389 ...
           379.167740114933	379.530912000852	379.340485178409	379.277056940884 ...
           32.7896634477947	63.8129104269148	14.1421560040464	-3.93339000437345 ...
           -1.60186543288746	5.15294773851377	12.9892279062272	-1.81793902663058 ...
           -1.75848368670712	2.67088402588984	6.66712017328114	-1.03409476637000 ...
           -4.53964857345814	-1.56166894232836	32.3649477996494	1.26856475033374 ...
           3.66304132364389	-5.31605570356793	-4.35518374697698	-3.38141715740433]';

fuzhu_ini=x_ini(1:N);

upsilon_ini=[172.552238898892	172.552238898892	172.552238898892	172.552238898892 ...
             172.552238898892	172.552238898892	172.552238898892	172.552238898892 ...
             172.552238898892	172.552238898892	172.552238898892	172.552238898892 ...
             172.552238898892	172.552238898892	172.552238898892	172.552238898892]';                                                    %Initial fast system states

nu_ini=[-0.00758212493845042	-0.0182492317603708	0.00714905390846842	0.0186908981861182 ...
        -0.00758212493845042	-0.0182492317603708	0.00714905390846842	0.0186908981861182 ...
        -0.00758212493845042	-0.0182492317603708	0.00714905390846842	0.0186908981861182 ...
        -0.00758212493845042	-0.0182492317603708	0.00714905390846842	0.0186908981861182]';                                                          %Initial fast system states

lambda_ini=zeros(N*(BN+N+1),1);

theta_ini=zeros(N*(BN+N+1),1);

lambda_v_low_ini=zeros(N,1);

lambda_v_up_ini=zeros(N,1);

lambda_I_low_ini=zeros(BN,1);    % is I_l

lambda_I_up_ini=zeros(BN,1);    % is I_l



%=========================================Power network parameters====================================

Dw=0.01;

alpha_w=0.8;


R=diag([2.0 1.8 1.6 1.5 ...
        2.2 1.9 2.6 2.8 ...
        1.7 1.9 2.4 2.3 ...
        2.6 2.1 2.0 1.9]*10^-2);

L=diag([1.8 2.0 3.0 2.2 ...
        2.2 2.3 1.9 2.5 ...
        1.7 2.2 2.3 2.9 ...
        2.3 2.5 1.9 3.0]*10^-3);

C=diag([2.2 1.9 2.5 1.7 ...
        1.8 2.2 2.9 3.1 ...
        1.7 2.5 2.6 1.9 ...
        1.9 1.9 2.2 2.3]*10^-3);

L_l=diag([2.1 2.0 3.0 2.2 ...
          2.3 2.6 2.0 1.9 ...
          2.3 2.5 2.8 2.4 ...
          2.3 2.5 2.6 2.4 ...
          2.1 2.3 2.5 2.6]*10^-6);

R_l=diag([7 5 8 6 ...
          8 8 5 6 ...
          8 7 6 7 ...
          8 7 8 5 ...
          7 6 7 8]*10^-2);

R_l_inverse=R_l^-1;

changetime=5;

LoadIchangeno=3;

LoadZchangeno=3;

changerate=3;

changevariance=0.1;

LoadI=[30.0;15.0;30.0;26.0; ...
       18.9;22.8;29.6;32.1; ...
       22.6;27.6;16.5;19.9; ...
       26.5;19.4;22.7;24.9];

LoadIchange=[30.0;15.0;30.0;26.0; ...
       18.9;22.8;29.6;32.1; ...
       22.6;27.6;16.5;19.9; ...
       26.5;19.4;22.7;24.9]+3*ones(N,1);

LoadZ=diag([16.7;50.0;16.7;20.0; ...
            25.0;16.0;14.2;29.4; ...
            18.5;19.6;22.6;25.6; ...
            23.4;26.4;25.2;33.7]);

LoadZchange=diag([16.7;50.0;16.7;20.0; ...
            25.0;16.0;14.2;29.4; ...
            18.5;19.6;22.6;25.6; ...
            23.4;26.4;25.2;33.7])+3*diag(ones(N,1));


%LoadI=zeros(4,1);
%LoadI=ones(4,1);

Lineconstrants=10*ones(BN,1);                                              %Transmission line limits

Line_loss=150;                                                               %Transmission line power loss limits

V_r=[380.16  381.31  380.27  380.13 ...
     380.12  381.24  380.15  380.35 ...
     380.15  381.13  380.09  380.22 ...
     380.19  381.26  380.21  380.19]';                             %Reference voltage 1    

V_hat_low=(380-3)*ones(N,1);                                                 %Voltage lower bound  

V_hat_up=(380+3)*ones(N,1);                                                  %Voltage upper bound

pr=0.005;                                                                  %Curnot model parameter

l=2.5;                                                                       %Curnot model parameter

r=[1.06 1.03 1.05 1.04 ...
   1.02 1.03 1.04 1.01 ...
   1.03 1.04 1.06 1.07 ...
   1.01 1.03 1.04 1.02];                                  %initial r_i

A_r=kron(diag(r),eye(N+BN+1));                                              %initial kenonect r

alpha_I=[10.6   10.1   10.2   10.5 ...
         10.7   10.9   10.1   10.8 ...
         10.6   10.5    9.6   11.3 ...
         10.3   10.9   10.2   10.5]'; 
 
alpha_V=[0.7    0.6    0.8    0.9 ...
         0.7    1.2    0.8    1.24 ...
         0.3    0.3    0.0    1.70 ...
         1.0    0.5    1.1    1.0]';                            %Optimization model parameter  xxxxxxxxxxxxxxxxxx

alpha_u=[1.1    1.4    1.8    1.1 ...
         1.1    1.9    1.9    1.2 ...
         1.3    1.4    1.8    0.6 ...
         1.8    2.1    1.5    0.9]';                            %Optimization model parameter  xxxxxxxxxxxxx

alpha_I_l=[1.4    1.9    1.7    1.9 ...
           0.5    1.2    1.6    1.6 ...
           0.7    1.3    1.8    1.4 ...
           0.3    1.4    1.7    1.4 ...
           1.1    1.9    2.0    1.8]';                         %Optimization model parameter  xxxxxxxxxx

%=========================================  Microgrid matrix  ===================================

DGU_inverse=blkdiag(L^-1,C^-1,L_l^-1);

M=[-R,-eye(N),zeros(N,BN); eye(N),-LoadZ^-1,B;zeros(BN,N),-B',-R_l];

A=[eye(N),     -LoadZ^-1,   B; zeros(BN,N),  R_l^-1*B',   eye(BN); zeros(1,N),  zeros(1,N),   V_r'*B ];

Mchange=[-R,-eye(N),zeros(N,BN); eye(N),-LoadZchange^-1,B;zeros(BN,N),-B',-R_l];

Achange=[eye(N),     -LoadZchange^-1,   B; zeros(BN,N),  R_l^-1*B',   eye(BN); zeros(1,N),  zeros(1,N),   V_r'*B ];

s_A=zeros(N*(N+BN+1),1);

for i=1:N

    s_A((i-1)*(N+BN+1)+1:i*(N+BN+1))=[zeros(N,1);zeros(BN,1);Line_loss/N];

    s_A((i-1)*(N+BN+1)+i)=LoadI(i);

end

s_Achange=zeros(N*(N+BN+1),1);

for i=1:N

    s_Achange((i-1)*(N+BN+1)+1:i*(N+BN+1))=[zeros(N,1);zeros(BN,1);Line_loss/N];

    s_Achange((i-1)*(N+BN+1)+i)=LoadIchange(i);

end

%%=========================================  Loacal constraints  ===================================

D=[R,eye(N),zeros(N,BN)];

F_x_I=zeros(BN,N);  %xxxxxxxxxxxxxxxxxxxxxxxx

        for i=1:N 

            if NCLN(i,1)~=0

                choose_line=find(Lineindex(:,1)==i);

                for j=1:NCLN(i,1)                                          %The numbers of the controlled line of agent i

                    aa=Lineindex(choose_line(j),4);                                  %choose which line
        
                    F_x_I(aa,aa)=r(i).*alpha_I_l(aa);    
        
                end

            end
            
        end

F_x=blkdiag(diag(r)*diag(alpha_I)+diag(r)*pr*diag(V_r),diag(r)*diag(alpha_V),F_x_I);

%---------------------------------ojection local constraint setting--------------------------------

normal_vector_nonstandar=zeros(3*N,1);

normal_vector=zeros(3*N,1);

reference_vector_low_1=zeros(3*N,1);

reference_vector_low_2=zeros(3*N,1);

bound_vector_low_nonstandar=zeros(3*N,1);


ound_low_vector=zeros(3*N,1);

reference_vector_up_1=zeros(3*N,1);

reference_vector_up_2=zeros(3*N,1);

bound_vector_up_nonstandar=zeros(3*N,1);

bound_up_vector=zeros(3*N,1);


for i=1:N                                                                  % for all nodes calculate the normal vector
    
    normal_vector_nonstandar((i-1)*3+1:i*3)=[-1/R(i,i);1;1/R(i,i)];                    % caulculate the Projection normal vector
    
    normal_vector((i-1)*3+1:i*3)=(1/norm(normal_vector_nonstandar((i-1)*3+1:i*3),2))...            % normalized the Projection normal vector
                     *[-1/R(i,i);1;1/R(i,i)];
    
end    

for i=1:N                                                                               % for all nodes calculate the bound vector
    
    reference_vector_low_1((i-1)*3+1:i*3)=[V_hat_low(i);0;V_hat_low(i)];                            % first reference node located at the lower bound
    
    reference_vector_low_2((i-1)*3+1:i*3)=[0;-V_hat_low(i)/R(i,i);V_hat_low(i)];                    % second reference node located at the lower bound
    
    bound_vector_low_nonstandar((i-1)*3+1:i*3)=reference_vector_low_1((i-1)*3+1:i*3)-reference_vector_low_2((i-1)*3+1:i*3); % vector along the lower bound
    
    bound_low_vector((i-1)*3+1:i*3)=(1/norm(bound_vector_low_nonstandar((i-1)*3+1:i*3),2))...                   % normalized the vector along the lower bound
                         *bound_vector_low_nonstandar((i-1)*3+1:i*3);
                    
    reference_vector_up_1((i-1)*3+1:i*3)=[V_hat_up(i);0;V_hat_up(i)];                               % first reference node located at the upper bound
    
    reference_vector_up_2((i-1)*3+1:i*3)=[0;-V_hat_up(i)/R(i,i);V_hat_up(i)];                       % second reference node located at the upper bound
    
    bound_vector_up_nonstandar((i-1)*3+1:i*3)=reference_vector_up_1((i-1)*3+1:i*3)-reference_vector_up_2((i-1)*3+1:i*3);    % vector along the upper bound
    
    bound_up_vector((i-1)*3+1:i*3)=(1/norm(bound_vector_up_nonstandar((i-1)*3+1:i*3),2))...                     % normalized the vector along the upper bound
                       *bound_vector_up_nonstandar((i-1)*3+1:i*3);
end  
