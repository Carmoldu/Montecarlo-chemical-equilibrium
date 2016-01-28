%This solver computes the equilibrium of chemical reactions at different
%temperatures and pressures through the use of MonteCarlo Method. For this
%solver to work, it is necessary to include the HGS functions, which can be
%found at https://github.com/OpenLlop/HGS .


%Data
P=10; %bar
T=3000; %K

tolerance=0.01;

%Species
Spec={'H2','O2','H2O','H','O','OH'};
   %H2 O2 H2O H  0 OH
A=[ 2 0 2 1 0 1  %H
    0 2 1 0 1 1];%O

A1=[ 2 0   %H
    0 2];%O

A2=[2 1 0 1  %H
    1 0 1 1];%O


b=[5    %H
   1];  %O


%Starting point of the iteration
ndaux=[0
       0
       0
       0];


%A*n=b   A has range 2. With 6 variables, this means 4DoF (thus, 4 elements 
%must be defined in order to have a solution for the system). 
%These 4 DoF will be set in vector nd=[nH20, nH,nO,nOH] and we eill give them
%random initial values.
%Then we will compute nd-->n:
%      A                    n
%  A1 |    A2
%[2 0 | 2 1 0 1]     *     [nH2 ]       =B
%[0 2 | 1 0 1 1]           [nO2 ]
%     |-->variables       -------
%         fix for     nd<-|[nH2O]
%                         |[nH  ]
%                         |[nO  ]
%                         |[nOH ]

%A1*[nH2 ]+A2*nd=B ==>       |[nH2 ]=A1^(-1)*(B-A2*nd)      
%   [nO2 ]               ni<-|[nO2 ]

%ni=A1\b-A1\A2*ndaux;

%n=[ndaux
%  ni];

%Then we can compute Free Gibbs Energy
%[~,~,~,~,~,~,~,Gaux,~]=hgsprop(Spec,n,T,P);



%Iteration Parameters
here=1;
iti=0;
err=999;
iter=0;
ni=[-999
    -999];
Gaux=9999999;


figure

while err>tolerance
    while ni(1)<0 || ni(2)<0 || nd(1)<0 || nd(2)<0 || nd(3)<0 || nd(4)<0
        
        nd=ndaux+(-1+2*rand(4,1))./(here*rand(4,1));

        ni=A1\(b-A2*nd);
        
        iter=iter+1;
    end

    n=[nd
       ni];

    [~,~,~,~,~,~,~,G,~]=hgsprop(Spec,n,T,P);
    
    
    if G<Gaux
        ndaux=nd
        err=abs(Gaux-G)
        Gaux=G
        here=here+1;
        
    end
    
    ni=[-999
        -999];
    iti=iti+1;
    
    plot(iti,n)
    hold on
    drawnow
    %plot(iti,G)
    %hold on
end
