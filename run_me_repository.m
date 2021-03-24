close ALL
numb_traj=250;
k_stp=80;
v_delay=9;
v_noise_multiplier=2.3;
[L,A,B,K,H,C,C_2,D0_p,D0_v,eta,nu,Q,R]=calculate_feedback_gains(k_stp,v_delay,v_noise_multiplier);
[x_all,x_hat_all,y_all,control_all] = generate_reach_trajectories(numb_traj,L,A,B,K,H,C,C_2,D0_p,D0_v,eta,nu,Q,R,v_delay);

%% plot cursor and hand paths
numb_cond=size(x_all,1);
plt_titles={'cursor position','hand position','hand velocity','applied force','external force'};
cols={'k','k--','b','r','c','r--'};
time=linspace(-200,k_stp*10-200-10,k_stp);
%hand paths
figure(1)
hold on 
for cond=1:numb_cond
    subplot(2,3,cond)
    hold on 
    xtmp=squeeze(x_all(cond,:,1,:));
    ytmp=squeeze(x_all(cond,:,8,:));
    for i=1:10:numb_traj
        plot(xtmp(i,1:2:end),ytmp(i,1:2:end),cols{cond})
    end
end

%cursor paths
figure(2)
hold on 
for cond=1:numb_cond
    subplot(2,3,cond)
    hold on 
    xtmp=squeeze(x_all(cond,:,2,:));
    ytmp=squeeze(x_all(cond,:,9,:));
    for i=1:10:numb_traj
        plot(xtmp(i,1:2:end),ytmp(i,1:2:end),cols{cond})
    end    
end

%plot managements
for f=1:2
   figure(f)
   for cond=1:numb_cond
      subplot(2,3,cond)
      xlim([-0.18,0.18])
      ylim([-0.06,0.3])
      axis square
      plot(0.02*cos(0:0.1:2*pi),0.02*sin(0:0.1:2*pi),'k')
      plot(0.02*cos(0:0.1:2*pi),0.02*sin(0:0.1:2*pi)+0.25,'k')
   end
end


%%
for state_variable=1:5
    figure
    hold on 
    for cond=1:numb_cond
        xtmp=squeeze(x_all(cond,:,state_variable,:));
        for i=1:numb_traj
            plot(xtmp(i,:),cols{cond})
        end    
    end
    title(plt_titles{state_variable})
    
end


for state_variable=1:5
    figure
    hold on 
    for cond=1:numb_cond
        xtmp=squeeze(x_all(cond,:,state_variable,:));
        shadedErrorBar(time,mean(xtmp,1),std(xtmp,[],1),cols{cond})
    end
    title(plt_titles{state_variable})
end


figure
hold on 
subplot(1,2,1)
hold on 
for cond=1:2
    xtmp=squeeze(x_all(cond,:,1,:));
    plot(time,std(xtmp,[],1),cols{cond})
end
xlabel('Time (ms)')
ylabel('Standard deviation position (m)')


subplot(1,2,2)
hold on 
for cond=3:4
    xtmp=squeeze(x_all(cond,:,1,:));
    plot(time,std(xtmp,[],1),cols{cond})
end
xlabel('Time (ms)')
ylabel('Standard deviation position (m)')

% plot bar plot for control output from 0-200ms post-perturbation
figure
hold on 
arr_control=[];
for cond=3:6
    if cond<6
        xtmp=squeeze(mean(x_all(cond,:,4,21+9:21+20),4));
        bar(cond,mean(xtmp),cols{cond}(1))
        errorbar(cond,mean(xtmp),std(xtmp),cols{cond}(1))
        arr_control(1,cond-2)=mean(xtmp);
    else
        xmech=squeeze(mean(x_all(4,:,4,21+9:21+20),4));
        xvis=squeeze(mean(x_all(5,:,4,21+9:21+20),4));
        bar(cond,mean(xmech)+mean(xvis),'k')
        arr_control(1,cond-2)=mean(xmech)+mean(xvis);
    end
end
xticks(3:6)
xticklabels({'Mech with vision','Mech without vision','Cursor Pert','Additive Model'})
xtickangle(45)
ylabel('Control Output (a.u.)')
title('Control from 90-200ms epoch')




function [x_all,x_hat_all,y_all,control_all] = generate_reach_trajectories(numb_traj,L,A,B,K,H,C,C_2,D0_p,D0_v,eta,nu,Q,R,v_delay)
    dim=size(A,1);
    states=int32(size(H,1)/2);
    k_stp=size(L,3)+1;
    H_p=H(1:states,:);
    H_v=H(states+1:end,:);

    %plots hand trajectories for 4 condition:
    % condition 1: unpert reach
    % condition 2: unpert reach without vision
    % condition 3: mech pert reach
    % condition 4: mech pert without vision
    % condition 5: cursor pert

    numb_cond=6;
    %numb_traj=50;   %number of trajectories per condition

    %matrices containing the state vector (x_all), estimated state vector
    %(x_hat_all), observations (y_all) and control signal (control_all)
    x_all=zeros(numb_cond,numb_traj,dim,k_stp);
    x_hat_all=zeros(numb_cond,numb_traj,dim,k_stp);
    y_all=zeros(numb_cond,numb_traj,states*2,k_stp);
    control_all=zeros(numb_cond,numb_traj,k_stp);

    for cond=1:numb_cond
    Cost_avg=0;
    template=squeeze(mean(x_all(3,:,:,:),2));   %get hand position from condition 2 during mechanical pert

    for j=1:numb_traj
        %for each trajectory create vectors to hold the state, estimated state
        %observations and control signal
        x=zeros(dim,k_stp);
        x_hat=zeros(dim,k_stp);
        y_vec=zeros(states*2,k_stp);
        control=zeros(1,k_stp);

        %set up target posiiton
        x(13,:)=.25;   %hand tweek here
        x_hat(13,:)=.25;%ditto

        Cost=0;

        for i=2:k_stp
            %generate control signal
            u=-L(:,:,i-1)*x_hat(:,i-1);
            control(1,i)=u(1);

            %calculate cost function
            Cost=Cost+(u.'*R*u);
            Cost=Cost+x(:,i-1).'*Q(:,:,i-1)*x(:,i-1);

            %forward propagate state equation
            x(:,i)=A*x(:,i-1)+B*u+eta.*randn(dim,1)+randn(dim,1).*(C*u)+randn(dim,1).*(C_2*u);

            %calculate observation
            y=H*x(:,i-1)+[D0_p.*randn(states,1);D0_v.*randn(states,1)];  

            %% is visual feedback provided?
            if cond==2 || cond == 4 
                if i>20
                    H_no_vfdbck=zeros(states,dim);
                    H=[H_p;H_no_vfdbck];
                else
                    H=[H_p;H_v];
                end
            else
                H=[H_p;H_v];
            end

    %         % is prop feedback provided?
    %         if cond==5
    %             if i>20
    %                 H_no_pfdbck=zeros(states,dim);
    %                 H=[H_no_pfdbck;H_v];
    %             else
    %                 H=[H_p;H_v];
    %             end
    %         else
    %             H=[H_p;H_v];
    %         end

            %% is an external perturbation applied?
            if i>20
                if cond==3 || cond == 4 || cond == 6
                    x(5,i)=-5; %size of external load
                else
                    x(5,i)=0;
                end
            end

            %% is there a cursor perturbation applied
            if i>20 && i<50
                if cond == 5
                     
                    x(1,i)=x(1,20)+template(1,i);   %update the cursor position

                    if i>20+v_delay
                        y([15,17])=[x(1,20);y(17)]+template([1,3],i-v_delay);   %update vision's estimate for cursor velocity
    %                     y([15,16,17,19])=[x(1,20);y([16,17,19])]+template([1,2,3,5],i-v_delay);   %update vision's estimate for cursor velocity
                    end  
                end

                if cond == 6
%                     template=squeeze(mean(x_all(3,:,:,:),2));   %get hand position from condition 2 during mechanical pert 
                    x(1,i)=x(1,20)-1*template(1,i);   %update the cursor position

                    if i>20+v_delay
                        y([15,17])=-1*template([1,3],i-v_delay);   %update vision's estimate for cursor velocity
                    end  
                end
            end

            %forward propagate estimated hand position
            x_hat(:,i)=A*x_hat(:,i-1)+B*u+K(:,:,i-1)*(y-H*x_hat(:,i-1))+nu.*randn(dim,1);

            %store observation matrix
            y_vec(:,i)=y;

        end
        %calculate cost function
        Cost=Cost+x_hat(:,k_stp).'*Q(:,:,k_stp)*x_hat(:,k_stp);

        Cost_avg=Cost_avg+Cost;
        x_all(cond,j,:,:)=x(:,:);
        x_hat_all(cond,j,:,:)=x_hat(:,:);
        y_all(cond,j,:,:)=y_vec;
        control_all(cond,j,:)=control;    
    end

    fprintf('Average Cost: %f \n',Cost_avg/numb_traj);



    end


end

function [L,A,B,K,H,C,C_2,D0_p,D0_v,eta,nu,Q,R]=calculate_feedback_gains(k_stp,v_delay,v_noise_multiplier)

%state vector x: element number
%                     x position elements 1-8
%                     y position elements 9-16
%                     1,8: cursor position
%                     2,9: hand position
%                     3,10: hand velocity
%                     4,11: Applied Force
%                     5,12: External Forces
%                     6,13: target hand position
%                     7,14: target hand velocity


dt=.010; %time step
G=.15;  % viscous constant
tau=.066;   %time constant of linear filter
%k_stp=80;  %number of time steps (i.e. time horizon)
states=14;  %number of elements in state vector (8 unique elements for X and Y)
%v_delay=9; %number of time steps for visual delay
dim=states*v_delay; %total state vector's dimension once augmented for delays








%state update matrix for state vector x
%A matrix
%X Direction
A=zeros(dim,dim);
A(1,1)=1;               
A(1,3)=dt;
A(2,2)=1;               
A(2,3)=dt;
A(3,3)=1-G*dt;
A(3,4)=dt;
A(3,5)=dt;
A(4,4)=1-dt/tau;
A(5,5)=1;
A(6,6)=1;
A(7,7)=1;

%Y Direction
A(8,8)=1;
A(8,10)=dt;
A(9,9)=1;
A(9,10)=dt;
A(10,10)=1-G*dt;
A(10,11)=dt;
A(10,12)=dt;
A(11,11)=1-dt/tau;
A(12,12)=1;
A(13,13)=1;
A(14,14)=1;

%accounting for the augmented state vectors
for i=(states+1):dim
    A(i,i-states)=1;
end



%how control signals update state vector
%B matrix
B=zeros(dim,2);
B(4,1)=dt/tau;
B(11,2)=dt/tau;

%H matrix; Observability Matrices
%Proprioceptive Observability Matrix
H_p=zeros(states,dim);
for i=1:states
    if i~=1 && i~=8 
        H_p(i,states*5+i)=1;         
    end
end


H_v=zeros(states,dim);
for i=1:states
    if i~=2 && i~=9 && i~=4 && i~=11 && i~=5 && i~=12 
        H_v(i,dim-states+i)=1;         
    end
end

H=[H_p;H_v];



%sensory noise from Crevecoeur et al 2016

%proprioception
D0_p=zeros(states,1);
for i=1:states
    D0_p(i,1)=1;
end
D0_p=D0_p*sqrt(10^(-5));
omega_omega_p=diag(D0_p).^2;


%vision

D0_v=zeros(states,1);
for i=1:states
    D0_v(i,1)=1;
end
D0_v=D0_v*sqrt(10^(-5))/sqrt(v_noise_multiplier);
omega_omega_v=(diag(D0_v).^2);

% defining the total sensory noise term
temp=cat(1,D0_p,D0_v);
omega_omega=diag(temp).^2;






%motor noise additive from Nashed et al 2012
eta=zeros(dim,1);
eta(4,1)=10^(-1);
eta(11,1)=10^(-1);
omega_eta=eta*eta.';





%motor noise multiplicative (signal dependent-noise); 
%from Nashed et al 2012
C=B*[1 0; 0 1]*0.15;
C_2=B*[0 1;-1 0]*0.05;



%state dependent-noise
D=zeros(states,states);



%additive noise to estimation/prediction
%from Crevecoeur et al 2016
nu=zeros(dim,1);
for i=1:states
    nu(i,1)=1;
end

nu=nu*sqrt(0.8)*10^(-3);
omega_nu=diag(nu).^2;








%cost matrix from Nashed et al 2012
Q=zeros(dim,dim,k_stp);
p=zeros(dim,1);
p(1,1)=1;   %cursor position
p(6,1)=-1;



v=zeros(dim,1);
v(3,1)=1;   %hand velocity
v(7,1)=-1;

Q(:,:,k_stp)=p*p.'+v*v.';



p=zeros(dim,1);
p(8,1)=1;   %cursor position
p(13,1)=-1;




v=zeros(dim,1);
v(10,1)=1;  %hand velocity
v(14,1)=-1;

Q(:,:,k_stp)=Q(:,:,k_stp)+p*p.'+v*v.';




%define R matrix for control cost
r=10^(-6);
R=eye(2,2)*r;%/k_stp;



%Kalman filter for proprioception
K_p=zeros(dim,states,k_stp-1);


%Kalman filter for vision
K_v=zeros(dim,states,k_stp-1);


K=[K_p, K_v];


%Calculations below follow from Todorov 2005
for j=1:30
    %feedback gain L
    L=zeros(2,dim,k_stp-1);
    S_x=Q(:,:,k_stp);
    S_e=0;
    

    for i=k_stp-1:-1:1
        L(:,:,i)=pinv(R+B.'*S_x*B+C.'*(S_x+S_e)*C+C_2.'*(S_x+S_e)*C_2)*B.'*S_x*A;
        S_x=Q(:,:,i)+A.'*S_x*(A-B*L(:,:,i));
        S_e=A.'*S_x*B*L(:,:,i)+(A-K(:,:,i)*H).'*S_e*(A-K(:,:,i)*H);

    end

    
    S_e=omega_eta;
    S_x_hat=omega_nu;
   
    S_e_x=0;
    
    for i=1:k_stp
      
        K(:,:,i)=A*(S_e*H.')*pinv(H*S_e*H.'+omega_omega);
                
        if i<k_stp
            S_e=omega_eta+omega_nu+(A-K(:,:,i)*H)*S_e*A.'+C*L(:,:,i)*S_x_hat*L(:,:,i).'*C.'+C_2*L(:,:,i)*S_x_hat*L(:,:,i).'*C_2.';
            S_x_hat=omega_nu+K(:,:,i)*H*S_e*A.'+(A-B*L(:,:,i))*S_x_hat*(A-B*L(:,:,i)).'+(A-B*L(:,:,i))*S_e_x*H.'*K(:,:,i).'+K(:,:,i)*H*S_e_x.'*(A-B*L(:,:,i)).';
            S_e_x=(A-B*L(:,:,i))*S_e_x*(A-K(:,:,i)*H).'-omega_nu;
        end
    end

end


end