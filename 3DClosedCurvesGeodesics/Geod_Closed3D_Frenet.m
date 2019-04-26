function [dist,geod,frames_geod,pushoff_geod,gamI]=Geod_Closed3D_Frenet(p1,p2)

% Input:  Two closed space curves p1 and p2, entered as 3xn matrices of
%         sampled points
% Output: dist = geodesic distance between the (relative) frenet framed curves, 
%                optimized over rotations, reparameterizations.   
%         geod = The full geodesic between the base curves.
%         pushoff_geod = geodesic between the frames, viewed as pushoffs of
%                        the curves
%         gamI = Optimal reparameterization   

% What displays you want to see, set to 1 if you want to see figures. Set
% Align_center_mass=1 if you wish to display the registered curves as
% aligned with centers of mass at the origin. Otherwise they will both
% start at the origin.
    Disp_geodesic_between_curves = 0;
    Disp_initial_registration = 0;
    Disp_registration_between_curves = 0;
    Align_center_mass=0;
    Disp_optimal_reparameterization = 1;
    Disp_geodesic_tubed_curves = 1;
    
% Choose some parameters. These can be adjusted by the user.
    N=200; % Number of samples to take of each curve.
    k=20;   % Number of optimal seeds to take for the initial parameterizations.

% Resample the curves to have N points
    p1 = ReSampleCurve(p1,N);
    p2 = ReSampleCurve(p2,N);
    
% Closed curves may have last sample point equal to the first, which will
% cause problems with the assumption that the curves are immersed.
    if p1(:,N)==p1(:,1)
        p1(:,N)=1/2*p1(:,N)+1/2*p1(:,N-1);
    end
    
    if p2(:,N)==p2(:,1)
        p2(:,N)=1/2*p2(:,N)+1/2*p2(:,N-1);
    end
    
% Preprocess the curves to have the same length.
    p1=preprocess_curve(p1);
    p2=preprocess_curve(p2);
    
% Find the best k seeds for p2, with this initial parameterization.
    seeds=best_seeds(p1,p2,k);

% Find the optimal reparameterization from the given seeds.
    [z1best,w1best,z2best,w2best,gamI] = Optimal_Reparam(p1,p2,seeds);
    
% Find the optimized geodesic distance and compute the geodesic.

    [dist,~,~,~,~,z_geod,w_geod]=Grassmannian_Geodesic(z1best,w1best,z2best,w2best); 

% Normalize the geodesic distance. The diameter of the Grassmannian is 
% sqrt(2)*pi. We normalize it to have radius=1.
    if (dist>0.000001)
        dist = 2*dist/(sqrt(2)*pi);
    else
        dist=0;
    end

% Send the Grassmannian geodesic to a geodesic of framed loops.

[~,T]=size(z_geod);

% Write the geodesic as a complex matrix. 
Z_geod=zeros(2,N,T);

for l=1:N
    for j=1:T
        Z_geod(:,l,j)=[z_geod(l,j);w_geod(l,j)];
    end
end

% Send the complex matrix to a quaternionic representation.
quat_geod=zeros(4,N,T);

for j=1:T
    quat_geod(:,:,j)=complex_curve_to_quat(Z_geod(:,:,j));
end

% Define framed curve geodesic via the Hopf map.
geod=zeros(3,N,T);
pushoff_geod=zeros(3,N,T);
frames_geod=zeros(3,N,T);

for j=1:T
    [geod(:,:,j),frames_geod(:,:,j),pushoff_geod(:,:,j)]=quat_to_framed_curve(quat_geod(:,:,j));
end

for j=1:k
    geod(:,end,j)=[0;0;0];
end

% Create figures
    if Disp_optimal_reparameterization
        figure(100); clf;
        t=linspace(0,1,N);
        plot(gamI,'LineWidth',2)
        hold on
        plot(t)
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end
    
    if Disp_registration_between_curves
        if Align_center_mass
            p1new = p1 - repmat(mean(p1')',1,size(p1,2));
            p2new = p2new - repmat(mean(p2new')',1,size(p2new,2));
        else
            p1new=p1;
        end
        figure(102); clf;
        plot3(p1new(1,:),p1new(2,:),p1new(3,:),'LineWidth',3);
        hold on
        plot3(p2new(1,:),p2new(2,:),p2new(3,:),'LineWidth',3);
        axis equal
    end
    
    if Disp_geodesic_between_curves
        figure(103); clf;
        movie3D(geod)
    end
    
    if Disp_geodesic_tubed_curves
        figure(104); clf;
        pics_to_show = [1,10,25,35,50,50,75,85,100];
        for j=1:length(pics_to_show)
            tubeplot(geod(:,:,pics_to_show(j)),frames_geod(:,:,pics_to_show(j)),.01,20);
            axis equal
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'ztick',[])
            set(gca,'Visible','off')
            daspect([1,1,1]); view(39,48)
            pause
        end
    end
            
            
            
