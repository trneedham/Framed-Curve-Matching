function [dist,geod,pushoff_geod,q1,q2new]=OpenGeod_Opt(p1,p2)
 
%input: two curves p1 and p2. Entered as samples along the curves, as 3xn1
%       and 3xn2 matrices. It is not required for n1=n2.
%output: dist=geodesic distance between framed curve. The first curve has 
%             its Frenet framing and the second has a framing which has 
%             been optimized to minimize geodesic distance. Curves have
%             been optimally registered over rotations and
%             reparameterizations.
%        geod=The full geodesic joining the curves. This can be viewed as a
%             movie using movie(geod), or it will automatically display if
%             figs=1. 
%        pushoff_geod=geodesic joining the framings of the curves. The
%                     geodesic joining the full framed curves can be viewed
%                     using movie_frames(geod,pushoff_geod).
        
% Allow the program to access subfolders
    addContainingDirAndSubDir()

    tic
% What displays you want to see, set to 1 if you want to see figures. Set
% Align_center_mass=1 if you wish to display the registered curves as
% aligned with centers of mass at the origin. Otherwise they will both
% start at the origin.
    Disp_geodesic_between_curves = 0;
    Disp_initial_registration = 0;
    Disp_registration_between_curves = 0;
    Align_center_mass=0;
    Disp_optimal_reparameterization = 0;
    Disp_geodesic_tubed_curves = 1;

% Display initial curve registration, if desired.
    if Disp_initial_registration
        figure(101); clf;
        plot3(p1(1,:),p1(2,:),p1(3,:),'LineWidth',3);
        hold on
        plot3(p2(1,:),p2(2,:),p2(3,:),'LineWidth',3);
        axis equal
    end
    
% Resample the curves to have N points. N can be changed freely.
    N = 200;
    dt = 1/N;
    p1 = ReSampleCurve(p1,N);
    p2 = ReSampleCurve(p2,N);
    
% Preprocess the curves to have the same length=2
    p1 = preprocess_curve(p1);
    p2 = preprocess_curve(p2);
    
% In this version we will align the curves pointwise, since we don't want
% our initial alignment to depend on framing. This is done by a simple SVD
% decomposition.
    A = p1*p2';
    [U,~,V] = svd(A);
    if det(A)> 0
        Ot = U*V';
    else
        Ot = U*([V(:,1) V(:,2) -V(:,3)])';
    end
    p2 = Ot*p2;
    
% Form the quaternionic curve representation
    q1 = curve2quat_frenet(p1);
    q2 = curve2quat_frenet(p2);
    
% Renormalize the quaternionic curves
    sq_norms1=zeros(1,N);
    sq_norms2=zeros(1,N);
    
    for i=1:N
        sq_norms1(i)=norm(q1(:,i))^2;
        sq_norms2(i)=norm(q2(:,i))^2;
    end

    q1=q1/sqrt(sum(sq_norms1)*dt);
    q2=q2/sqrt(sum(sq_norms2)*dt);

% Find the optimal framing. Alternates between optimization over frame
% twisting and registration over rotations. Number of iterations can be
% changed.
    for i=1:3
       [~,q2]=optimal_framing(q1,q2);
       q2 = opt_quat_rot(q1,q2);
    end
    
% Applying optimal re-parameterization to the second curve
    [gam] = DynamicProgrammingQ2(q1,q2,0,0);
    gamI = invertGamma(gam);
    gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
    q2new = Group_Action_by_Gamma_Coord_Quaternions(q2,gamI);
     
% Renormalize q2new to lie on the L2 sphere.
    sq_norms2new=zeros(1,N);
    
    for i=1:N
        sq_norms2new(i)=norm(q2new(:,i))^2;
    end

    q2new=q2new/sqrt(sum(sq_norms2new)*dt);

% Find optimal framing and rotation once again. This is iterated until
% convergence of geodesic distance. The number of iterations can be freely
% changed.
    for i=1:3
        [~,q2new]=optimal_framing(q1,q2new);
        q2new=opt_quat_rot(q1,q2new);
    end 

% Define the newly parameterized and rotated curve and calculate the
% geodesic.
    [p2new,~,~]=quat_to_framed_curve(q2new);
    
    [geod,frames_geod,pushoff_geod,dist]=open_geodesic_frames(q1,q2new);
    
% We normalize geodesic distance by pi/2, since that is the diameter of the 
% L2 real projective space.
if dist>0.000001
    dist = dist/(pi/2);
  else
    dist=0;
end
    toc
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
        movie(geod)
    end
    
    if Disp_geodesic_tubed_curves
        figure(104); clf;
        pics_to_show = [1,25,50,75,100];
        for j=1:length(pics_to_show)
            tubeplot(geod(:,:,pics_to_show(j)),frames_geod(:,:,pics_to_show(j)),.01,8);
            axis equal
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            set(gca,'ztick',[])
            set(gca,'Visible','off')
            daspect([1,1,1]); view(39,48)
            pause
        end
    end
            
            
            