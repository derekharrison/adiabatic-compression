%Adiabatic compression of an ideal monatomic gas

%Simulation parameters
R = 0.02;                           %Particle radii
L = 4.0;                            %One half the length of one of the cube edges
N = 60;                             %Number of particles in the system

max_vel = 9.0;                      %Max value particle velocity
mass = 1.0;                         %Mass of a particle
max_t = 15;                         %Max simulation time
start_comp_t = 5;                   %Start time compression
end_comp_t = 12;                    %End time compression
v_wall = 0.2;                       %Compression rate
dt = 5e-3;                          %Timestep size

e  = 1.0;                           %Normal restitution coefficient
mu = 0.00;                          %Friction coefficient
Bo = 1.0;                           %Coefficient of tangential restitution

writevideo = false;                  %Generate video of simulation
sim_with_cube_boundary = false;     %Generate video of simulation with compression boundary

%Initialization
tic
err = 1e-10;                        %Small number
h = L - R - err;                    %Domain limit particle injection

l_wall_east = L;
l_wall_west = -L;
l_wall_south = -L;
l_wall_north = L;
l_wall_bottom = -L;
l_wall_top = L;

east_wall = N + 1;
west_wall = N + 2;
north_wall = N + 3;
south_wall = N + 4;
top_wall = N + 5;
bottom_wall = N + 6;

x = zeros(N,1);
y = zeros(N,1);
z = zeros(N,1);

I =2/5*mass*R^2;
Wx = zeros(N,1);
Wy = zeros(N,1);
Wz = zeros(N,1);

avg_kin_energy_before_compression = 0;
avg_kin_energy_after_compression = 0;

frame_counter = 0;

frameps = 80;                       %Set framerate (fps) for video
if writevideo == true
    writerObj = VideoWriter('C:\Users\d-w-h\Desktop\Home\testing_video_1.avi','Motion JPEG AVI');
    writerObj.FrameRate = frameps;
    open(writerObj);
end

%Generating initial position of first particle
x(1) = 2 * h * rand - h;
y(1) = 2 * h * rand - h;
z(1) = 2 * h * rand - h;

%Generating initial positions of remaining particles
for n = 2:N
    x(n) = 2 * h * rand - h;
    y(n) = 2 * h * rand - h;
    z(n) = 2 * h * rand - h;
    overlap = true;
    while overlap == true
        overlap = false;
        for i = 1:(n-1)
            dx = x(i) - x(n);
            dy = y(i) - y(n);
            dz = z(i) - z(n);
            if dx*dx + dy*dy + dz*dz < 1.0001*(2*R)*(2*R)
                overlap = true;
            end
        end
        if overlap == true
            x(n) = 2 * h * rand - h;
            y(n) = 2 * h * rand - h;
            z(n) = 2 * h * rand - h;
        end        
    end
end

%Check overlap
there_is_overlap = false;
for n = 1:N
    for i = (n+1):N 
        dx = x(i) - x(n);
        dy = y(i) - y(n);
        dz = z(i) - z(n);
        if dx*dx + dy*dy + dz*dz < (2*R)*(2*R)
             there_is_overlap = true;
        end    
    end
end
there_is_overlap;

%Check if in boundaries
in_boundaries = false;
for n = 1:N
    in_boundaries = (x(n) <= l_wall_east) && (x(n) >= l_wall_west && ... 
                     y(n) <= l_wall_north) && (y(n) >= l_wall_south && ...
                     z(n) <= l_wall_top) && (z(n) >= l_wall_bottom);
    if in_boundaries
        in_boundaries = true;
    end
end
in_boundaries;

%Generate initial velocities
v_x = 2*max_vel*rand(N,1) - max_vel;
v_y = 2*max_vel*rand(N,1) - max_vel;
v_z = 2*max_vel*rand(N,1) - max_vel;

%Start simulation
time = 0.0;
while time <= max_t
    
    %Calculating average kinetic energy of particle
    avg_kin_energy = 0.0;
    for n = 1:N
        v2 = v_x(n) * v_x(n) + v_y(n) * v_y(n) + v_z(n) * v_z(n);
        avg_kin_energy  = 0.5 * mass * v2 + avg_kin_energy;  
    end
    total_kin_energy = avg_kin_energy;
    avg_kin_energy = avg_kin_energy / N;

    if time < start_comp_t
        avg_kin_energy_before_compression = avg_kin_energy;
    end
    
    if time > end_comp_t
        avg_kin_energy_after_compression = avg_kin_energy;
    end
    
    %Calculate minimum collision time
    collision_with_wall = false;
    collision_with_particle = false;
    coll_time = 1e+10;
    
    for n = 1:N
        %Checking collision time between particles
        for i = (n+1):N
            rab = [x(n) - x(i);y(n) - y(i);z(n) - z(i)];
            vab = [v_x(n) - v_x(i);v_y(n) - v_y(i);v_z(n) - v_z(i)];
            Disc = (rab'*vab)^2-vab'*vab*(rab'*rab-(2*R)^2);
            coll_time_particle = dt + err;

            if Disc > 0
                coll_time_particle = ((-rab'*vab) - sqrt(Disc)) / (vab'*vab);
            end

            if coll_time_particle < coll_time && coll_time_particle >= 0
                collision_with_particle = true;
                coll_time = coll_time_particle;
                coll_partner_1 = n;
                coll_partner_2 = i;
            end
        end

        if start_comp_t <= time && time <= end_comp_t
            v_wall_east = -v_wall;
            v_wall_west = v_wall;
            v_wall_north = -v_wall;
            v_wall_south = v_wall;
            v_wall_top = -v_wall;
            v_wall_bottom = v_wall;
        else
            v_wall_east = -0;
            v_wall_west = 0;
            v_wall_north = -0;
            v_wall_south = 0;
            v_wall_top = -0;
            v_wall_bottom = 0;
        end
        
        coll_time_east   = (x(n) - l_wall_east + R) / (v_wall_east - v_x(n));
        coll_time_west   = (x(n) - l_wall_west - R) / (v_wall_west - v_x(n));
        coll_time_north  = (y(n) - l_wall_north + R) / (v_wall_north - v_y(n));
        coll_time_south  = (y(n) - l_wall_south - R) / (v_wall_south - v_y(n));
        coll_time_top    = (z(n) - l_wall_top + R) / (v_wall_top - v_z(n));
        coll_time_bottom = (z(n) - l_wall_bottom - R) / (v_wall_bottom - v_z(n));

       %Checking collision with east face        
        if coll_time_east > 0 && coll_time_east <= coll_time
            coll_time = coll_time_east;
            coll_partner_1 = n;
            coll_partner_2 = east_wall;
            collision_with_wall = true;
            collision_with_particle = false;
        end
        
        if coll_time_west > 0 && coll_time_west <= coll_time
            coll_time = coll_time_west;
            coll_partner_1 = n;
            coll_partner_2 = west_wall;
            collision_with_wall = true;
            collision_with_particle = false;        
        end
        
        if coll_time_north > 0 && coll_time_north <= coll_time
            coll_time = coll_time_north;
            coll_partner_1 = n;
            coll_partner_2 = north_wall;
            collision_with_wall = true;
            collision_with_particle = false;        
        end         
        
        if coll_time_south > 0 && coll_time_south <= coll_time
            coll_time = coll_time_south;
            coll_partner_1 = n;
            coll_partner_2 = south_wall;
            collision_with_wall = true;
            collision_with_particle = false;        
        end 
        
        if coll_time_top > 0 && coll_time_top <= coll_time
            coll_time = coll_time_top;
            coll_partner_1 = n;
            coll_partner_2 = top_wall;
            collision_with_wall = true;
            collision_with_particle = false;        
        end 

        if coll_time_bottom > 0 && coll_time_bottom <= coll_time
            coll_time = coll_time_bottom;
            coll_partner_1 = n;
            coll_partner_2 = bottom_wall;
            collision_with_wall = true;
            collision_with_particle = false;        
        end 
    end
    
    %Update position and velocity    
    if dt <= coll_time
        x = v_x * dt * (1 - err) + x;
        y = v_y * dt * (1 - err) + y;
        z = v_z * dt * (1 - err) + z;
        time = time + dt;
        
        l_wall_east = v_wall_east * dt + l_wall_east;
        l_wall_west = v_wall_west * dt + l_wall_west;
        l_wall_north = v_wall_north * dt + l_wall_north;        
        l_wall_south = v_wall_south * dt + l_wall_south;
        l_wall_top = v_wall_top * dt + l_wall_top;        
        l_wall_bottom = v_wall_bottom * dt + l_wall_bottom;

    elseif coll_time < dt    
        x = v_x * coll_time * (1 - err) + x;
        y = v_y * coll_time * (1 - err) + y;
        z = v_z * coll_time * (1 - err) + z;
        time = time + coll_time;
        
        l_wall_east = v_wall_east * coll_time + l_wall_east;
        l_wall_west = v_wall_west * coll_time + l_wall_west;
        l_wall_north = v_wall_north * coll_time + l_wall_north;        
        l_wall_south = v_wall_south * coll_time + l_wall_south;
        l_wall_top = v_wall_top * coll_time + l_wall_top;        
        l_wall_bottom = v_wall_bottom * coll_time + l_wall_bottom;

        %Update velocities of colliding particles
        if collision_with_particle == true
            ra = [x(coll_partner_1);y(coll_partner_1);z(coll_partner_1)];
            rb = [x(coll_partner_2);y(coll_partner_2);z(coll_partner_2)];
            va = [v_x(coll_partner_1);v_y(coll_partner_1);v_z(coll_partner_1)];
            vb = [v_x(coll_partner_2);v_y(coll_partner_2);v_z(coll_partner_2)];              
            n_vec = (ra-rb)/sqrt((ra-rb)'*(ra-rb));
            vab=[v_x(coll_partner_1)-v_x(coll_partner_2);
                 v_y(coll_partner_1)-v_y(coll_partner_2);
                 v_z(coll_partner_1)-v_z(coll_partner_2)];
            wa = [Wx(coll_partner_1);
                  Wy(coll_partner_1);
                  Wz(coll_partner_1)];
            wb = [Wx(coll_partner_2);
                  Wy(coll_partner_2);
                  Wz(coll_partner_2)]; 

            RiWi = R*wa + R*wb;
            crossRiWin = cross(RiWi, n_vec);
            vab = vab - crossRiWin;

            B1 = 7/2*(1/mass+1/mass);
            B2 = 1/mass+1/mass;    

            t = (vab - n_vec*(vab'*n_vec)) / sqrt((vab - n_vec*(vab'*n_vec))'*(vab - n_vec*(vab'*n_vec)) + err);        
            Jn = -(1+e)*(vab'*n_vec)/B2;

            stickyslide = (1+Bo)*(vab'*t)/Jn/B1;

            if mu < stickyslide %Sliding 
                Jt = -mu*Jn;                
            elseif mu >= stickyslide %Sticking
                Jt = -(1+Bo)*(vab'*t)/B1;      
            end

            J = Jn*n_vec+Jt*t;

            v_x(coll_partner_1) = J(1)/mass+v_x(coll_partner_1);
            v_y(coll_partner_1) = J(2)/mass+v_y(coll_partner_1);
            v_z(coll_partner_1) = J(3)/mass+v_z(coll_partner_1);

            v_x(coll_partner_2) = -J(1)/mass+v_x(coll_partner_2);
            v_y(coll_partner_2) = -J(2)/mass+v_y(coll_partner_2);
            v_z(coll_partner_2) = -J(3)/mass+v_z(coll_partner_2);

            crossnJ = cross(n_vec,-J);

            Wx(coll_partner_1) = -crossnJ(1)*R/I+Wx(coll_partner_1);
            Wy(coll_partner_1) = -crossnJ(2)*R/I+Wy(coll_partner_1);
            Wz(coll_partner_1) = -crossnJ(3)*R/I+Wz(coll_partner_1);

            Wx(coll_partner_2) = -crossnJ(1)*R/I+Wx(coll_partner_2);
            Wy(coll_partner_2) = -crossnJ(2)*R/I+Wy(coll_partner_2);      
            Wz(coll_partner_2) = -crossnJ(3)*R/I+Wz(coll_partner_2);

        elseif collision_with_wall == true   
            if coll_partner_2 == east_wall
                v_x(coll_partner_1) = 2 * v_wall_east - v_x(coll_partner_1);
            elseif coll_partner_2 == west_wall
                v_x(coll_partner_1) = 2 * v_wall_west - v_x(coll_partner_1);
            elseif coll_partner_2 == north_wall
                v_y(coll_partner_1) = 2 * v_wall_north - v_y(coll_partner_1);
            elseif coll_partner_2 == south_wall
                v_y(coll_partner_1) = 2 * v_wall_south - v_y(coll_partner_1);           
            elseif coll_partner_2 == top_wall
                v_z(coll_partner_1) = 2 * v_wall_top - v_z(coll_partner_1);
            elseif coll_partner_2 == bottom_wall
                v_z(coll_partner_1) = 2 * v_wall_bottom - v_z(coll_partner_1);
            end   
        end
    end
    
    %Save frames at regular intervals spaced at dt
    if (frame_counter == floor(time/dt)) && (writevideo == true)
        if sim_with_cube_boundary == true
            close all
            hold on

            patch([l_wall_west,   l_wall_east,   l_wall_east,   l_wall_west], ...
                  [l_wall_south,  l_wall_south,  l_wall_north,  l_wall_north], ...
                  [l_wall_bottom, l_wall_bottom, l_wall_bottom, l_wall_bottom], ...
                  'blue', ...
                  'FaceAlpha', 0.1)

            patch([l_wall_west,   l_wall_east,   l_wall_east,  l_wall_west], ...
                  [l_wall_south,  l_wall_south,  l_wall_south, l_wall_south], ...
                  [l_wall_bottom, l_wall_bottom, l_wall_top,   l_wall_top], ...
                  'red', ...
                  'FaceAlpha', 0.1)

            patch([l_wall_east,   l_wall_east,  l_wall_east,  l_wall_east], ...
                  [l_wall_south,  l_wall_south, l_wall_north, l_wall_north], ...
                  [l_wall_bottom, l_wall_top,   l_wall_top,   l_wall_bottom], ...
                  'blue', ...
                  'FaceAlpha', 0.1)  

            patch([l_wall_east,   l_wall_east,  l_wall_west,  l_wall_west], ...
                  [l_wall_north,  l_wall_north, l_wall_north, l_wall_north], ...
                  [l_wall_bottom, l_wall_top,   l_wall_top,   l_wall_bottom], ...
                  'red', ...
                  'FaceAlpha', 0.1) 

            patch([l_wall_west,   l_wall_west,  l_wall_west,  l_wall_west], ...
                  [l_wall_south,  l_wall_south, l_wall_north, l_wall_north], ...
                  [l_wall_bottom, l_wall_top,   l_wall_top,   l_wall_bottom], ...
                  'blue', ...
                  'FaceAlpha', 0.1) 

            patch([l_wall_west,  l_wall_west,  l_wall_east,  l_wall_east], ...
                  [l_wall_south, l_wall_north, l_wall_north, l_wall_south], ...
                  [l_wall_top,   l_wall_top,   l_wall_top,   l_wall_top], ...
                  'red', ...
                  'FaceAlpha', 0.1)

            rotate3d
            grid on
            view(3)

            scatter3(x,y,z, 2 * R * 4000, 'k', 'filled')
            xlim([-L L])
            ylim([-L L])
            zlim([-L L])
            hold off 
        else
            scatter3(x,y,z, 2 * R * 4000, 'k', 'filled')
            xlim([-L L])
            ylim([-L L])
            zlim([-L L])
        end
        
        if writevideo == true 
            set(gca,'nextplot','replacechildren');
            set(gcf,'Renderer','zbuffer') 
            frame = getframe(gcf); 
            writeVideo(writerObj,frame);
        end
        frame_counter = frame_counter + 1;
    end
    TIME = time
end
        
l_wall_east
temperature_ratio_theory = (((2*L)^3) / ((2*l_wall_east)^3))^(2/3)
temperature_ratio_simulation = avg_kin_energy_after_compression / avg_kin_energy_before_compression

if writevideo == true
    close(writerObj);
end

toc
