# Export a function mesti_subpixel_smoothing()
export mesti_subpixel_smoothing

    """
    MESTI_SUBPIXEL_SMOOTHING utilize GeometryPrimitives module and Kottke's algorithm (PRE 77, 036611 (2008) and Opt. Lett. 34, 2778(2009)) to do subpixel smoothing.
        Now it can be applied to 2D (TM or TE) and 3D cases. This is the 3D case
        === Input Arguments ===
        delta_x (positive scalar; required):
            Discretization grid size.
        domain (Cuboid structre):
            A 3D cuboid object defined the domain for subpixel smoothing
        domain_epsilon (numeric scalar; required):
            Relative permittivity of the (background) domain
        object_list::Vector{<:Shape}:
            A vector of objects, which put in the domain
        object_epsilon_list::Union{Vector{<:Int},Vector{<:Real},Vector{<:Complex}}:
            A vector of stored permittivity corresponding to the object_list
        xBC (string; required):
            Boundary condition (BC) at the two ends in x direction, effectively
            specifying Ex(n,m,l) at n=0 and n=nx_Ex+1,
                        Ey(n,m,l) at n=0 and n=nx_Ey+1,
                        Ez(n,m,l) at n=0 and n=nx_Ez+1.    
            one pixel beyond the computation domain. Available choices are:
            "Bloch"    - Ex(n+nx_Ex,m,l) = Ex(n,m,l)*exp(1i*kx_B*nx_Ex*delta_x),
                            Ey(n+nx_Ey,m,l) = Ey(n,m,l)*exp(1i*kx_B*nx_Ey*delta_x),
                            Ez(n+nx_Ez,m,l) = Ez(n,m,l)*exp(1i*kx_B*nx_Ez*delta_x).   
            "periodic" - equivalent to "Bloch" with kx_B = 0
            "PEC"      - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = Ez(nx_Ex,m,l),
                            Ey(0,m,l) = Ey(nx_Ey+1,m,l) = 0,
                            Ez(0,m,l) = Ez(nx_Ez+1,m,l) = 0.   
            "PMC"      - Ex(0,m,l) = Ex(nx_Ex+1,m,l) = 0,
                            Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                            Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l).    
            "PECPMC"   - Ex(0,m,l) = Ex(1,m,l); Ex(nx_Ex+1,m,l) = 0,
                            Ey(0,m,l) = 0; Ey(nx_Ey+1,m,l) = Ey(nx_Ey,m,l),
                            Ez(0,m,l) = 0; Ez(nx_Ez+1,m,l) = Ez(nx_Ez,m,l),    
            "PMCPEC"   - Ex(0,m,l) = 0; Ex(nx_Ex+1,m,l) = Ex(nx_Ex,m,l),
                            Ey(0,m,l) = Ey(1,m,l); Ey(nx_Ey+1,m,l) = 0,
                            Ez(0,m,l) = Ez(1,m,l); Ez(nx_Ez+1,m,l) = 0.
            where PEC stands for perfect electric conductor and PMC stands for perfect
            magnetic conductor.
        yBC (string; optional):
            Boundary condition in y direction, analogous to xBC.
        zBC (string; optional):
            Boundary condition in z direction, analogous to xBC.
        without_sb (boolean scalar; optional, defaults to false):
            Whether or not to turn off subpixel smoothing and output permittivity tensor without smoothing.
        === Output Arguments ===
        epsilon_xx (numeric array or matrix, real or complex):
            An nx_Ex-by-ny_Ex-by-nz_Ex array discretizing the relative permittivity
            profile epsilon_xx(x,y,z). Specifically, epsilon_xx(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the point (x_{n+0.5}, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_xy (numeric array or nothing, real or complex)
            An nx_Ez-by-ny_Ez-by-nz_Ex array discretizing the relative permittivity
            profile epsilon_xy(x,y,z). Specifically, epsilon_xy(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_xz (numeric array or nothing, real or complex):    
            An nx_Ey-by-ny_Ex-by-nz_Ey array discretizing the relative permittivity
            profile epsilon_xz(x,y,z). Specifically, epsilon_xz(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_yx (numeric array or nothing, real or complex):    
            An nx_Ez-by-ny_Ez-by-nz_Ey array discretizing the relative permittivity
            profile epsilon_yx(x,y,z). Specifically, epsilon_yx(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_yy (numeric array or nothing, real or complex):    
            An nx_Ey-by-ny_Ey-by-nz_Ey array discretizing the relative permittivity
            profile epsilon_yy(x,y,z). Specifically, epsilon_yy(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the point (x_n, y_{m+0.5}, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_yz (numeric array or nothing, real or complex):    
            An nx_Ey-by-ny_Ex-by-nz_Ex array discretizing the relative permittivity
            profile epsilon_yz(x,y,z). Specifically, epsilon_yz(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_zx (numeric array or nothing, real or complex):    
            An nx_Ey-by-ny_Ez-by-nz_Ey array discretizing the relative permittivity
            profile epsilon_zx(x,y,z). Specifically, epsilon_zx(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_zy (numeric array or nothing, real or complex):    
            An nx_Ez-by-ny_Ex-by-nz_Ex array discretizing the relative permittivity
            profile epsilon_zy(x,y,z). Specifically, epsilon_zy(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the low corner (x_n, y_m, z_l) on the 
            Yee lattice (n,m,l).
        epsilon_zz (numeric array or nothing, real or complex):    
            An nx_Ez-by-ny_Ez-by-nz_Ez array discretizing the relative permittivity
            profile epsilon_zz(x,y,z). Specifically, epsilon_zz(n,m,l) is averaged over a 
            cube with volume (delta_x)^3 centered at the point (x_n, y_m, z_{l+0.5}) on the 
            Yee lattice (n,m,l).
    """
    function mesti_subpixel_smoothing(delta_x::Real, domain::Cuboid{3, 9}, domain_epsilon::Union{Int,Real,Complex}, object_list::Vector{<:Shape}, object_epsilon_list::Union{Vector{<:Int},Vector{<:Real},Vector{<:Complex}}, xBC::String, yBC::String, zBC::String, without_sb::Bool=false)
        
        if ~isapprox((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x, round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)) 
            # If domain length along x is not multiple of delta_x, we adjust the domain size based on delta_x.
            @warn "Domain length along x is not multiple of delta_x. We adjust the upper bound in the domain from $(bounds(domain)[2][1]) to $(max(round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x,delta_x))"
            # The max() standing for length of domain along x should at least be delta_x        
            domain = Cuboid(
         [max((round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x)/2,delta_x/2),(bounds(domain)[2][2]-bounds(domain)[1][2])/2,(bounds(domain)[2][3]-bounds(domain)[1][3])/2]
            , [max((round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x),delta_x),(bounds(domain)[2][2]-bounds(domain)[1][2]),(bounds(domain)[2][3]-bounds(domain)[1][3])])
        end
        
        if ~isapprox((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x, round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)) 
            # If domain length along y is not multiple of delta_x, we adjust the domain size based on delta_x.
            @warn "Domain length along y is not multiple of delta_x. We adjust the upper bound in the domain from $(bounds(domain)[2][2]) to $(max(round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x,delta_x))"
            # The max() standing for length of domain along y should at least be delta_x        
            domain = Cuboid(
     [(bounds(domain)[2][1]-bounds(domain)[1][1])/2,max((round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x)/2,delta_x/2),(bounds(domain)[2][3]-bounds(domain)[1][3])/2]
        , [(bounds(domain)[2][1]-bounds(domain)[1][1]),max((round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x),delta_x),(bounds(domain)[2][3]-bounds(domain)[1][3])])
        end
    
            if ~isapprox((bounds(domain)[2][3] - bounds(domain)[1][3])/delta_x, round((bounds(domain)[2][3] - bounds(domain)[1][3])/delta_x)) 
            # If domain length along z is not multiple of delta_x, we adjust the domain size based on delta_x.
            @warn "Domain length along z is not multiple of delta_x. We adjust the upper bound in the domain from $(bounds(domain)[2][3]) to $(max(round((bounds(domain)[2][3] - bounds(domain)[1][3])/delta_x)*delta_x,delta_x))"
            # The max() standing for length of domain along z should at least be delta_x        
            domain = Cuboid(
     [(bounds(domain)[2][1]-bounds(domain)[1][1])/2,(bounds(domain)[2][2]-bounds(domain)[1][2])/2,max((round((bounds(domain)[2][3] - bounds(domain)[1][3])/delta_x)*delta_x)/2,delta_x/2)]
        , [(bounds(domain)[2][1]-bounds(domain)[1][1]),(bounds(domain)[2][2]-bounds(domain)[1][2]),max((round((bounds(domain)[2][3] - bounds(domain)[1][3])/delta_x)*delta_x),delta_x)])
        end

        # If the domain does not start from (0,0,0), translate the coordinates of the domain and obejcts.
        if bounds(domain)[1][1] != 0 || bounds(domain)[1][2] != 0 || bounds(domain)[1][3] != 0
            for obj_ind = 1:length(object_list)
                object_list[obj_ind] = translate(object_list[obj_ind], [-bounds(domain)[1][1], -bounds(domain)[1][2], -bounds(domain)[1][3]])
            end
            domain = translate(domain, [-bounds(domain)[1][1], -bounds(domain)[1][2], -bounds(domain)[1][3]])
        end
            
        # Check the boundary conditions in all directions
        xBC = convert_BC_sbpsm(xBC, "x")
        yBC = convert_BC_sbpsm(yBC, "y")
        zBC = convert_BC_sbpsm(zBC, "z")

        # If the boundary conidition is Bloch periodic along x, check x coordinate in every object. When starting  coordinate or ending coordinate of a object is outside the domain, but the other is not, insert the periodic image along x of this object into the object list.
        if xBC == "periodic" || xBC == "Bloch"
            object_list_temp = deepcopy(object_list)        
            num_new_periodic_object = 0
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                if bounds(object)[1][1] <= bounds(domain)[1][1] && bounds(object)[2][1] >= bounds(domain)[1][1] && bounds(object)[2][1] < bounds(domain)[2][1]
                    insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [bounds(domain)[2][1]-bounds(domain)[1][1],0,0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])
                    num_new_periodic_object = num_new_periodic_object + 1
                elseif bounds(object)[1][1] <= bounds(domain)[2][1] && bounds(object)[2][1] >= bounds(domain)[2][1] && bounds(object)[1][1] > bounds(domain)[1][1]
                     insert!(object_list_temp,  translate(object, obj_ind+num_new_periodic_object+1, [-(bounds(domain)[2][1]-bounds(domain)[1][1]),0,0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                
                    num_new_periodic_object = num_new_periodic_object + 1               
                end
            end
            object_list = deepcopy(object_list_temp)
            object_list_temp = nothing
        end

        # If the boundary conidition is Bloch periodic along y, check y coordinate in every object. When starting  coordinate or ending coordinate of a object is outside the domain, but the other is not, insert the periodic image along y of this object into the object list.
        if yBC == "periodic" || yBC == "Bloch"
            object_list_temp = deepcopy(object_list)        
            num_new_periodic_object = 0        
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                if bounds(object)[1][2] <= bounds(domain)[1][2] && bounds(object)[2][2] >= bounds(domain)[1][2] && bounds(object)[2][2] < bounds(domain)[2][2]
                    insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [0, bounds(domain)[2][2]-bounds(domain)[1][2], 0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                
                    num_new_periodic_object = num_new_periodic_object + 1
                elseif bounds(object)[1][2] <= bounds(domain)[2][2] && bounds(object)[2][2] >= bounds(domain)[2][2] && bounds(object)[1][2] > bounds(domain)[1][2]
                     insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [0, -(bounds(domain)[2][2]-bounds(domain)[1][2]), 0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                
                    num_new_periodic_object = num_new_periodic_object + 1               
                end
            end
            object_list = deepcopy(object_list_temp)
            object_list_temp = nothing        
        end

        # If the boundary conidition is Bloch periodic along z, check z coordinate in every object. When starting  coordinate or ending coordinate of a object is outside the domain, but the other is not, insert the periodic image along z of this object into the object list.    
        if zBC == "periodic" || zBC == "Bloch"
            object_list_temp = deepcopy(object_list)        
            num_new_periodic_object = 0        
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                if bounds(object)[1][3] <= bounds(domain)[1][3] && bounds(object)[2][3] >= bounds(domain)[1][3] && bounds(object)[2][3] < bounds(domain)[2][3]
                    insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [0, 0, bounds(domain)[2][3]-bounds(domain)[1][3]]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                
                    num_new_periodic_object = num_new_periodic_object + 1
                elseif bounds(object)[1][3] <= bounds(domain)[2][3] && bounds(object)[2][3] >= bounds(domain)[2][3] && bounds(object)[1][3] > bounds(domain)[1][3]
                     insert!(object_list_temp, object, obj_ind+num_new_periodic_object+1, translate(object, [0, 0, -(bounds(domain)[2][3]-bounds(domain)[1][3])]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                
                    num_new_periodic_object = num_new_periodic_object + 1               
                end
            end
            object_list = deepcopy(object_list_temp)
            object_list_temp = nothing        
        end
    
        # 3D case
        # Coordinate for Ex-site
        Ex_x_coord = (bounds(domain)[1][1]+delta_x/2):delta_x:(bounds(domain)[2][1])
        Ex_y_coord = (bounds(domain)[1][2]):delta_x:(bounds(domain)[2][2]-delta_x/2)
        Ex_z_coord = (bounds(domain)[1][3]):delta_x:(bounds(domain)[2][3]-delta_x/2)

        # Coordinate for Ey-site   
        Ey_x_coord = (bounds(domain)[1][1]):delta_x:(bounds(domain)[2][1]-delta_x/2)
        Ey_y_coord = (bounds(domain)[1][2]+delta_x/2):delta_x:(bounds(domain)[2][2])
        Ey_z_coord = Ex_z_coord

        # Coordinate for Ez-site
        Ez_x_coord = Ey_x_coord
        Ez_y_coord = Ex_y_coord
        Ez_z_coord = (bounds(domain)[1][3]+delta_x/2):delta_x:(bounds(domain)[2][3])

        # Coordinate for Eo-site (lower corner of a Yee cell)    
        Eo_x_coord = Ey_x_coord
        Eo_y_coord = Ex_y_coord
        Eo_z_coord = Ex_z_coord        

        # the number of Yee grid                
        Nx = length(Eo_x_coord)
        Ny = length(Eo_y_coord)
        Nz = length(Eo_z_coord)
        
        common_type = promote_type(typeof(domain_epsilon*1.0), eltype(object_epsilon_list))
    
        # initializing epsilon on Eo, Ex, Ey, and Ez sites.    
        epsilon_Eo_site = Array{common_type}(undef,Nx,Ny,Nz,3,3)
        epsilon_Ex_site = Array{common_type}(undef,Nx,Ny,Nz,3,3)
        epsilon_Ey_site = Array{common_type}(undef,Nx,Ny,Nz,3,3)
        epsilon_Ez_site = Array{common_type}(undef,Nx,Ny,Nz,3,3)
    
        for k = 1:Nz, j = 1:Ny, i = 1:Nx
            epsilon_Eo_site[i,j,k,:,:] = domain_epsilon * Matrix(I,3,3)
            epsilon_Ex_site[i,j,k,:,:] = domain_epsilon * Matrix(I,3,3)
            epsilon_Ey_site[i,j,k,:,:] = domain_epsilon * Matrix(I,3,3)
            epsilon_Ez_site[i,j,k,:,:] = domain_epsilon * Matrix(I,3,3)        
        end

        for obj_ind = 1:length(object_epsilon_list)
            object = object_list[obj_ind]        
            # looping over the range of the object in Eo grids.
            # Note that adding eps() can help avoid non-ideal situations where the last pixel might not be properly counted due to floating-point precision issues.
            for kk = max((Int(div(bounds(object)[1][3]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][3]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Ny), ii = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Nx)
                ix = Eo_x_coord[ii]; iy = Eo_y_coord[jj]; iz = Eo_z_coord[kk]
                if without_sb
                    if [ix,iy,iz] ∈ object
                        epsilon_Eo_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)
                    end
                else
                    if all([[ix+dx,iy+dy,iz+dz] ∈ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dx = (-delta_x/2,delta_x/2)])
                        # If the eight corner of the voxel all are inside the object, the whole pixel should be occupied by the object.
                        epsilon_Eo_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)
                    elseif ~(all([[ix+dx,iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dx = (-delta_x/2,delta_x/2)]))
                        # The voxel is partially occupied by the object
                        # Calculating the surface point
                        # r0: position of the surface point
                        # n0: unit vector from point [ix,iy,iz] to r0
                        (r0,n0) = surfpt_nearby([ix,iy,iz], object)
                        vxl = (SVector(ix-delta_x/2,iy-delta_x/2,iz-delta_x/2), 
                               SVector(ix+delta_x/2,iy+delta_x/2,iz+delta_x/2))
                        # vol_frac: volume filling fraction of the object
                        vol_frac = volfrac(vxl, n0, r0)
                        epsilon_Eo_site[ii,jj,kk,:,:] = Kottke_smoothing(vol_frac, n0, object_epsilon_list[obj_ind] * Matrix(I,3,3), epsilon_Eo_site[ii,jj,kk,:,:])
                    end 
                end
            end

            # looping over the range of the object in Ex grids.                                    
            for kk = max((Int(div(bounds(object)[1][3]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][3]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Ny), ii = max((Int(div(bounds(object)[1][1]-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+eps(), delta_x)) + 1),Nx)               
                ix = Ex_x_coord[ii]; iy = Ex_y_coord[jj]; iz = Ex_z_coord[kk]
                if without_sb
                    if [ix,iy,iz] ∈ object
                        epsilon_Ex_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)
                    end
                else            
                    if all([[ix+dx,iy+dy,iz+dz] ∈ object for dx = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dz = (-delta_x/2,delta_x/2)])
                        epsilon_Ex_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)                
                    elseif ~(all([[ix+dx,iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dx = (-delta_x/2,delta_x/2)]))
                        (r0,n0) = surfpt_nearby([ix,iy,iz], object)
                        vxl = (SVector(ix-delta_x/2,iy-delta_x/2,iz-delta_x/2), 
                               SVector(ix+delta_x/2,iy+delta_x/2,iz+delta_x/2))
                        vol_frac = volfrac(vxl, n0, r0)
                        epsilon_Ex_site[ii,jj,kk,:,:] = Kottke_smoothing(vol_frac, n0, object_epsilon_list[obj_ind] * Matrix(I,3,3), epsilon_Ex_site[ii,jj,kk,:,:])
                    end 
                end
            end

            # looping over the range of the object in Ey grids.
            for kk = max((Int(div(bounds(object)[1][3]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][3]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][2]-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+eps(), delta_x)) + 1),Ny), ii = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Nx)            
                ix = Ey_x_coord[ii]; iy = Ey_y_coord[jj]; iz = Ey_z_coord[kk]
                if without_sb
                    if [ix,iy,iz] ∈ object
                        epsilon_Ey_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)
                    end
                else          
                    if all([[ix+dx,iy+dy,iz+dz] ∈ object for dx = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dz = (-delta_x/2,delta_x/2)])
                        epsilon_Ey_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)                
                    elseif ~(all([[ix+dx,iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dx = (-delta_x/2,delta_x/2)]))
                        (r0,n0) = surfpt_nearby([ix,iy,iz], object)
                        vxl = (SVector(ix-delta_x/2,iy-delta_x/2,iz-delta_x/2), 
                               SVector(ix+delta_x/2,iy+delta_x/2,iz+delta_x/2))
                        vol_frac = volfrac(vxl, n0, r0)
                        epsilon_Ey_site[ii,jj,kk,:,:] = Kottke_smoothing(vol_frac, n0, object_epsilon_list[obj_ind] * Matrix(I,3,3), epsilon_Ey_site[ii,jj,kk,:,:])
                    end 
                end
            end

            # looping over the range of the object in Ez grids.
            for kk = max((Int(div(bounds(object)[1][3]-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][3]+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Ny), ii = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Nx)
                ix = Ez_x_coord[ii]; iy = Ez_y_coord[jj]; iz = Ez_z_coord[kk]
                if without_sb
                    if [ix,iy,iz] ∈ object
                        epsilon_Ez_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)
                    end
                else            
                    if all([[ix+dx,iy+dy,iz+dz] ∈ object for dx = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dz = (-delta_x/2,delta_x/2)])
                        epsilon_Ez_site[ii,jj,kk,:,:] = object_epsilon_list[obj_ind] * Matrix(I,3,3)                
                    elseif ~(all([[ix+dx,iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2), dx = (-delta_x/2,delta_x/2)]))
                        (r0,n0) = surfpt_nearby([ix,iy,iz], object)
                        vxl = (SVector(ix-delta_x/2,iy-delta_x/2,iz-delta_x/2), 
                               SVector(ix+delta_x/2,iy+delta_x/2,iz+delta_x/2))
                        vol_frac = volfrac(vxl, n0, r0)
                        epsilon_Ez_site[ii,jj,kk,:,:] = Kottke_smoothing(vol_frac, n0, object_epsilon_list[obj_ind] * Matrix(I,3,3), epsilon_Ez_site[ii,jj,kk,:,:])                
                    end
                end
            end    
        end
    
        # Pick the pixel to be returned based on boundary condiions.
        return pick_epsilon_3d(epsilon_Eo_site, epsilon_Ex_site, epsilon_Ey_site, epsilon_Ez_site, xBC, yBC, zBC)               
    end    
    
    """
    MESTI_SUBPIXEL_SMOOTHING utilize GeometryPrimitives module  and Kottke's algorithm (PRE 77, 036611 (2008) and Opt. Lett. 34, 2778(2009)) to do subpixel smoothing.
        Now it can be applied to 2D (TM or TE) and 3D cases. This is the 2D case
        === Input Arguments ===
        delta_x (positive scalar; required):
            Discretization grid size.
        domain (Cuboid structre):
            A 2D cuboid object defined the domain for subpixel smoothing
        domain_epsilon (numeric scalar; required):
            Relative permittivity of the (background) domain
        object_list::Vector{<:Shape}:
            A vector of objects, which put in the domain
        object_epsilon_list::Union{Vector{<:Int},Vector{<:Real},Vector{<:Complex}}:
            A vector of stored permittivity corresponding to the object_list
        yBC (string; optional):
            Boundary condition (BC) at the two ends in y direction, in 2D coordinate (y,z) effectively
            specifying Ex(m,l) at m=0 and m=ny_Ex+1,
                    Ey(m,l) at m=0 and m=ny_Ey+1,
                    Ez(m,l) at m=0 and m=ny_Ez+1.    
            one pixel beyond the computation domain. Available choices are:
            "Bloch"    - Ex(m+ny_Ex,l) = Ex(m,l)*exp(1i*ky_B*ny_Ex*delta_x),
                        Ey(m+ny_Ey,l) = Ey(m,l)*exp(1i*ky_B*ny_Ey*delta_x),
                        Ez(m+ny_Ez,l) = Ez(m,l)*exp(1i*ky_B*ny_Ez*delta_x).   
            "periodic" - equivalent to "Bloch" with ky_B = 0
            "PEC"      - Ex(0,l) = Ex(ny_Ex+1,l) = 0,
                        Ey(0,l) = Ey(1,l); Ey(ny_Ey+1,l) = Ey(ny_Ey,l),
                        Ez(0,l) = Ez(ny_Ez+1,l) = 0.   
            "PMC"      - Ex(0,l) = Ex(1,l); Ex(ny_Ex+1,l) = Ex(ny_Ex,l),
                        Ey(0,l) = Ey(ny_Ey+1,l) = 0,
                        Ez(0,l) = Ez(1,l); Ez(ny_Ez+1,l) = Ez(ny_Ez,l).   
            "PECPMC"   - Ex(0,l) = 0; Ex(ny_Ex+1,l) = Ex(ny_Ex,l),
                        Ey(0,l) = Ey(1,l); Ey(ny_Ey+1,l) = 0,
                        Ez(0,l) = 0; Ez(ny_Ez+1,l) = Ez(ny_Ez,l),   
            "PMCPEC"   - Ex(0,l) = Ex(1,l); Ex(ny_Ex+1,l) = 0,
                        Ey(0,l) = 0; Ey(ny_Ey+1,l) = Ey(ny_Ey,l),
                        Ez(0,l) = Ez(1,l); Ez(ny_Ez+1,l) = 0.
            where PEC stands for perfect electric conductor and PMC stands for perfect
            magnetic conductor.
        zBC (string; optional):
            Boundary condition in z direction, analogous to yBC.

        use_2D_TM (boolean scalar; required, defaults to false):
            Whether to output permittivity for 2D TM waves.
        use_2D_TE (boolean scalar; required, defaults to false):
            Whether to output inverse permittivity for 2D TE waves.
        without_sb (boolean scalar; optional, defaults to false):
            Whether or not to turn off subpixel smoothing and output permittivity tensor without smoothing.
        === Output Arguments ===
        epsilon_xx (numeric array or matrix, real or complex, when use_2D_TM = true):
            An ny_Ex-by-nz_Ex array discretizing the relative permittivity
            profile epsilon_xx(y,z). Specifically, epsilon_xx(m,l) is averaged over a 
            square with area (delta_x)^2 centered at the point (y_m, z_l) on the 
            Yee lattice (m,l).
        inv_epsilon_yy (numeric array or nothing, real or complex, when use_2D_TE = true):    
            An ny_Ey-by-nz_Ey array discretizing the yy components of the discretized inverse 
            relative permittivity 1/epsilon(y,z) tensor. Specifically, inv_epsilon_yy(m,l) centered 
            at the point (y_{m+0.5}, z_l) on the Yee lattice (m,l).
        inv_epsilon_zz (numeric array or nothing, real or complex, when use_2D_TE = true):    
            An ny_Ez-by-nz_Ez array discretizing the zz components of the discretized inverse 
            relative permittivity 1/epsilon(y,z) tensor. Specifically, inv_epsilon_zz(m,l) centered 
            at the point (y_m, z_{l+0.5}) on the Yee lattice (m,l).
        inv_epsilon_yz (numeric array or nothing, real or complex, when use_2D_TE = true):    
            An ny_Ex-by-nz_Ex array discretizing the yz components of the discretized inverse 
            relative permittivity 1/epsilon(y,z) tensor. Specifically, inv_epsilon_yz(m,l) centered 
            at the point (y_m, z_l) on the Yee lattice (m,l).
    """
    function mesti_subpixel_smoothing(delta_x::Real, domain::Cuboid{2, 4}, domain_epsilon::Union{Int,Real,Complex}, object_list::Vector{<:Shape}, object_epsilon_list::Union{Vector{<:Int},Vector{<:Real},Vector{<:Complex}}, yBC::String, zBC::String, use_2D_TM::Bool=true, use_2D_TE::Bool=false, without_sb::Bool=false)
    
        if ~isapprox((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x, round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)) 
            # If domain length along y is not multiple of delta_x, we adjust the domain size based on delta_x.
            @warn "Domain length along x is not multiple of delta_x. We adjust the upper bound in the domain from $(bounds(domain)[2][1]) to $(max(round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x,delta_x))"
            # The max() standing for length of domain along y should at least be delta_x
            domain = Cuboid(
         [max((round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x)/2,delta_x/2),(bounds(domain)[2][2]-bounds(domain)[1][2])/2]
            , [max((round((bounds(domain)[2][1] - bounds(domain)[1][1])/delta_x)*delta_x),delta_x),(bounds(domain)[2][2]-bounds(domain)[1][2])])
        end
        
        if ~isapprox((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x, round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)) 
            # If domain length along z is not multiple of delta_x, we adjust the domain size based on delta_x.
            @warn "Domain length along y is not multiple of delta_x. We adjust the upper bound in the domain from $(bounds(domain)[2][2]) to $(max(round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x,delta_x))"
            # The max() standing for length of domain along z should at least be delta_x        
            domain = Cuboid(
     [(bounds(domain)[2][1]-bounds(domain)[1][1])/2,max((round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x)/2,delta_x/2)]
        , [(bounds(domain)[2][1]-bounds(domain)[1][1]),max((round((bounds(domain)[2][2] - bounds(domain)[1][2])/delta_x)*delta_x),delta_x)])
        end
        
        if (~use_2D_TM && ~use_2D_TE)
            throw(ArgumentError("In 2D case, use_2D_TM and/or use_2D_TE should be true"));
        end
    
        # If the domain does not start from (0,0,0), translate the coordinates of the domain and obejcts.
        if bounds(domain)[1][1] != 0 || bounds(domain)[1][2] != 0
            for obj_ind = 1:length(object_list)
                object_list[obj_ind] = translate(object_list[obj_ind], [-bounds(domain)[1][1], -bounds(domain)[1][2]])
            end
            domain = translate(domain, [-bounds(domain)[1][1], -bounds(domain)[1][2]])
        end
    
        # Check the boundary conditions in all directions
        yBC = convert_BC_sbpsm(yBC, "y")
        zBC = convert_BC_sbpsm(zBC, "z")
     
        # If the boundary conidition is Bloch periodic along y, check y coordinate in every object. When starting  coordinate or ending coordinate of a object is outside the domain, but the other is not, insert the periodic image along y of this object into the object list.    
        if yBC == "periodic" || yBC == "Bloch"
            object_list_temp = deepcopy(object_list)        
            num_new_periodic_object = 0
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                if bounds(object)[1][1] <= bounds(domain)[1][1] && bounds(object)[2][1] >= bounds(domain)[1][1] && bounds(object)[2][1] < bounds(domain)[2][1]
                    insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [bounds(domain)[2][1]-bounds(domain)[1][1],0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                                
                    num_new_periodic_object = num_new_periodic_object + 1
                elseif bounds(object)[1][1] <= bounds(domain)[2][1] && bounds(object)[2][1] >= bounds(domain)[2][1] && bounds(object)[1][1] > bounds(domain)[1][1]
                     insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [-(bounds(domain)[2][1]-bounds(domain)[1][1]),0]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                
                    num_new_periodic_object = num_new_periodic_object + 1               
                end
            end
            object_list = deepcopy(object_list_temp)
            object_list_temp = nothing
        end

        # If the boundary conidition is Bloch periodic along z, check z coordinate in every object. When starting  coordinate or ending coordinate of a object is outside the domain, but the other is not, insert the periodic image along z of this object into the object list.        
        if zBC == "periodic" || zBC == "Bloch"
            object_list_temp = deepcopy(object_list)        
            num_new_periodic_object = 0        
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                if bounds(object)[1][2] <= bounds(domain)[1][2] && bounds(object)[2][2] >= bounds(domain)[1][2] && bounds(object)[2][2] < bounds(domain)[2][2]
                    insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [0, bounds(domain)[2][2]-bounds(domain)[1][2]]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                                
                    num_new_periodic_object = num_new_periodic_object + 1
                elseif bounds(object)[1][2] <= bounds(domain)[2][2] && bounds(object)[2][2] >= bounds(domain)[2][2] && bounds(object)[1][2] > bounds(domain)[1][2]
                     insert!(object_list_temp, obj_ind+num_new_periodic_object+1, translate(object, [0, -(bounds(domain)[2][2]-bounds(domain)[1][2])]))
                    insert!(object_epsilon_list, obj_ind+num_new_periodic_object+1, object_epsilon_list[obj_ind+num_new_periodic_object])                                                
                    num_new_periodic_object = num_new_periodic_object + 1               
                end
            end
            object_list = deepcopy(object_list_temp)
            object_list_temp = nothing        
        end
    
        common_type = promote_type(typeof(domain_epsilon*1.0), eltype(object_epsilon_list))
    
        if use_2D_TM
            # 2D TM   
            # Coordinate for Ex-site
            Ex_y_coord = (bounds(domain)[1][1]):delta_x:(bounds(domain)[2][1]-delta_x/2)
            Ex_z_coord = (bounds(domain)[1][2]):delta_x:(bounds(domain)[2][2]-delta_x/2)
        
            # the number of Yee grid along x-direction and y-direction
            Nz = length(Ex_z_coord)
            Ny = length(Ex_y_coord)
    
            # initializing epsilon_xx
            epsilon_xx = Array{common_type}(undef,Ny,Nz)
            for k = 1:Nz, j = 1:Ny
                epsilon_xx[j,k] = domain_epsilon
            end
        
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                # looping over the range of the object in Ex grids. 
                # Note that adding eps() can help avoid non-ideal situations where the last pixel might not be properly counted due to floating-point precision issues.            
                for kk = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Ny)
                    iz = Ex_z_coord[kk]; iy = Ex_y_coord[jj]
                    if without_sb
                        if [iy,iz] ∈ object
                            epsilon_xx[jj,kk] = object_epsilon_list[obj_ind]
                        end
                    else                                
                        if all([[iy+dy,iz+dz] ∈ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)]) # Note that this line make it different from build_epsilon_disorder.m
                            # If the eight corner of the pixel all are inside the object, the whole pixel should be occupied by the object. 
                            epsilon_xx[jj,kk] = object_epsilon_list[obj_ind]
                        elseif ~(all([[iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)]))
                            # The voxel is partially occupied by the object
                            # Calculating the surface point
                            # r0: position of the surface point
                            # n0: unit vector from point [ix,iy,iz] to r0
                            (r0,n0) = surfpt_nearby([iy,iz], object)
                            vxl = (SVector(iy-delta_x/2,iz-delta_x/2), 
                                   SVector(iy+delta_x/2,iz+delta_x/2))
                            # vol_frac: volume filling fraction of the object
                            vol_frac = volfrac(vxl, n0, r0)
                            epsilon_xx[jj,kk] = (Kottke_smoothing(vol_frac, SVector(0,n0[1],n0[2]), object_epsilon_list[obj_ind] * Matrix(I,3,3), epsilon_xx[jj,kk] * Matrix(I,3,3)))[1,1]
                        end
                    end
                end
            end
        end
        
        
        if use_2D_TE
            # 2D TE 
            # Coordinate for Ez-site        
            Ez_y_coord = (bounds(domain)[1][1]):delta_x:(bounds(domain)[2][1]-delta_x/2)
            Ez_z_coord = (bounds(domain)[1][2]+delta_x/2):delta_x:(bounds(domain)[2][2])
        
            # Coordinate for Ey-site        
            Ey_y_coord = (bounds(domain)[1][1]+delta_x/2):delta_x:(bounds(domain)[2][1])
            Ey_z_coord = (bounds(domain)[1][2]):delta_x:(bounds(domain)[2][2]-delta_x/2)

            # Coordinate for Eo-site        
            Eo_y_coord = Ez_y_coord
            Eo_z_coord = Ey_z_coord
                
            # the number of Yee grid for x-direction and y-direction        
            Nz = length(Ez_z_coord)
            Ny = length(Ez_y_coord)     
        
            inv_epsilon_Eo_site = Array{common_type}(undef,Ny,Nz,3,3)
            inv_epsilon_Ey_site = Array{common_type}(undef,Ny,Nz,3,3)
            inv_epsilon_Ez_site = Array{common_type}(undef,Ny,Nz,3,3)
    
            for k = 1:Nz, j = 1:Ny
                inv_epsilon_Eo_site[j,k,:,:] = (domain_epsilon^(-1)) * Matrix(I,3,3)
                inv_epsilon_Ey_site[j,k,:,:] = (domain_epsilon^(-1)) * Matrix(I,3,3)
                inv_epsilon_Ez_site[j,k,:,:] = (domain_epsilon^(-1)) * Matrix(I,3,3)      
            end
            
            for obj_ind = 1:length(object_list)
                object = object_list[obj_ind]
                # looping over the range of the object in Eo grids.            
                # Note that adding eps() can help avoid non-ideal situations where the last pixel might not be properly counted due to floating-point precision issues.        
                for kk = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Ny)                           iz = Ey_z_coord[kk]; iy = Ez_y_coord[jj]
                    if without_sb
                        if [iy,iz] ∈ object
                            inv_epsilon_Eo_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        end
                    else            
                        if all([[iy+dy,iz+dz] ∈ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)])
                            # If the eight corner of the pixel all are inside the object, the whole pixel should be occupied by the object.                                        
                            inv_epsilon_Eo_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        elseif ~(all([[iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)]))                    
                            # The voxel is partially occupied by the object
                            # Calculating the surface point
                            # r0: position of the surface point
                            # n0: unit vector from point [ix,iy,iz] to r0
                            (r0,n0) = surfpt_nearby([iy,iz], object)
                            vxl = (SVector(iy-delta_x/2,iz-delta_x/2), 
                                   SVector(iy+delta_x/2,iz+delta_x/2))
                            # vol_frac: volume filling fraction of the object
                            vol_frac = volfrac(vxl, n0, r0)
                            inv_epsilon_Eo_site[jj,kk,:,:] = inv(Kottke_smoothing(vol_frac, SVector(0,n0[1],n0[2]), object_epsilon_list[obj_ind] * Matrix(I,3,3), inv(inv_epsilon_Eo_site[jj,kk,:,:])))
                        end
                    end
                end
                # looping over the range of the object in Ey grids.                        
                for kk = max((Int(div(bounds(object)[1][2]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+delta_x/2+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][1]-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+eps(), delta_x)) + 1),Ny)                   
                    iz = Ey_z_coord[kk]; iy = Ey_y_coord[jj]
                    if without_sb
                        if [iy,iz] ∈ object
                            inv_epsilon_Ey_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        end
                    else                            
                        if all([[iy+dy,iz+dz] ∈ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)])
                            inv_epsilon_Ey_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        elseif ~(all([[iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)]))
                            (r0,n0) = surfpt_nearby([iy,iz], object)
                            vxl = (SVector(iy-delta_x/2,iz-delta_x/2), 
                                   SVector(iy+delta_x/2,iz+delta_x/2))
                            # vol_frac: volume filling fraction of the object
                            vol_frac = volfrac(vxl, n0, r0)
                            inv_epsilon_Ey_site[jj,kk,:,:] = inv(Kottke_smoothing(vol_frac, SVector(0,n0[1],n0[2]), object_epsilon_list[obj_ind] * Matrix(I,3,3), inv(inv_epsilon_Ey_site[jj,kk,:,:])))

                        end
                    end
                end
                # looping over the range of the object in Ez grids.                            
                for kk = max((Int(div(bounds(object)[1][2]-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][2]+eps(), delta_x)) + 1),Nz), jj = max((Int(div(bounds(object)[1][1]+delta_x/2-eps(), delta_x)) + 1),1):min((Int(div(bounds(object)[2][1]+delta_x/2+eps(), delta_x)) + 1),Ny)            
                    iz = Ez_z_coord[kk]; iy = Ez_y_coord[jj]
                    if without_sb
                        if [iy,iz] ∈ object
                            inv_epsilon_Ez_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        end
                    else            
                        if all([[iy+dy,iz+dz] ∈ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)])
                            inv_epsilon_Ez_site[jj,kk,:,:] = ((object_epsilon_list[obj_ind]^(-1)) * Matrix(I,3,3))
                        elseif ~(all([[iy+dy,iz+dz] ∉ object for dz = (-delta_x/2,delta_x/2), dy = (-delta_x/2,delta_x/2)]))
                            (r0,n0) = surfpt_nearby([iy,iz], object)
                            vxl = (SVector(iy-delta_x/2,iz-delta_x/2), 
                                   SVector(iy+delta_x/2,iz+delta_x/2))
                            # vol_frac: volume filling fraction of the object
                            vol_frac = volfrac(vxl, n0, r0)
                            inv_epsilon_Ez_site[jj,kk,:,:] = inv(Kottke_smoothing(vol_frac, SVector(0,n0[1],n0[2]), object_epsilon_list[obj_ind] * Matrix(I,3,3), inv(inv_epsilon_Ez_site[jj,kk,:,:])))
                        end
                    end
                end
            end
        end
    
        if use_2D_TM && use_2D_TE
            # Return epsilon_xx, inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz 
            # Pick the pixel to be returned based on boundary condiions.        
            return pick_epsilon_2d_TM(epsilon_xx, yBC, zBC), pick_inv_epsilon_2d_TE(inv_epsilon_Ey_site, inv_epsilon_Ez_site, inv_epsilon_Eo_site, yBC, zBC)
        elseif use_2D_TM
            # Return epsilon_xx 
            # Pick the pixel to be returned based on boundary condiions.
            return pick_epsilon_2d_TM(epsilon_xx, yBC, zBC)
        elseif use_2D_TE
            # Return inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz
            # Pick the pixel to be returned based on boundary condiions.    
            return pick_inv_epsilon_2d_TE(inv_epsilon_Ey_site, inv_epsilon_Ez_site, inv_epsilon_Eo_site, yBC, zBC)
        end    
    end
    
    
    function tau_inverse_trans(tau)
        tau_11, tau_21, tau_31, tau_12, tau_22, tau_32, tau_13, tau_23, tau_33 = tau
        
        return [     -1/tau_11         -tau_12/tau_11                  -tau_13/tau_11;
                -tau_21/tau_11 (tau_22 - tau_21*tau_12/tau_11) (tau_23 - tau_21*tau_13/tau_11);
                -tau_31/tau_11 (tau_32 - tau_31*tau_12/tau_11) (tau_33 - tau_31*tau_13/tau_11)
               ]
    end
    
    function tau_trans(epsilon)
        epsilon_11, epsilon_21, epsilon_31, epsilon_12, epsilon_22, epsilon_32, epsilon_13, epsilon_23, epsilon_33 = epsilon
    
        return [-1/epsilon_11 epsilon_12/epsilon_11 epsilon_13/epsilon_11;
                epsilon_21/epsilon_11 (epsilon_22 - epsilon_21*epsilon_12/epsilon_11) (epsilon_23 - epsilon_21*epsilon_13/epsilon_11);
                epsilon_31/epsilon_11 (epsilon_32 - epsilon_31*epsilon_12/epsilon_11) (epsilon_33 - epsilon_31*epsilon_13/epsilon_11)
               ]
    end

    function convert_BC_sbpsm(BC::Union{String,Int64,Float64,ComplexF64},direction::String)
        if lowercase(BC) == "pec"
            return "PEC"
        elseif lowercase(BC) == "pmc"
            return "PMC"
        elseif lowercase(BC) == "pecpmc"
            return "PECPMC"
        elseif lowercase(BC) == "pmcpec"
            return "PMCPEC"
        elseif lowercase(BC) == "periodic"
            return "periodic"
        elseif lowercase(BC) == "bloch"
            return "Bloch"
        else
            throw(ArgumentError("Input argument $(direction)BC = \"$(BC)\" is not a supported option."))         
        end
    end

    function pick_epsilon_3d(epsilon_Eo::Array{ComplexF64}, epsilon_Ex::Array{ComplexF64}, epsilon_Ey::Array{ComplexF64}, epsilon_Ez::Array{ComplexF64}, xBC::String, yBC::String, zBC::String)
        # Pick the right site of components of epsilon based on Yee cell and boundary condition
    
        epsilon_xx = epsilon_Ex[:,:,:,1,1]; epsilon_xy = epsilon_Eo[:,:,:,1,2]; epsilon_xz = epsilon_Eo[:,:,:,1,3]
        epsilon_yx = epsilon_Eo[:,:,:,2,1]; epsilon_yy = epsilon_Ey[:,:,:,2,2]; epsilon_yz = epsilon_Eo[:,:,:,2,3]
        epsilon_zx = epsilon_Eo[:,:,:,3,1]; epsilon_zy = epsilon_Eo[:,:,:,3,2]; epsilon_zz = epsilon_Ez[:,:,:,3,3]
    
        if xBC == "PEC"
            epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[2:end,:,:]; epsilon_xz = epsilon_xz[2:end,:,:]
            epsilon_yx = epsilon_yx[2:end,:,:]; epsilon_yy = epsilon_yy[2:end,:,:]; epsilon_yz = epsilon_yz[2:end,:,:]        
            epsilon_zx = epsilon_zx[2:end,:,:]; epsilon_zy = epsilon_zy[2:end,:,:]; epsilon_zz = epsilon_zz[2:end,:,:]
        elseif xBC == "PMC"
            #epsilon_xx = epsilon_xx[1:end-1,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,:]       
            epsilon_xx = epsilon_xx[1:end-1,:,:];         
        elseif xBC == "PECPMC"
            epsilon_xx = epsilon_xx[1,end-1,:]; epsilon_xy = epsilon_xy[2:end,:,:]; epsilon_xz = epsilon_xz[2:end,:,:]
            epsilon_yx = epsilon_yx[2:end,:,:]; epsilon_yy = epsilon_yy[2:end,:,:]; epsilon_yz = epsilon_yz[2:end,:,:]        
            epsilon_zx = epsilon_zx[2:end,:,:]; epsilon_zy = epsilon_zy[2:end,:,:]; epsilon_zz = epsilon_zz[2:end,:,:]
        else
            #epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,:]
        end
    
        if yBC == "PEC"
            epsilon_xx = epsilon_xx[:,2:end,:]; epsilon_xy = epsilon_xy[:,2:end,:]; epsilon_xz = epsilon_xz[:,2:end,:]
            epsilon_yx = epsilon_yx[:,2:end,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,2:end,:] 
            epsilon_zx = epsilon_zx[:,2:end,:]; epsilon_zy = epsilon_zy[:,2:end,:]; epsilon_zz = epsilon_zz[:,2:end,:]
        elseif yBC == "PMC"
            #epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,1:end-1,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,:] 
            epsilon_yy = epsilon_yy[:,1:end-1,:];
        elseif yBC == "PECPMC"
            epsilon_xx = epsilon_xx[:,2:end,:]; epsilon_xy = epsilon_xy[:,2:end,:]; epsilon_xz = epsilon_xz[:,2:end,:]
            epsilon_yx = epsilon_yx[:,2:end,:]; epsilon_yy = epsilon_yy[:,1:end-1,:]; epsilon_yz = epsilon_yz[:,2:end,:] 
            epsilon_zx = epsilon_zx[:,2:end,:]; epsilon_zy = epsilon_zy[:,2:end,:]; epsilon_zz = epsilon_zz[:,2:end,:]
        else
            #epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,:]
        end    

        if zBC == "PEC"
            epsilon_xx = epsilon_xx[:,:,2:end]; epsilon_xy = epsilon_xy[:,:,2:end]; epsilon_xz = epsilon_xz[:,:,2:end]
            epsilon_yx = epsilon_yx[:,:,2:end]; epsilon_yy = epsilon_yy[:,:,2:end]; epsilon_yz = epsilon_yz[:,:,2:end] 
            epsilon_zx = epsilon_zx[:,:,2:end]; epsilon_zy = epsilon_zy[:,:,2:end]; epsilon_zz = epsilon_zz[:,:,:]
        elseif zBC == "PMC"
            #epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,1:end-1] 
            epsilon_zz = epsilon_zz[:,:,1:end-1];
        elseif zBC == "PECPMC"
            epsilon_xx = epsilon_xx[:,:,2:end]; epsilon_xy = epsilon_xy[:,:,2:end]; epsilon_xz = epsilon_xz[:,:,2:end]
            epsilon_yx = epsilon_yx[:,:,2:end]; epsilon_yy = epsilon_yy[:,:,2:end]; epsilon_yz = epsilon_yz[:,:,2:end] 
            epsilon_zx = epsilon_zx[:,:,2:end]; epsilon_zy = epsilon_zy[:,:,2:end]; epsilon_zz = epsilon_zz[:,:,1:end-1]
        else
            #epsilon_xx = epsilon_xx[:,:,:]; epsilon_xy = epsilon_xy[:,:,:]; epsilon_xz = epsilon_xz[:,:,:]
            #epsilon_yx = epsilon_yx[:,:,:]; epsilon_yy = epsilon_yy[:,:,:]; epsilon_yz = epsilon_yz[:,:,:]
            #epsilon_zx = epsilon_zx[:,:,:]; epsilon_zy = epsilon_zy[:,:,:]; epsilon_zz = epsilon_zz[:,:,:]
        end        
        return epsilon_xx, epsilon_xy, epsilon_xz, epsilon_yx , epsilon_yy, epsilon_yz, epsilon_zx, epsilon_zy, epsilon_zz 
    end


    function pick_epsilon_2d_TM(epsilon_xx::Union{Array{Int64},Array{Float64},Array{ComplexF64}}, yBC::String, zBC::String)        
        # Pick the right site of components of epsilon based on Yee cell and boundary condition
        
        if yBC == "PEC"
            epsilon_xx = epsilon_xx[2:end,:]
        elseif yBC == "PMC"
            #epsilon_xx = epsilon_xx[:,:]
        elseif yBC == "PECPMC"
            epsilon_xx = epsilon_xx[2:end,:]
        else
            #epsilon_xx = epsilon_xx[:,:]
        end    

        if zBC == "PEC"
            epsilon_xx = epsilon_xx[:,2:end]
        elseif zBC == "PMC"
            #epsilon_xx = epsilon_xx[:,:]
        elseif zBC == "PECPMC"
            epsilon_xx = epsilon_xx[:,2:end]
        else
            #epsilon_xx = epsilon_xx[:,:]
        end        
        return epsilon_xx
    end


    function pick_inv_epsilon_2d_TE(inv_epsilon_Ey_site::Union{Array{Int64},Array{Float64},Array{ComplexF64}}, inv_epsilon_Ez_site::Union{Array{Int64},Array{Float64},Array{ComplexF64}}, inv_epsilon_Eo_site::Union{Array{Int64},Array{Float64},Array{ComplexF64}}, yBC::String, zBC::String)
        # Pick the right site of components of epsilon based on Yee cell and boundary condition

        inv_epsilon_yy = inv_epsilon_Ey_site[:,:,2,2]
        inv_epsilon_zz = inv_epsilon_Ez_site[:,:,3,3]
        inv_epsilon_yz = inv_epsilon_Eo_site[:,:,2,3]
    
        if yBC == "PEC"
            inv_epsilon_yy = inv_epsilon_yy[:,:]; 
            inv_epsilon_yz = inv_epsilon_yz[2:end,:] 
            inv_epsilon_zz = inv_epsilon_zz[2:end,:]
        elseif yBC == "PMC"
            inv_epsilon_yy = inv_epsilon_yy[1:end-1,:]
            inv_epsilon_yz = inv_epsilon_yz[:,:] 
            inv_epsilon_zz = inv_epsilon_zz[:,:]
        elseif yBC == "PECPMC"
            inv_epsilon_yy = inv_epsilon_yy[1:end-1,:] 
            inv_epsilon_yz = inv_epsilon_yz[2:end,:] 
            inv_epsilon_zz = inv_epsilon_zz[2:end,:]
        else
            #inv_epsilon_yy = inv_epsilon_yy[:,:]
            #inv_epsilon_yz = inv_epsilon_yz[:,:] 
            #inv_epsilon_zz = inv_epsilon_zz[:,:]
        end    

        if zBC == "PEC"
            inv_epsilon_yy = inv_epsilon_yy[:,2:end]
            inv_epsilon_yz = inv_epsilon_yz[:,2:end] 
            inv_epsilon_zz = inv_epsilon_zz[:,:]
        elseif zBC == "PMC"
            inv_epsilon_yy = inv_epsilon_yy[:,:] 
            inv_epsilon_yz = inv_epsilon_yz[:,:]        
            inv_epsilon_zz = inv_epsilon_zz[:,1:end-1]
        elseif zBC == "PECPMC"
            inv_epsilon_yy = inv_epsilon_yy[:,2:end] 
            inv_epsilon_yz = inv_epsilon_yz[:,2:end] 
            inv_epsilon_zz = inv_epsilon_zz[:,1:end-1]
        else
            #inv_epsilon_yy = inv_epsilon_yy[:,:]
            #inv_epsilon_yz = inv_epsilon_yz[:,:] 
            #inv_epsilon_zz = inv_epsilon_zz[:,:]
        end        
        return inv_epsilon_yy, inv_epsilon_zz, inv_epsilon_yz
    end

    function Kottke_smoothing(vol_frac::Float64, n0::SVector{3, Float64}, epsilon_object::Union{Array{<:Int},Array{<:Real},Array{<:Complex}}, epsilon_voxel::Union{Array{<:Real},Array{<:Complex}})
        Scomp = @SMatrix rand(3,3-1)  # directions complementary to n12; works even for K = 1
        Stemp = [n0 Scomp]  # SMatrix{K,K}
        S = qr(Stemp).Q  # nonallocating; 1st column is normalized n12

        tau_1 = tau_trans(transpose(S) * epsilon_voxel * S)  # express param1 in S coordinates, and apply tau transform
        tau_2 = tau_trans(transpose(S) * epsilon_object * S)  # express param2 in S coordinates, and apply tau transform
        tau_avg = tau_2 .* vol_frac + tau_1 .* (1-vol_frac)  # volume-weighted average             
        epsilon = S * tau_inverse_trans(tau_avg) * transpose(S)
        return epsilon
    end